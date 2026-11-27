"""Animate the exact paths and typeset expressions drawn on a lesson frame."""
from dataclasses import dataclass
from math import ceil, hypot
import re

from PIL import Image, ImageDraw

FPS = 10
AMBER = "#ffd27d"
CYAN = "#80dcff"


@dataclass
class Trace:
    points: list
    kind: str
    source: str = ""
    target: str = ""
    color: str = AMBER


def record(canvas, points, kind, source="", target="", color=AMBER):
    if len(points) >= 2:
        canvas.info.setdefault("traces", []).append(Trace(points, kind, source, target, color))


def names(text):
    lhs = re.search(r"([A-Za-z]+)\s*=", text)
    if lhs:
        return [lhs[1]]
    return re.findall(r"(?<![A-Za-z])[A-Za-z]+(?=\s*[(:：]|\s*$)", text)



def timeline(canvas, chapter_id, slide_num, sentence, cue_index):
    traces = canvas.info.get("traces", [])
    graph = [t for t in traces if t.kind == "edge"]
    pairs = [t for t in traces if t.kind == "arc"]
    regions = [t for t in traces if t.kind == "region"]
    equations = [t for t in traces if t.kind == "equation"]
    mode = sentence.get("motion", {})
    backward = {"00": {6, 9}, "03": {1, 2, 3, 4}, "05": {3},
                "06": {3, 5}, "07": {5}, "08": {2, 4, 5, 6}}
    reverse = mode.get("direction", "down" if not pairs and slide_num in backward.get(chapter_id, ()) else "up") == "down"

    if reverse:
        # Start at the root, then visit the edges at each successive depth.
        graph = sorted(graph, key=lambda t: t.points[-1][1])
        graph = [Trace(t.points[::-1], t.kind, t.target, t.source, CYAN) for t in graph]
    targets = mode.get("nodes", [])
    if targets and graph:
        # A fork x -> (a,b) -> y has two routes. Preserve all edges between
        # named vertices, including the two different incoming edges at y.
        order = {v:i for i,v in enumerate(targets)}
        selected = [t for t in graph if any(a in order and b in order and order[a]<order[b]
                           for a in names(t.source) for b in names(t.target))]
        if selected:
            graph = sorted(selected,key=lambda t: min(order[b] for b in names(t.target) if b in order))
    if chapter_id == "00" and slide_num in (5,6,7) and targets:
        forward = {"x":"x = 2","u":"u = 3","v":"v = 12"}
        tangent = {"x":"dx/dx = 1","u":"du/dx = 1","v":"dv/dx = 4"}
        adjoint = {"v":"dv/dv = 1","u":"dv/du = 4","x":"dv/dx = 4"}
        if slide_num == 7:
            forward = {"x":"x = 3","a":"a = 9","b":"b = 6","y":"y = 15"}
            adjoint = {"y":"dy/dy = 1","a":"dy/da = 1","b":"dy/db = 1","x":"dy/dx = 5"}
        canvas.info["evaluation"] = adjoint if reverse else tangent if slide_num==5 and cue_index>0 else forward
    canvas.info["reverse"] = reverse
    elements = pairs + regions if pairs or regions else graph
    if not elements:
        elements = [t for t in traces if t.kind == "line"]
    if equations:
        # Carry the result back to the equation after following the geometry.
        elements = elements + [equations[min(cue_index, len(equations)-1)]]
    return elements


def animate(still, steps, frame_number, frames, cue_index):
    if not steps:
        return still
    # Build the graph in successive stages. Later operands stay readable in the
    # context view; the active node switches from its defining expression to
    # its evaluated value (or its adjoint on a backward pass).
    passes = max(1, ceil(frames / FPS / (len(steps)*3.2)))
    position = ((frame_number + .5) / frames * passes * len(steps)) % len(steps)
    index = int(position)
    step = steps[index]
    arrived = position % 1 >= .82
    phase = (index,arrived)
    if still.info.get("phase") != phase:
        scene = still.copy()
        scene.info.clear()
        pen = ImageDraw.Draw(scene)
        for j,t in enumerate(steps):
            if t.kind in ("arc","edge","region") and j>=index:
                # Future paths remain as graph context in charcoal.
                pen.line(t.points,fill="#383a40",width=4)
            if t.kind == "equation" and j==index and len(steps)>1:
                pen.line(t.points,fill="#474032",width=2)
        nodes=still.info.get("nodes",{})
        active=steps[:index]+[step]
        reached={t.source for t in active if t.kind=="edge"}
        reached.update(t.target for t in steps[:index+int(arrived)] if t.kind=="edge")
        for label,box in nodes.items():
            involved = any(t.kind=="edge" and label in (t.source,t.target) for t in steps)
            if not involved:
                continue
            current = step.kind=="edge" and label==(step.target if arrived else step.source)
            if label not in reached:
                bounds=tuple(map(round,box))
                tile=scene.crop(bounds)
                tile=Image.blend(Image.new("RGB",tile.size),tile,.44)
                scene.paste(tile,bounds[:2])
            elif label in reached:
                color=step.color if current else "#77958a"
                expression=next((still.info.get("evaluation",{}).get(n) for n in names(label)),None)
                if expression == "dy/dx = 5":
                    incoming=[t for t in steps[:index+int(arrived)] if t.kind=="edge" and t.target==label]
                    if len(incoming)==1:
                        expression="dy/dx += " + ("3" if "a" in names(incoming[0].source) else "2")
                if expression:
                    x,y,r,b=box
                    pen.rounded_rectangle(box,8,fill="#172023",outline=color,width=3)
                    from typeset import math_font
                    size=28
                    while size>12 and pen.textlength(expression,font=math_font(size))>r-x-6:
                        size-=1
                    pen.text(((x+r)/2,(y+b)/2),expression,font=math_font(size),
                             fill="#ffffff",anchor="mm")
                elif current:
                    pen.rounded_rectangle(box,8,outline=color,width=3)
        # Equations unfold row by row in the first explanation and remain
        # available when subsequent sentences revisit the same diagram.
        rows=still.info.get("rows",[])
        if cue_index==0 and len(rows)>1:
            upto = 1+int((index+int(arrived))*(len(rows)-1)/max(1,len(steps)-1))
            for row in rows[upto:]:
                pen.rectangle(row,fill="#000000")
        still.info["phase"], still.info["scene"] = phase,scene
    canvas = still.info["scene"].copy()
    draw = ImageDraw.Draw(canvas)
    progress = position % 1
    points = step.points
    # Parametrise the polyline to keep the marker's speed continuous on arcs.
    segment = progress * (len(points)-1)
    i = min(int(segment),len(points)-2)
    x0, y0 = points[i]
    x1, y1 = points[i+1]
    fraction = segment-i
    tip = (x0+(x1-x0)*fraction, y0+(y1-y0)*fraction)
    line = [*points[:i+1], tip]
    draw.line(line, fill=step.color, width=4)
    x, y = tip
    if step.kind == "edge":
        dx,dy=x1-x0,y1-y0
        length=max(.01,hypot(dx,dy))
        dx,dy=dx/length,dy/length
        draw.polygon([(x+6*dx,y+6*dy),(x-7*dx-5*dy,y-7*dy+5*dx),
                      (x-7*dx+5*dy,y-7*dy-5*dx)], fill=step.color)
    else:
        draw.ellipse((x-4,y-4,x+4,y+4), fill=step.color)
    return canvas
