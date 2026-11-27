"""Draw an RNA along its 5′–3′ axis using the actual pairs of a lesson witness."""
from functools import lru_cache
from pathlib import Path
import math
import re

from PIL import Image, ImageDraw

from motion import record
from typeset import formulas, layout_tree


WITNESSES = Path(__file__).resolve().parents[2] / "documents/lessons/derivations"
FIXED = "#6ecde9"
ADDED = "#e3a8ff"
REGIONS = ("#d3adff", "#ffe28e", "#97e2a9")


def pairs(brackets, opening, closing):
    stack, result = [], []
    for i, sign in enumerate(brackets, 1):
        if sign == opening:
            stack.append(i)
        elif sign == closing:
            result.append((stack.pop(), i))
        elif sign != ".":
            raise ValueError(f"Unknown bracket {sign!r}")
    if stack:
        raise ValueError("Unbalanced structure")
    return result


@lru_cache(maxsize=24)
def witness(name):
    source = (WITNESSES / (name + ".md")).read_text()
    g = re.search(r"^固定対　([^\n]+)$", source, re.M)[1]
    extra = re.search(r"^追加対　([^\n]+)$", source, re.M)[1]
    assert len(g) == len(extra)
    return len(g), pairs(g, "(", ")"), pairs(extra, "[", "]")


def draw_arc(canvas, spec, font, box):
    draw = ImageDraw.Draw(canvas)
    left, top, right, bottom = box
    if "case" in spec:
        length, fixed, extra = witness(spec["case"])
    else:
        length, fixed, extra = spec["length"], spec["fixed"], spec["added"]
    focus = {tuple(a) for a in spec.get("focus", [])}
    if spec.get("fixed_only"):
        extra = []
    a, b = spec.get("view", [1, length])
    def x(pos):
        return round(left + 22 + (right-left-44)*(pos-a)/(b-a))
    bottom_arc = bottom - 67
    upper_arc = top + 60
    y = (upper_arc + bottom_arc)//2 if extra else bottom_arc
    draw.text((left,top),spec["heading"],font=font(21),fill="#f7f7eb")
    draw.text((left,top+32),"━ 固定対 G",font=font(18),fill=FIXED)
    draw.text((left+145,top+32),"━ 追加対 G′",font=font(18),fill=ADDED)
    draw.text((right-128,top+32),"5′ → 3′",font=font(18),fill="#cad9dc")
    visible = [(i,j,col) for group,col in ((fixed,FIXED),(extra,ADDED))
               for i,j in group if a <= i < j <= b]
    longest = max((j-i for i,j,_ in visible),default=1)
    height = y-upper_arc
    assert height>20
    for i,j,color in sorted(visible, key=lambda p:p[0]-p[1]):
        active = not focus or (i,j) in focus
        ink = color if active else ("#47788a" if color==FIXED else "#9667a9")
        h = height*math.sqrt((j-i)/longest)
        points = []
        for t in range(81):
            angle = math.pi*t/80
            points.append(((x(i)+x(j))/2-(x(j)-x(i))/2*math.cos(angle),
                           y+(-1 if color==FIXED else 1)*h*math.sin(angle)))
        draw.line(points,fill=ink,width=4 if active else 2)
        if active:
            record(canvas,points,"arc",str(i),str(j),ink)
    draw.line((x(a),y,x(b),y),fill="#bed0d1",width=2)
    for pos in range(a,b+1):
        xx=x(pos)
        draw.ellipse((xx-2,y-2,xx+2,y+2),fill="#e5f3ee")
    ticks = sorted(set(spec.get("ticks",[a,b])))
    for j,pos in enumerate(ticks):
        # Alternate vertical positions when adjacent residues need both labels.
        offset = 13 if j and x(pos)-x(ticks[j-1]) < 22 else 0
        draw.text((x(pos),bottom_arc+3+offset),str(pos),font=font(15),anchor="ma",fill="#f1f6ef")
    for i,region in enumerate(spec.get("regions",[])):
        tag,start,end = region
        col=REGIONS[i%len(REGIONS)]
        record(canvas,[(x(start)-3,bottom_arc+33),(x(end)+3,bottom_arc+33)],"region",color=col)
        draw.line((x(start)-3,bottom_arc+33,x(end)+3,bottom_arc+33),fill=col,width=5)
        draw.text((left+i*((right-left)/max(1,len(spec["regions"]))),bottom-22),
                  tag,font=font(16),fill=col)


def structure_panel(canvas, slide, font, subtitle):
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((40,177,1240,512),fill="#000000")
    draw_arc(canvas,slide["arc"],font,(56,183,769,433))
    draw.line((794,188,794,430),fill="#34535b",width=2)
    if slide.get("original_figure") and subtitle in [v["text"] for v in slide["sentences"][:2]]:
        source = WITNESSES.parent / "textbook-arcs" / slide["original_figure"]
        with Image.open(source) as fig:
            fig.thumbnail((410,218),Image.Resampling.LANCZOS)
            canvas.paste(fig,(807+(410-fig.width)//2,200+(218-fig.height)//2))
        draw.text((815,179),"教科書の原図",font=font(16),fill="#c2d5ce")
    elif slide.get("tree"):
        layout_tree(canvas,slide["tree"],(802,183,1225,433),font)
    elif slide.get("definition"):
        for i,line in enumerate(slide["definition"]):
            draw.text((818,215+i*41),line,font=font(23),fill="#e1f3e9")
    formulas(canvas,slide["latex"],(55,446,1222,503))
