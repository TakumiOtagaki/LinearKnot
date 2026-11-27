"""Compose a slide's mathematical panel from real pdfLaTeX output and a labelled tree."""
import functools
import hashlib
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from PIL import Image, ImageDraw, ImageOps, ImageFont
import re
import os

from motion import record

CACHE = Path(__file__).resolve().parents[2] / "build/lessons/latex-cache"
PREAMBLE = r"""\documentclass{article}
\usepackage{amsmath,amssymb}
\pagestyle{empty}
\begin{document}
\setbox0=\hbox{\fontsize{22}{27}\selectfont$\displaystyle FORMULA$}
\pdfpagewidth=\dimexpr\wd0+4pt\relax
\pdfpageheight=\dimexpr\ht0+\dp0+4pt\relax
\pdfhorigin=2pt\pdfvorigin=2pt
\shipout\vbox{\offinterlineskip\box0}
\end{document}
"""


@functools.lru_cache(maxsize=512)
def equation(formula):
    source = PREAMBLE.replace("FORMULA", formula)
    CACHE.mkdir(parents=True, exist_ok=True)
    target = CACHE / (hashlib.sha256(source.encode()).hexdigest() + ".png")
    if not target.exists():
        with TemporaryDirectory(dir=CACHE) as tmp:
            folder = Path(tmp)
            (folder/"formula.tex").write_text(source)
            proc = subprocess.run(["pdflatex","-no-shell-escape","-interaction=batchmode",
                                   "-halt-on-error","formula.tex"],cwd=folder,capture_output=True)
            if proc.returncode:
                raise ValueError(formula + "\n" + (folder/"formula.log").read_text()[-1600:])
            subprocess.run(["gs","-q","-dSAFER","-dNOPAUSE","-dBATCH","-dTextAlphaBits=4",
                            "-sDEVICE=pnggray","-r180","-sOutputFile=formula.png","formula.pdf"],
                           cwd=folder,check=True)
            with Image.open(folder/"formula.png") as im:
                mask = ImageOps.invert(im.convert("L"))
                bbox = mask.getbbox()
                if bbox is None:
                    raise ValueError("Empty formula: " + formula)
                mask.crop(bbox).save(folder/"mask.png")
            (folder/"mask.png").replace(target)
    with Image.open(target) as im:
        return im.copy()


def formulas(canvas, lines, box):
    x,y,right,bottom = box
    masks = [equation(line) for line in lines]
    gap = 21
    available_height = bottom-y - (len(masks)-1)*gap
    scale = min((right-x)/max(m.width for m in masks),
                available_height/sum(m.height for m in masks), .89)
    sizes = [(max(1,round(m.width*scale)),max(1,round(m.height*scale))) for m in masks]
    h = sum(b for _,b in sizes) + (len(sizes)-1)*gap
    top = y + (bottom-y-h)//2
    for mask,(w,h) in zip(masks,sizes):
        ink = mask.resize((w,h),Image.Resampling.LANCZOS)
        left = x+(right-x-w)//2
        canvas.paste("#ddf9eb",(left,top,left+w,top+h),ink)
        record(canvas, [(left, top+h+4), (left+w, top+h+4)], "equation")
        canvas.info.setdefault("rows",[]).append((left,top,left+w,top+h))
        top += h+gap


@functools.lru_cache(maxsize=1)
def math_face():
    override = os.environ.get("LESSON_MATH_FONT")
    if override:
        return override
    return subprocess.check_output(["kpsewhich","latinmodern-math.otf"],text=True).strip()


@functools.lru_cache(maxsize=64)
def math_font(size):
    return ImageFont.truetype(math_face(), size)


def label_runs(label, size, japanese):
    # Latin Modern Math includes the uppercase grammar identifiers, digits,
    # Greek letters and operators. Use the Japanese font for Japanese text.
    for run in re.findall(r"[ -~α-ωΑ-Ω×÷＋−→←∂∅□₀-₉]+|[^ -~α-ωΑ-Ω×÷＋−→←∂∅□₀-₉]+",label):
        math = run[0].isascii() or run[0] in "αβγΔ×÷＋−→←∂∅□₀₁₂₃₄₅₆₇₈₉"
        yield run.replace("＋","+") if math else run, math_font(size) if math else japanese(size)


def label_width(draw, label, size, japanese):
    return sum(draw.textlength(run,font=f) for run,f in label_runs(label,size,japanese))


def layout_tree(canvas, tree, box, font):
    draw = ImageDraw.Draw(canvas)
    nodes, edges, leaves = [], [], []
    shared = {}
    def visit(item, depth):
        child_indices = [visit(c, depth+1) for c in item.get("children", [])]
        index = len(nodes)
        if child_indices:
            first,last = nodes[child_indices[0]][0],nodes[child_indices[-1]][0]
            position = (first+last)/2
        else:
            position = len(leaves)
            leaves.append(index)
        if "shared" in item:
            shared.setdefault(item["shared"], []).append(index)
        nodes.append((position,depth,item["label"]))
        edges.extend((index,child) for child in child_indices)
        return index
    visit(tree,0)
    x,y,right,bottom = box
    width = (right-x)/len(leaves)
    levels = max(depth for _,depth,_ in nodes)+1
    level_h = (bottom-y)/levels
    vertices = [(round(x+(p+.5)*width), round(y+(d+.5)*level_h)) for p,d,_ in nodes]
    # Preserve both incoming uses but place their common input at a single vertex.
    aliases = {}
    for copies in shared.values():
        anchor = copies[0]
        assert all(nodes[i][2] == nodes[anchor][2] for i in copies)
        vertices[anchor] = (round(sum(vertices[i][0] for i in copies)/len(copies)),
                            max(vertices[i][1] for i in copies))
        aliases.update({i:anchor for i in copies[1:]})
    edges = [(aliases.get(a,a), aliases.get(b,b)) for a,b in edges]
    boxes=[]
    for (cx,cy),(_,_,text) in zip(vertices,nodes):
        max_width = min(width-6, right-x-10)
        size = min(34, max(19, round(level_h*.48)))
        while label_width(draw,text,size,font)>max_width-16 and size>17:
            size-=1
        lines=[""]
        for ch in text:
            if ch=="\n" or (lines[-1] and label_width(draw,lines[-1]+ch,size,font)>max_width-12):
                lines.append("")
            if ch!="\n":
                lines[-1]+=ch
        line_height = size+4
        w = min(max_width, max(label_width(draw,t,size,font) for t in lines)+18)
        h = len(lines)*line_height+6
        if h>level_h-4:
            raise ValueError("Tree label too long for diagram: "+text)
        boxes.append((cx,cy,w,h,size,lines,line_height))
    for parent,child in edges:
        px,py,_,ph,*_ = boxes[parent]
        cx,cy,_,ch,*_ = boxes[child]
        draw.line((px,py+ph/2,cx,cy-ch/2),fill="#8ed2ad",width=4)
        record(canvas, [(cx,cy-ch/2),(px,py+ph/2)], "edge", nodes[child][2],nodes[parent][2])
    for index, (cx,cy,w,h,size,lines,line_height) in enumerate(boxes):
        if index in aliases:
            continue
        canvas.info.setdefault("nodes",{})[nodes[index][2]]=(cx-w/2,cy-h/2,cx+w/2,cy+h/2)
        draw.rounded_rectangle((cx-w/2,cy-h/2,cx+w/2,cy+h/2),radius=8,
                               fill="#111820",outline="#71938e",width=2)
        for j,line in enumerate(lines):
            left = cx-label_width(draw,line,size,font)/2
            baseline = cy+(j-(len(lines)-1)/2)*line_height+size*.34
            for run,f in label_runs(line,size,font):
                draw.text((left,baseline),run,font=f,fill="#f3fbf5",anchor="ls")
                left+=draw.textlength(run,font=f)



def math_panel(canvas, slide, font):
    draw = ImageDraw.Draw(canvas)
    draw.rectangle((40,177,1240,512),fill="#000000")
    tree = slide.get("tree")
    lines = slide.get("latex")
    if isinstance(lines,str):
        lines=[lines]
    if tree and lines:
        draw.line((655,193,655,493),fill="#34535b",width=2)
        formulas(canvas,lines,(60,195,640,495))
        layout_tree(canvas,tree,(670,190,1220,500),font)
    elif tree:
        layout_tree(canvas,tree,(60,190,1220,500),font)
    elif lines:
        formulas(canvas,lines,(65,195,1215,495))
