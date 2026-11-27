"""Base-course arc panels, using the textbook's 20/32/50-nt witnesses."""
from math import pi, sin

from PIL import ImageDraw

from motion import record
from typeset import formulas

BG = "#000000"
TEXT = "#f3f9fa"
FIXED = "#6ecde9"
ADDED = "#e3a8ff"
DIM = "#677f86"


def base_arc_panel(canvas, slide, font):
    spec = slide["base_arc"]
    d = ImageDraw.Draw(canvas)
    d.rectangle((40, 177, 1240, 512), fill=BG)
    d.text((60, 185), spec["headline"], fill=TEXT, font=font(29))
    d.text((61, 230), "● 固定対 G", fill=FIXED, font=font(20))
    d.text((260, 230), "● 追加対 G′", fill=ADDED, font=font(20))
    d.text((473, 230), "横軸：5′ → 3′", fill=TEXT, font=font(20))
    left, right, axis = 90, 777, 390
    length = spec["length"]
    gap = (right-left)/(length-1)
    def x(i):
        return left+(i-1)*gap
    fixed = spec["fixed"]
    added = spec["added"]
    all_pairs = fixed+added
    for index, (i,j) in enumerate(all_pairs):
        h = 105*(j-i)/(spec.get("max_span") or max(b-a for a,b in all_pairs))
        points = [(x(i)+(x(j)-x(i))*k/90, axis-h*sin(pi*k/90)) for k in range(91)]
        color = FIXED if index<len(fixed) else ADDED
        if "focus" in spec and [i,j] not in spec["focus"]:
            color = "#3e647a" if index<len(fixed) else "#80648a"
        d.line(points, fill=color,width=4)
        if "focus" not in spec or [i,j] in spec["focus"]:
            record(canvas, points, "arc", str(i),str(j), color)
    d.line((left,axis,right,axis),fill=DIM,width=2)
    for i in range(1,length+1):
        xx=x(i)
        d.ellipse((xx-3,axis-3,xx+3,axis+3),fill=TEXT)
    for i in spec["ticks"]:
        d.text((x(i),axis+10),str(i),fill=TEXT,font=font(19),anchor="mt")
    for i,j,color in spec.get("bars",[]):
        record(canvas,[(x(i)-gap*.33,axis+42),(x(j)+gap*.33,axis+42)], "region",color=color)
        d.rounded_rectangle((x(i)-gap*.33,axis+37,x(j)+gap*.33,axis+47),radius=3,fill=color)
    d.text((61,455),spec["caption"],fill=TEXT,font=font(22))
    d.line((814,228,814,496),fill=DIM,width=2)
    formulas(canvas,slide["latex"],(833,247,1214,470))
