"""Build the 72-page flight companion from the released motion-edition render.

Run with Python containing Pillow and pypdf from the checkout root. LuaLaTeX
(luatexja, BIZ UDGothic), pdfLaTeX, and the original rendered frames are required.
The optional sibling textbook packets stay in build/flight-print for local use.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess

from PIL import Image, ImageChops, ImageOps
from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / 'documents/lessons'
PUBLIC = DOC / 'print'
BUILD = ROOT / 'build/flight-print'
FRAMES = ROOT / 'build/motion-edition/complete'
mm = 72 / 25.4


def escape(text):
    special = {'&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_',
               '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
               '\\': r'\textbackslash{}'}
    return ''.join(special.get(c, c) for c in text)


def compile_tex(source, engine):
    cache = BUILD / 'texcache'
    cache.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TEXMFCACHE=str(cache), TEXMFVAR=str(cache),
               TEXINPUTS=f'{DOC}/textbook-arcs/:{BUILD}/:')
    with (BUILD / (source.stem + '-compile.log')).open('w') as log:
        subprocess.run([engine, '-interaction=nonstopmode', '-halt-on-error',
                        '-output-directory='+str(BUILD), str(source)],
                       cwd=BUILD, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    return BUILD / (source.stem+'.pdf')


def appendix():
    preamble=(PUBLIC/'outside-backprop.tex').read_text().split(r'\newcommand{\unit}')[0]
    pages=[preamble, r'\fancyfoot[L]{\small LinearKnot　動画①・紙上版}\fancyfoot[R]{\thepage\,/\,72}',
           r'\begin{document}\setcounter{page}{13}']
    chapters=json.loads((FRAMES/'complete-chapters.json').read_text())
    selected=[c for c in json.loads((DOC/'video.json').read_text()) if c['id'] not in ('10','11','12')]
    assert len(selected)==len(chapters)
    for chapter,(begin,_,title) in zip(selected,chapters):
        assert chapter['title']==title
        folder=FRAMES/('chapter-'+chapter['id'])
        cues=json.loads((folder/'cues.json').read_text())
        n=0
        for i,slide in enumerate(chapter['slides'],1):
            sec=int(begin/1000+cues[n]['start'])
            texts=[s['text'] for s in slide['sentences']]
            assert [s['text'] for s in cues[n:n+len(texts)]]==texts
            n+=len(texts)
            image=Image.open(folder/f'frame-{n:03d}.png').convert('RGB').crop((40,178,1241,512))
            # Preserve the visual category hues while mapping the dark screen to white paper.
            h,s,v=image.convert('HSV').split()
            # Bright coloured strokes would become almost black regardless of hue.
            # Retain enough value there to keep fixed and added arcs distinguishable.
            ink_value=ImageChops.lighter(ImageOps.invert(v),s.point(lambda x: int(.48*x)))
            paper=Image.merge('HSV',(h,s,ink_value)).convert('RGB')
            pic=f'panel-{chapter["id"]}-{i:02}.png'
            paper.save(BUILD/pic)
            pages.extend([r'\noindent{\small 動画①　第'+chapter['id']+f'章　画面 {i:02}　再生 {sec//60:02}:{sec%60:02}'+r'}\par',
                          r'{\Large\bfseries '+escape(slide['title'])+r'}\par\medskip',
                          r'\includegraphics[width=\linewidth]{'+pic+r'}\par\medskip'])
            for text in texts:
                pages.append(escape(text)+r'\par')
            pages.append(r'\vfill\noindent\textcolor{gray}{\rule{\linewidth}{0.3pt}}\par{\small 計算・メモ}\par\vspace{8mm}\clearpage')
    tex=BUILD/'flight-frames.tex'
    tex.write_text('\n'.join(pages+[r'\end{document}']))
    return compile_tex(tex,'lualatex')


def combine(dest,parts):
    writer=PdfWriter()
    for title,pdf,first,last in parts:
        writer.append(str(pdf),outline_item=title,pages=(first-1,last),import_outline=False)
    writer.write(str(dest))
    writer.close()


def impose(source,dest):
    total=len(PdfReader(source).pages)
    assert total%4==0
    lines=[r'\documentclass{article}\usepackage[paperwidth=297mm,paperheight=210mm,margin=0mm]{geometry}\usepackage{graphicx}\pagestyle{empty}\begin{document}']
    order=[]
    for sheet in range(total//4):
        for left,right in [(total-2*sheet,2*sheet+1),(2*sheet+2,total-2*sheet-1)]:
            order.extend([left,right])
            def slot(n):
                return r'\makebox[148.5mm][c]{\includegraphics[page='+str(n)+r',width=148mm,height=209mm,keepaspectratio]{'+str(source)+r'}}'
            lines.append(r'\noindent\hbox to \paperwidth{'+slot(left)+slot(right)+r'\hss}\newpage')
    assert sorted(order)==list(range(1,total+1))
    tex=BUILD/(dest.stem+'.tex')
    tex.write_text('\n'.join(lines+[r'\end{document}']))
    pdf=compile_tex(tex,'pdflatex')
    assert len(PdfReader(pdf).pages)==total//2
    if pdf!=dest:shutil.copy2(pdf,dest)


if __name__=='__main__':
    BUILD.mkdir(parents=True,exist_ok=True)
    guide=compile_tex(PUBLIC/'outside-backprop.tex','lualatex')
    panels=appendix()
    assert len(PdfReader(guide).pages)==12 and len(PdfReader(panels).pages)==60
    shutil.copy2(guide,PUBLIC/guide.name)
    combine(PUBLIC/"forward-backward-p6.pdf", [("4本の矢印と経路の重み", guide, 6, 6)])
    flight=PUBLIC/'flight-72.pdf'
    combine(flight,[('Outside と backprop',guide,1,12),('動画① 全60画面',panels,1,60)])
    impose(flight,PUBLIC/'flight-72-booklet-A4.pdf')
