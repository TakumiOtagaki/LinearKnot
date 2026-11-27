"""Render narrated PNG slides, synchronized subtitles and MP4 chapters.

Requires Python 3.10+, Pillow, imageio-ffmpeg and a Japanese TrueType font.
Set LESSON_FONT and LESSON_FONT_BOLD to point to the font files. Start the
official VOICEVOX Engine (tested with the 0.26.2 Linux CPU archive) with
    run --host 127.0.0.1 --port 58123 --cpu_num_threads 4
Then execute
    python scripts/lessons/render.py
An alternate engine address can be supplied in VOICEVOX_URL. The voice
is speaker 3 (ずんだもん・ノーマル). After changing the engine version or
synthesis settings, use a fresh --output directory to avoid old audio.
The verification utility additionally requires PyAV.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import urllib.parse
import urllib.request
import wave

from PIL import Image, ImageDraw, ImageFont
import imageio_ffmpeg
from typeset import math_panel
from structure_arcs import structure_panel
from base_arcs import base_arc_panel

ROOT = Path(__file__).resolve().parents[2]
FONT = os.environ.get("LESSON_FONT", str(Path.home()/".fonts/bizud/BIZUDGothic-Regular.ttf"))
BOLD = os.environ.get("LESSON_FONT_BOLD", str(Path.home()/".fonts/bizud/BIZUDGothic-Bold.ttf"))
ENGINE = os.environ.get("VOICEVOX_URL", "http://127.0.0.1:58123")
CREDIT = "VOICEVOX:ずんだもん"
FFMPEG = imageio_ffmpeg.get_ffmpeg_exe()


def post(path, params, body=b""):
    req = urllib.request.Request(ENGINE + "/" + path + "?" + urllib.parse.urlencode(params),
                                 data=body, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as response:
        return response.read()


def speech(reading, cache):
    key = hashlib.sha256(reading.encode("utf-8")).hexdigest()
    path = cache / (key + ".wav")
    if not path.exists():
        query = json.loads(post("audio_query", {"text": reading, "speaker": 3}))
        query["speedScale"] = 1.08
        query["prePhonemeLength"] = 0.12
        query["postPhonemeLength"] = 0.32
        audio = post("synthesis", {"speaker": 3}, json.dumps(query).encode())
        with wave.open(io.BytesIO(audio)) as w:
            assert w.getnchannels() == 1 and w.getsampwidth() == 2
        path.write_bytes(audio)
        (cache / (key + ".json")).write_text(json.dumps(query, ensure_ascii=False))
    return path


def wrapped(draw, text, font, width):
    lines, line = [], ""
    for c in text:
        if c == "\n" or (line and draw.textlength(line + c, font=font) > width):
            lines.append(line)
            line = ""
        if c != "\n":
            line += c
    if line:
        lines.append(line)
    return lines


def frame(chapter, idx, slide, subtitle, path):
    canvas = Image.new("RGB", (1280, 720), "#0f2026")
    d = ImageDraw.Draw(canvas)
    def label(x, y, string, size, color="#eaf1ed", bold=False):
        d.text((x, y), string, font=ImageFont.truetype(BOLD if bold else FONT, size), fill=color)
    label(58, 30, "LinearKnot    /    第" + chapter["id"] + "章    " + chapter["title"], 25, "#8ee4af")
    d.line((58, 83, 1222, 83), fill="#38545c", width=2)
    label(58, 116, slide["title"], 48, bold=True)
    y = 220
    for line in slide["lines"]:
        font = ImageFont.truetype(FONT, 32)
        for chunk in wrapped(d, line, font, 1130):
            d.text((77, y), chunk, font=font, fill="#b9edd9")
            y += 50
    if y > 510:
        raise ValueError("Slide content overlaps subtitles: " + slide["title"])
    # These diagrams illustrate the coordinates explicitly named on the slide.
    if chapter["id"] == "00" and idx in (0,1,2):
        # The local secant near x=2 approaches the tangent to y=x².
        origin = (125, 475)
        def graphpos(x):
            return (origin[0]+165*x, origin[1]-13*x*x)
        d.line((125,475,672,475),fill="#8da6ad",width=2)
        d.line((125,475,125,324),fill="#8da6ad",width=2)
        points = [graphpos(k/40) for k in range(128)]
        d.line(points,fill="#8ee4af",width=4)
        if idx in (1, 2):
            slope = 4.1 if idx == 1 else 4
            endpoints = [(origin[0]+165*x, origin[1]-13*(4+slope*(x-2))) for x in (1.15,3.1)]
            d.line(endpoints, fill="#efa66a", width=3)
            label(728, 455, "比で測る傾き 4.1" if idx==1 else "一点での傾き 4",23,"#ffc178")
        for x in ((2,3) if idx==0 else (2,)):
            px,py=graphpos(x)
            d.ellipse((px-6,py-6,px+6,py+6),fill="#ffc178")
            label(px+10,py-20,f"({x}, {x*x})",21,"#ffc178")
        label(732,365,"横軸　入力の数",26)
        label(732,417,"縦軸　二乗した出力",26)
    if chapter["id"] == "00" and idx in (4,5):
        tags = ["入力 2","1 を足す","4 を掛ける","出力 12"]
        for k,tag in enumerate(tags):
            x=78+296*k
            d.rounded_rectangle((x,358,x+221,422),radius=12,fill="#29414a")
            label(x+23,374,tag,27)
            if k<3:
                label(x+236,370,"→" if idx==4 else "←",32,"#8ee4af")
        label(84,456,"値を作る順番　　左から右" if idx==4 else "変化率を戻す順番　　右から左",25,"#8ee4af")
    if chapter["id"] == "00" and idx==6:
        points=[(207,414),(583,357),(583,468),(985,414)]
        for j,k in [(0,1),(0,2),(1,3),(2,3)]:
            x,y0=points[j];x2,y2=points[k]
            d.line((x+81,y0,x2-82,y2),fill="#8ee4af",width=4)
        for (x,y0),tag in zip(points,["中間 3","上 9","下 6","合計 15"]):
            d.rounded_rectangle((x-82,y0-24,x+82,y0+25),8,fill="#29414a")
            label(x-62,y0-17,tag,26)
    if chapter["id"] == "01" and idx == 0:
        xs = [234 + 156 * k for k in range(6)]
        d.line((xs[0], 453, xs[-1], 453), fill="#718c91", width=3)
        for left, right, top, color in [(0, 4, 337, "#8ee4af"), (2, 5, 366, "#ffc178")]:
            d.arc((xs[left], top, xs[right], 906 - top), 180, 360, fill=color, width=5)
        for k, x in enumerate(xs):
            d.ellipse((x-6, 447, x+6, 459), fill="#eaf1ed")
            label(x-7, 466, str(k+1), 21)
    if chapter["id"] == "04" and idx in (0, 1, 2):
        for k in range(14):
            x = 116 + 76*k
            arm = k < 4 or k >= 10
            d.rounded_rectangle((x, 370, x+63, 426), radius=9,
                                fill="#306450" if arm else "#24373d")
            label(x + (23 if k<9 else 17), 383, str(k+1), 24)
        label(208, 443, "左腕", 24, "#8ee4af")
        label(490, 443, "穴（別の部分が担当）", 24, "#ffc178")
        label(984, 443, "右腕", 24, "#8ee4af")
    if chapter["id"] == "06" and idx in range(5):
        nodes = [(214,415), (623,355), (623,461), (1040,415)]
        for start, end, weight in [(0,1,2),(1,3,5),(0,2,3),(2,3,7)]:
            x, y0 = nodes[start]
            x2, y2 = nodes[end]
            color = "#8ee4af" if 1 in (start,end) else "#ffc178"
            d.line((x+42,y0,x2-48,y2),fill=color,width=4)
            tipx,tipy = x2-44,y2
            d.polygon([(tipx,tipy),(tipx-12,tipy-7),(tipx-12,tipy+7)],fill=color)
            label((x+x2)//2-10, (y0+y2)//2-29 if 1 in (start,end) else (y0+y2)//2+2,
                  str(weight),22,color)
        for (x,y0), tag in zip(nodes,["始点","A","B","終点"]):
            d.rounded_rectangle((x-42,y0-18,x+42,y0+22),8,fill="#29414a")
            label(x-(21 if len(tag)==2 else 8),y0-12,tag,24)
    if chapter["id"] == "07" and idx in (1,2,3,4):
        tags = [("a1 × b1", "4"), ("a1 × b2", "8"), ("a2 × b1", "4"), ("a2 × b2", "8")]
        if idx >= 2:
            tags.append(("s", "6"))
        for k,(tag, value) in enumerate(tags):
            x = 79 + k*227
            show_q = idx >= 3 and k in (1,3)
            d.rounded_rectangle((x,358,x+202,476),9,fill="#306450" if show_q else "#24373d")
            label(x+17,371,tag,24)
            label(x+17,410,value,34,"#8ee4af" if show_q else "#eaf1ed")
            if show_q:
                label(x+100,422,"q",23,"#8ee4af")
    if chapter["id"] == "10" and idx < 5:
        n, pairs = (12,[(2,11),(3,10)]) if idx < 3 else (
            (13,[(2,13),(4,11)]) if idx == 3 else (18,[(1,18),(3,8),(10,16)]))
        def pos(i):
            return round(100+(i-1)*1050/(n-1))
        d.line((pos(1),450,pos(n),450),fill="#76929c",width=3)
        for j,(left,right) in enumerate(pairs):
            height = 106 * (right-left)/(n-1)
            d.arc((pos(left),450-height,pos(right),450+height),180,360,
                  fill=("#91e8b0","#ffc178","#83c8ee")[j],width=5)
        for k in range(1,n+1):
            x=pos(k)
            d.ellipse((x-6,444,x+6,456), fill="#fffef7")
            label(x-9,462,str(k),20)
    if chapter["id"] == "09" and idx in (0,1,2):
        label(78,344,"全体の重み 31",24)
        for k in range(31):
            x=79+36*k
            d.rounded_rectangle((x,397,x+32,440),4,fill="#8ee4af" if k<30 else "#ffc178")
        label(78,457,"印あり 30",23,"#8ee4af")
        label(1000,457,"印なし 1",23,"#ffc178")
    if chapter["id"] == "13" and idx in (0,1):
        keep = idx == 0
        for k,(name,value,col) in enumerate([("a",9,"#8ee4af"),("b",1,"#ffc178"),("c",2,"#85c1eb")]):
            top=347+k*46
            label(83,top,name,25)
            length=value*82
            if keep or k!=1:
                d.rectangle((140,top+1,140+length,top+29),fill=col)
                label(157+length,top-1,str(value),22,col)
            else:
                d.rectangle((140,top+1,140+length,top+29),outline="#7a8c8f",width=2)
                label(157+length,top-1,"候補から除外",22)
    if chapter["id"] == "11" and idx in (2, 3, 4):
        runs = [
            ("BE", [(1,4),(11,14)]), ("WI", [(5,5)]),
            ("VP", [(6,10),(15,15),(20,23)]),
            ("WI", [(24,26)]), ("BE", [(16,19),(32,35)]),
            ("VP", [(27,31),(36,36),(41,44)]),
            ("BE", [(37,40),(47,50)]), ("WI", [(45,46)])]
        owners = [""] * 51
        for owner, regions in runs[:{2:3, 3:6, 4:8}[idx]]:
            for left, right in regions:
                for pos in range(left, right+1):
                    assert not owners[pos]
                    owners[pos] = owner
        colors = {"BE":"#72d7a2", "WI":"#ffe49b", "VP":"#85c1eb", "":"#e8eceb"}
        # Two lines of 25 positions preserve 1-based numbering at a legible size.
        for pos in range(1, 51):
            x = 76 + ((pos-1)%25)*45
            top = 342 + ((pos-1)//25)*50
            d.rounded_rectangle((x, top, x+41, top+39), radius=4, fill=colors[owners[pos]])
            label(x+6, top+6, str(pos),21,"#14232a")
        for k,(tag,color) in enumerate(colors.items()):
            x = 77+k*260
            d.rectangle((x,463,x+24,486),fill=color)
            label(x+34,459,tag or "未担当",21)
    if slide.get("arc"):
        structure_panel(canvas, slide, lambda size: ImageFont.truetype(FONT, size), subtitle)
    elif slide.get("base_arc"):
        base_arc_panel(canvas, slide, lambda size: ImageFont.truetype(FONT, size))
    elif slide.get("latex") or slide.get("tree"):
        math_panel(canvas, slide, lambda size: ImageFont.truetype(FONT, size))
    d.rectangle((0, 513, 1280, 678), fill="#20343c")
    font = ImageFont.truetype(FONT, 30)
    rows = wrapped(d, subtitle, font, 1160)
    if len(rows) > 4:
        raise ValueError("Subtitle exceeds four lines: " + subtitle)
    for i, row in enumerate(rows):
        d.text((59, 525 + i * 36), row, font=font, fill="#fffef7")
    label(59, 684, CREDIT, 19, "#adc1bb")
    label(1070, 684, f"{idx+1:02d} / {len(chapter['slides']):02d}", 19, "#adc1bb")
    canvas.save(path)


def timestamp(t):
    total = round(t * 1000)
    sec, ms = divmod(total, 1000)
    minutes, sec = divmod(sec, 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours:02d}:{minutes:02d}:{sec:02d},{ms:03d}"


def render(chapter, output, cache):
    folder = output / ("chapter-" + chapter["id"])
    folder.mkdir(exist_ok=True)
    cues, segments, audio_frames, transcript = [], [], [], []
    sample_count = 0
    rate = 24000
    for i, slide in enumerate(chapter["slides"]):
        for sentence in slide["sentences"]:
            text = sentence["text"]
            reading = sentence.get("reading", text)
            wav = speech(reading, cache)
            with wave.open(str(wav)) as w:
                assert w.getframerate() == rate and w.getnchannels() == 1 and w.getsampwidth() == 2
                frames = w.readframes(w.getnframes())
            n = len(frames) // 2
            start, end = sample_count / rate, (sample_count + n) / rate
            num = len(cues) + 1
            png = folder / f"frame-{num:03d}.png"
            frame(chapter, i, slide, text, png)
            segments.append(f"file '{png.name}'\nduration {(end-start):.8f}\n")
            audio_frames.append(frames)
            cues.append(f"{num}\n{timestamp(start)} --> {timestamp(end)}\n{text}\n")
            transcript.append(dict(start=start, end=end, text=text, reading=reading, wav=wav.name))
            sample_count += n
    segments.append(f"file '{png.name}'\n")
    (folder/"images.ffconcat").write_text("ffconcat version 1.0\n" + "".join(segments))
    with wave.open(str(folder/"audio.wav"), "wb") as w:
        w.setparams((1, 2, rate, 0, "NONE", "not compressed"))
        w.writeframes(b"".join(audio_frames))
    (folder/"cues.json").write_text(json.dumps(transcript, ensure_ascii=False, indent=2))
    srt = output / f"chapter-{chapter['id']}.srt"
    srt.write_text("\n".join(cues))
    target = output / f"chapter-{chapter['id']}.mp4"
    subprocess.run([FFMPEG, "-hide_banner", "-loglevel", "error", "-y", "-f", "concat", "-safe", "0",
                    "-i", str(folder/"images.ffconcat"), "-i", str(folder/"audio.wav"), "-i", str(srt),
                    "-map","0:v","-map","1:a","-map","2:s","-c:v","libx264","-threads","3",
                    "-vf","fps=5,tpad=stop_mode=clone:stop_duration=4","-tune","stillimage","-preset","veryfast","-crf","19",
                    "-pix_fmt","yuv420p","-c:a","aac","-b:a","96k","-c:s","mov_text",
                    "-metadata:s:s:0","language=jpn","-metadata","title="+chapter["title"],
                    "-movflags","+faststart","-t",str(sample_count/rate),str(target)],check=True)
    print(target, round(sample_count/rate, 1), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--script", type=Path, default=ROOT/"documents/lessons/video.json")
    parser.add_argument("--output", type=Path, default=ROOT/"build/lessons")
    parser.add_argument("--chapter")
    parser.add_argument("--visuals", type=Path, default=None)
    parser.add_argument("--speech-cache",type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    cache = args.speech_cache or args.output/"speech-cache"
    cache.mkdir(parents=True, exist_ok=True)
    visuals = json.loads(args.visuals.read_text()) if args.visuals else {}
    for chapter in json.loads(args.script.read_text()):
        for i, slide in enumerate(chapter["slides"], 1):
            slide.update(visuals.get(chapter["id"], {}).get(str(i), {}))
        if not args.chapter or chapter["id"] == args.chapter:
            render(chapter,args.output,cache)
