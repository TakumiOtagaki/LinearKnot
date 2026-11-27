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
    args = parser.parse_args()
    cache = args.output/"speech-cache"
    cache.mkdir(parents=True, exist_ok=True)
    for chapter in json.loads(args.script.read_text()):
        if not args.chapter or chapter["id"] == args.chapter:
            render(chapter,args.output,cache)
