"""Join rendered chapters on a common clock, with Japanese subtitles and bookmarks."""
from fractions import Fraction
from pathlib import Path
import json
import subprocess
import tempfile
import wave

import av
import imageio_ffmpeg

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "build/lessons"
RATE = 24000
FFMPEG = imageio_ffmpeg.get_ffmpeg_exe()


def ffmpeg(*args):
    subprocess.run([FFMPEG, "-hide_banner", "-loglevel", "error", "-y", *args], check=True)


def timestamp(ms):
    seconds, ms = divmod(ms, 1000)
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours:02}:{minutes:02}:{seconds:02},{ms:03}"


def build(folder, chapters):
    positions = []
    subtitles = []
    inputs = ["ffconcat version 1.0"]
    elapsed = 0
    with wave.open(str(folder / "sound.wav"), "wb") as sound:
        sound.setparams((1, 2, RATE, 0, "NONE", ""))
        for chapter in chapters:
            name = "chapter-" + chapter["id"]
            source = OUT / (name + ".mp4")
            clip = folder / (name + ".mp4")
            ffmpeg("-i", str(source), "-map", "0:v:0", "-c", "copy", str(clip))
            with av.open(str(clip)) as video:
                stream = video.streams.video[0]
                assert stream.average_rate == 5
                video_samples = stream.frames * (RATE // 5)
            with wave.open(str(OUT / name / "audio.wav"), "rb") as audio:
                assert (audio.getframerate(), audio.getsampwidth(), audio.getnchannels()) == (RATE, 2, 1)
                samples = audio.getnframes()
                sound.writeframesraw(audio.readframes(samples))
            # Keep every speech sample and let the last still span any rounding gap.
            length = max(video_samples, samples)
            sound.writeframesraw(bytes(2 * (length - samples)))
            start = round(1000 * elapsed / RATE)
            end = round(1000 * (elapsed + length) / RATE)
            positions.append((start, end, chapter["title"]))
            inputs.extend([f"file '{clip.name}'", f"duration {length / RATE:.9f}"])
            cues = json.loads((OUT / name / "cues.json").read_text())
            for cue in cues:
                begin = round(1000 * (elapsed / RATE + cue["start"]))
                finish = round(1000 * (elapsed / RATE + cue["end"]))
                subtitles.append((begin, finish, cue["text"]))
            elapsed += length
    (folder / "clips.txt").write_text("\n".join(inputs) + "\n")
    srt = folder / "complete.srt"
    srt.write_text("\n".join(f"{i}\n{timestamp(s)} --> {timestamp(e)}\n{t}\n"
                              for i, (s, e, t) in enumerate(subtitles, 1)))
    target = folder / "complete.mp4"
    ffmpeg("-f", "concat", "-safe", "0", "-i", str(folder / "clips.txt"),
           "-i", str(folder / "sound.wav"), "-i", str(srt), "-map", "0:v", "-map", "1:a", "-map", "2:s",
           "-map_chapters", "-1", "-video_track_timescale", "10000", "-c:v", "copy", "-c:a", "aac", "-b:a", "96k",
           "-c:s", "mov_text", "-enc_time_base:s", "1:1000", "-metadata:s:s:0", "language=jpn", "-metadata:s:a:0", "language=jpn",
           str(folder / "tracks.mp4"))
    with av.open(str(folder / "tracks.mp4")) as tracks, av.open(
            str(target), "w", options={"movflags": "+faststart"}) as bookmarked:
        mapping = {s.index: bookmarked.add_stream_from_template(s) for s in tracks.streams}
        for source in tracks.streams:
            mapping[source.index].metadata.update(source.metadata)
        bookmarked.set_chapters([
            dict(id=i, time_base=Fraction(1, 1000), start=start, end=end,
                 metadata={"title": c["id"] + " " + title})
            for i, ((start, end, title), c) in enumerate(zip(positions, chapters))
        ])
        for packet in tracks.demux():
            if packet.dts is not None:
                packet.stream = mapping[packet.stream.index]
                bookmarked.mux(packet)
    with av.open(str(target)) as movie:
        actual = []
        for packet in movie.demux(subtitles=0):
            data = bytes(packet)
            size = int.from_bytes(data[:2], "big")
            if size:
                actual.append((round(1000 * packet.pts * packet.time_base),
                               data[2:2+size].decode()))
        assert actual == [(s,t) for s,_,t in subtitles]
        bookmarks = movie.chapters()
        assert len(bookmarks) == len(positions)
        for entry, (start, end, title), chapter in zip(bookmarks, positions, chapters):
            assert round(entry["start"] * entry["time_base"] * 1000) == start
            assert abs(entry["end"] * entry["time_base"] - end / 1000) < 0.2
            assert entry["metadata"]["title"] == chapter["id"] + " " + title
    with av.open(str(target)) as movie:
        duration = sum(f.samples / f.sample_rate for f in movie.decode(audio=0))
        assert abs(duration - elapsed / RATE) < 0.1
    with av.open(str(target)) as movie:
        last = None
        for frame in movie.decode(video=0):
            last = frame
        assert last is not None and abs(last.time + 0.2 - elapsed / RATE) < 0.35
    target.replace(OUT / target.name)
    srt.replace(OUT / srt.name)
    (OUT / "complete-chapters.json").write_text(json.dumps(positions, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    plan = json.loads((ROOT / "documents/lessons/video.json").read_text())
    with tempfile.TemporaryDirectory(dir=OUT) as work:
        build(Path(work), plan)
