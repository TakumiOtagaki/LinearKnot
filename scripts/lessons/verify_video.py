"""Check the complete encoded tracks against the script and synthesis timings."""
import hashlib
import json
from pathlib import Path

import av

root = Path(__file__).resolve().parents[2]
output = root / "build/lessons"
chapters = json.loads((root / "documents/lessons/video.json").read_text())
results = []
for chapter in chapters:
    name = "chapter-" + chapter["id"]
    items = [v for s in chapter["slides"] for v in s["sentences"]]
    cues = json.loads((output / name / "cues.json").read_text())
    assert len(items) == len(cues)
    for cue, item in zip(cues, items):
        assert cue["text"] == item["text"]
        reading = item.get("reading", item["text"])
        assert cue["reading"] == reading
        assert cue["wav"] == hashlib.sha256(reading.encode()).hexdigest() + ".wav"
    media = output / (name + ".mp4")
    with av.open(media) as c:
        assert [s.type for s in c.streams] == ["video", "audio", "subtitle"]
        assert c.streams.video[0].codec_context.name == "h264"
        frames = sum(1 for frame in c.decode(video=0))
        assert frames > 0 and abs(frames / 5 - cues[-1]["end"]) < 0.35
    with av.open(media) as c:
        frames_audio = c.decode(audio=0)
        duration = sum(f.samples / f.sample_rate for f in frames_audio)
        assert abs(duration - cues[-1]["end"]) < 0.1
    with av.open(media) as c:
        subtitles = []
        for packet in c.demux(subtitles=0):
            payload = bytes(packet)
            size = int.from_bytes(payload[:2], "big")
            if size:
                subtitles.append((float(packet.pts*packet.time_base),
                                  payload[2:2+size].decode("utf-8")))
        assert len(subtitles) == len(cues)
        for (start,text), cue in zip(subtitles,cues):
            assert abs(start - cue["start"]) < 0.002 and text == cue["text"]
    results.append(dict(chapter=chapter["id"], seconds=round(duration,2), sentences=len(cues),
                        video_frames=frames, sha256=hashlib.sha256(media.read_bytes()).hexdigest()))
print(json.dumps(results,indent=2))
