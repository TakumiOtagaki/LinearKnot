# 教材の動画を生成する

原稿は `documents/lessons/video.json` にあります。各画面の `lines` が画面に出す式と図の説明、各発話の `text` が字幕、任意の `reading` がずんだもんに読み上げさせる文字列です。日本語の字幕を校正するときは `documents/lessons/transcript.md` にも反映してください。

この教材の生成には VOICEVOX Engine 0.26.2 のずんだもん（ノーマル、話者 ID 3）を使いました。エンジンを起動した環境で、Python の Pillow と imageio-ffmpeg、日本語グリフを含む通常・太字のフォントを用意します。検証には PyAV（パッケージ名 `av`）も使います。以下はリポジトリ直下で、エンジンが標準のポート 50021 で動いている場合の例です。

```sh
export VOICEVOX_URL=http://127.0.0.1:50021
export LESSON_FONT=/path/to/BIZUDGothic-Regular.ttf
export LESSON_FONT_BOLD=/path/to/BIZUDGothic-Bold.ttf
python scripts/lessons/render.py
python scripts/lessons/verify_video.py
python scripts/lessons/join.py
```

結果は `build/lessons/` に出ます。一章だけを更新するには、例えば `render.py --chapter 11` を使います。台詞ごとの合成音声を再利用するので、別のエンジンや音声設定を試すときは `render.py --output build/lessons-new` のように出力先を替えます。配布済みの MP4 と SRT は [media フォルダー](../../documents/lessons/media/) にあります。

結合版は `build/lessons/complete.mp4`、通しの日本語字幕は `complete.srt` です。結合では各章の音声サンプル数と映像フレーム数から境界を計算し、章選択の目印と字幕時刻に同じ境界を使います。
