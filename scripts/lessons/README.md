# 教材の動画を生成する

各章の原稿には、発話と字幕の `text`、必要に応じて発音を指定する `reading`、画面に出す数式 `latex` と木 `tree` を収めています。基礎編の数式と木は `scripts/lessons/math_slides.json` から読み込みます。

生成には VOICEVOX Engine 0.26.2 のずんだもん（ノーマル、話者 ID 3）、Python 3.10 以降、Pillow、imageio-ffmpeg、PyAV、日本語を含む通常・太字フォントを使います。数式の組版には pdfLaTeX（amsmath と amssymb）と Ghostscript を使います。以下はリポジトリ直下で、エンジンがポート 50021 で待ち受けている場合の例です。

```sh
export VOICEVOX_URL=http://127.0.0.1:50021
export LESSON_FONT=/path/to/BIZUDGothic-Regular.ttf
export LESSON_FONT_BOLD=/path/to/BIZUDGothic-Bold.ttf
python scripts/lessons/render.py --visuals scripts/lessons/math_slides.json
python scripts/lessons/verify_video.py
python scripts/lessons/join.py --exclude 10 11 12
```

基礎編は `build/lessons/complete.mp4` と同名の SRT に出ます。個別動画は第00～13章のすべてを生成します。結合版には第00～09章と第13章を収め、文法の詳説には独立した次の二本を使います。

```sh
for course in nonterminals pk-workshop
do
  script="documents/lessons/${course}-video.json"
  output="build/lessons-${course}"
  python scripts/lessons/render.py --script "$script" --output "$output"
  python scripts/lessons/verify_video.py --script "$script" --output "$output"
  python scripts/lessons/join.py --script "$script" --output "$output"
  cp "$output/complete.mp4" "documents/lessons/media/${course}.mp4"
  cp "$output/complete.srt" "documents/lessons/media/${course}.srt"
done
```

一章だけ再描画するには render.py に `--chapter 03` のように番号を渡します。結合処理は章ごとの映像フレーム数と音声サンプル数を調べ、章選択と字幕の時刻を揃えます。音声設定やエンジンを変更した場合は新しい出力先・音声キャッシュを指定して再合成します。[配布用の動画と字幕](../../documents/lessons/media/)
