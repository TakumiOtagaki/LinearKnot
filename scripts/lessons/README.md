# 教材の動画を生成する

各章の原稿には、発話と字幕の `text`、必要に応じて発音を指定する `reading`、画面に出す数式 `latex` と木 `tree` を収めています。基礎編の数式と木は `scripts/lessons/math_slides.json` から読み込みます。

非終端編・擬似ノット演習編には配列上の弧を定める `arc` があります。`case` は図譜の実例名で、そのページの固定対・追加対の括弧列から対を復元します。`focus` は強調する対の座標、`regions` は区間の凡例・左端・右端です。両端をともに含めて左から 1、2、…と数えます。`view` は拡大する区間で、両端が入る対が描かれます。模式例では `fixed` と `added` に対を明記します。`original_figure` を指定した画面は最初の二発話に教科書の元図、それ以降は定義または導出木を並べます。

計算グラフでは値を子から親へ、逆伝播の寄与を親から子へ順に描きます。分岐した同じ入力には各出現で同一の `shared` 名を付けると、1つの頂点につながります。発話の `motion.nodes` で変数名を順に指定し、`motion.direction` に `up` または `down` を指定すると、その発話でたどる枝を選べます。RNA の図では塩基対の端から相手の端へ弧を伸ばし、腕や隙間を示す帯を描きます。字幕のタイミングは音声サンプル数、アニメーションは毎秒10フレームで計算します。

木の英字・演算子・座標には Latin Modern Math を使います。TeX の `kpsewhich latinmodern-math.otf` で取得し、別の数学フォントを使う場合は `LESSON_MATH_FONT` に OpenType フォントのパスを指定します。

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
