# 機内で読む LinearKnot

**[動画①＋backprop 補講を冊子印刷する（A4横・面付け済み）](flight-72-booklet-A4.pdf)** ｜ [同じ72ページを順に読む（A5）](flight-72.pdf) ｜ [backprop 補講だけ読む（12ページ）](outside-backprop.pdf)

A4横で両面・短辺とじ、1面につきPDFの1ページ、倍率100%で印刷する。面付け済みPDFは36面で、紙は18枚。先頭の表は「72｜1」、裏は「2｜71」。最初の1枚で向きを確認してから印刷し、重ねて中央で折る。プリンター自身の「冊子」機能を使う場合は、通常順のA5 PDFを選び、左とじで面付けを任せる。

補講（冊子 pp. 1–12）は、同じ計算を順方向・逆方向に手で解くところから始める。実数の inside を微分する場合と log inside を微分する場合を区別し、BE の2本の腕、塩基対確率、有限 beam の順に適用する。p.12 に解答付きの総合演習がある。図版編（pp. 13–72）は動画①の全60画面を1ページずつ収録し、各ページに再生開始時刻とメモ欄を付けた。動画と同じ色の役割を保ち、背景は紙用の白とした。

| 学びたいこと | 最初に開く紙面 | 動画①の再生時刻 |
|---|---:|---:|
| 途中変数を定義して逆伝播する | 2 | 03:10 |
| 木の穴と outside の関係 | 5 | 10:09 |
| BE の腕と、兄弟に残る部分 | 9 | 11:41 |
| 出現確率と平均出現回数の違い | 7 | 24:33 |
| 塩基対に印を付けて確率を取り出す | 10 | 27:05 |

**オフライン用動画**　出発前にファイルを端末に保存して再生を確認する。GitHubのファイルページでは「Download raw file」を選ぶ。[① 基礎編](../media/complete.mp4)・[② 非終端記号編](../media/nonterminals.mp4)・[③ 演習編](../media/pk-workshop.mp4)

原理の文献は Jason Eisner (2016), *Inside-Outside and Forward-Backward Algorithms Are Just Backprop*。本補講は LinearKnot の数値例を使った日本語の解説である。[原論文（ACL Anthology）](https://aclanthology.org/W16-5901/)

元の教科書は姉妹リポジトリ `LinearCPartyPaper/textbook` の2026年10月16日版 `vol1.pdf`～`vol4.pdf`（順に98、68、43、78ページ）。実装の基準は `e2cecbb4`。単巻 `main.pdf` は10月10日の旧版。教科書の log 表と動画の実数表の換算は補講 p.8 にまとめた。補講と図版編は教材ブランチのアニメーション版動画①（33分50秒）に対応する。

再生成には動画①の各章の描画結果を含む `build/motion-edition/complete/` を用意し、リポジトリ直下で `python3 scripts/lessons/flight_print.py` を実行する。
