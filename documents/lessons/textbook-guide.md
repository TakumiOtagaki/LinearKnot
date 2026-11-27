# 教科書の弧図と実装の式を行き来する

教科書第1章の20塩基例には、同じ二次構造の弧図と、その導出で各部分木が受け持つ位置の図がある。まず [20 塩基の弧図の読み方](grammar.md#20塩基の弧図でバンドを見つける) で黒と赤の対を読み、次に[教科書の図版](textbook-arcs/README.md)で、同じ規則を別の幾何に適用したときに変わる位置を比較できる。章名と図のラベルは教材の元になった `LinearCPartyPaper/textbook` のものを用いる。

| 問い | 教科書の章・定義または図 | この教材の手計算・導出 |
|---|---|---|
| 弧はどこで交差し、どれが同じ band に属すか | 第1章「ループ、閉領域、バンド」定義 `prelim-band`、図 `prelim-ownership20` | [黒の4本・赤の3本を端点から決める](grammar.md#20塩基の弧図でバンドを見つける) |
| どの構造の重みを足すか | 第2章「条件付き ensemble と density-2」定義 `quantities-density` | [2つのバンドの density と導出の和](grammar.md) |
| 未対合の列と枝のある領域をどの子が引き受けるか | 第3章「WI」「WIP」「VPL」「VPR」、図版 `03_wi_wip_context` と `03_vp_flanks` | [15状態・57識別子の規則表](nonterminals.md)、[VP の左右の全9配置](derivations/README.md) |
| 同じ座標が2つの子の範囲に入ってもよいのはなぜか | 第4章「区間と受け持ちを分ける」図 `pkblock-f1-ownership` | [50塩基の木](derivations/chain50.md)、[1段の積を微分する](equations.md#5穴のある-be-で同じ微分を実行する) |
| 3本の固定バンドをどう順に接続するか | 第4章「最初に鎖の全体を見る」、定義 `pkblock-pending` と `pkblock-wmb` | [右端の保留と次の固定対を座標でたどる](nonterminals.md#実際の木を読む4つの問題) |
| 何が数値の inside で、その逆向き計算は何を返すか | 第5章「実装における exact inside 値」、第6章「outside 値、穴のある文脈、事後値」 | [計算木の和・積と偏微分](inside-outside.md)、[数式編](equations.md) |

## 実数の重みと教科書の log 表

本教材の $\alpha(v)$（アルファ、inside）は、項目 $v$ を根とする導出の重みを実数で足した値。$\beta(v)$（ベータ、outside）はその項目の値が全体の和 $Z$ に届く倍率である。根を $r$ と書くと、

$$
Z=\alpha(r),\qquad \beta(v)=\frac{\partial Z}{\partial\alpha(v)},\qquad \beta(r)=1.
$$

教科書第5章の擬似コードは、実数の和を $A(v)$、その自然対数を $\alpha(v)$ と記している。第4章で子の値を指数関数で実数に戻す式を読むときも同じ換算を用いる。2つを並べて読むため、対数の表の値に $\ell$ を添えよう。

$$
\ell_\alpha(v)=\log\alpha(v),\qquad \ell_\beta(v)=\log\beta(v),\qquad
\alpha(v)=e^{\ell_\alpha(v)},\quad \beta(v)=e^{\ell_\beta(v)}.
$$

例えば2つの子の実数 inside が 2 と 3、局所重みが 4 なら、この規則の寄与は $4\cdot 2\cdot 3=24$。対数表では $\log 4+\log 2+\log 3=\log 24$ を得る。別の規則の寄与が 6 なら親の実数 inside は 30、対数表に置くのは $\log(e^{\log 24}+e^{\log 6})=\log 30$。重みが 0 のときの対数表の値は $-\infty$ とする。

対数 outside を保存することと、出力の対数を入力の対数で微分することには別の式がある。$Z>0,\alpha(v)>0$ のとき、連鎖律から

$$
\frac{\partial\log Z}{\partial\ell_\alpha(v)}
=\frac{\beta(v)\alpha(v)}{Z}.
$$

右辺は導出中でのその項目の出現回数の期待値である。項目が導出あたり高々1度なら出現確率にもなる。

## 教科書の状態一覧と現在の実装

教科書第3章は、当時の実装特有の補助値 `VP_CLOSED` と `VP_DIRECT` を含む17記号を表にしている。現行版の [NonTerminal の定義](../../src/scfg2/grammar/nonterminal.hh)にはこの2名がなく、[非終端記号編](nonterminals.md)は現行版の15記号を扱う。導出図譜に示した座標付きの木には、教科書に保存された実行記録の版を付した。教科書が構造に対して与える意味と、その例で実測した導出を、各節で順に照合できる。
