# VP の左右ポケット EB

固定対を丸括弧、追加対を角括弧で示す。同じ種類の開き括弧と閉じ括弧を、内側から順に対応させて読む。入力で決まっている固定対は、結合時に再び選び直さない。

```text
配列　　GGGAAACAAACGAAACC
固定対　(.....)..........
追加対　.[[.......][...]]
完成形　([[...)...][...]]
```

![具体的な導出木](vp_EB.svg)

[図を単独で拡大](vp_EB.svg)

| # | 状態 | 選んだ規則 | 子の番号（順序付き） |
|---:|---|---|---|
| 1 | `W(1,17)` | `W_SPLIT_WMB` | 2 |
| 2 | `WMB(1,17)` | `WMB_DIRECT_WMBP` | 3 |
| 3 | `WMBP(1,17)` | `WMBP_SPLIT_BE_WI_VP` | 4, 5 |
| 4 | `BE(1,7;1,7)` | `BE_BASE_SAMEPAIR` | ∅ |
| 5 | `VP(2,17)` | `VP_VP_WIP_RIGHT` | 6, 13 |
| 6 | `VP(3,11)` | `VP_WI_CASE1` | 7, 10 |
| 7 | `WI(4,6)` | `WI_EXTEND_UNPAIRED` | 8 |
| 8 | `WI(4,5)` | `WI_EXTEND_UNPAIRED` | 9 |
| 9 | `WI(4,4)` | `WI_BASE_SINGLE` | ∅ |
| 10 | `WI(8,10)` | `WI_EXTEND_UNPAIRED` | 11 |
| 11 | `WI(8,9)` | `WI_EXTEND_UNPAIRED` | 12 |
| 12 | `WI(8,8)` | `WI_BASE_SINGLE` | ∅ |
| 13 | `WIP(12,16)` | `WIP_BASE_V` | 14 |
| 14 | `V(12,16)` | `V_HAIRPIN` | ∅ |

教科書の対応する実行記録は `03_grammar/v58a2b796/vp_EB.txt`。
