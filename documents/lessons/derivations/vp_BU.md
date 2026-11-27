# VP の片側に枝がある場合

固定対を丸括弧、追加対を角括弧で示す。同じ種類の開き括弧と閉じ括弧を、内側から順に対応させて読む。入力で決まっている固定対は、結合時に再び選び直さない。

```text
配列　　GGGAAACGAAACAAACAC
固定対　(..........)......
追加対　.[[...][.......].]
完成形　([[...][...)...].]
```

![具体的な導出木](vp_BU.svg)

[図を単独で拡大](vp_BU.svg)

| # | 状態 | 選んだ規則 | 子の番号（順序付き） |
|---:|---|---|---|
| 1 | `W(1,18)` | `W_SPLIT_WMB` | 2 |
| 2 | `WMB(1,18)` | `WMB_DIRECT_WMBP` | 3 |
| 3 | `WMBP(1,18)` | `WMBP_SPLIT_BE_WI_VP` | 4, 5 |
| 4 | `BE(1,12;1,12)` | `BE_BASE_SAMEPAIR` | ∅ |
| 5 | `VP(2,18)` | `VP_WIP_VPR` | 6, 8 |
| 6 | `WIP(3,7)` | `WIP_BASE_V` | 7 |
| 7 | `V(3,7)` | `V_HAIRPIN` | ∅ |
| 8 | `VPR(8,17)` | `VPR_SPLIT_VP_BASEPAIR` | 9 |
| 9 | `VP(8,16)` | `VP_WI_CASE1` | 10, 13 |
| 10 | `WI(9,11)` | `WI_EXTEND_UNPAIRED` | 11 |
| 11 | `WI(9,10)` | `WI_EXTEND_UNPAIRED` | 12 |
| 12 | `WI(9,9)` | `WI_BASE_SINGLE` | ∅ |
| 13 | `WI(13,15)` | `WI_EXTEND_UNPAIRED` | 14 |
| 14 | `WI(13,14)` | `WI_EXTEND_UNPAIRED` | 15 |
| 15 | `WI(13,13)` | `WI_BASE_SINGLE` | ∅ |

教科書の対応する実行記録は `03_grammar/v58a2b796/vp_BU.txt`。
