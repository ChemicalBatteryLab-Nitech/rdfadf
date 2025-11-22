# rdfadf パッケージ（日本語版 README）

## 1. 概要（Summary）

`rdfadf` は、結晶構造およびアモルファス構造に対して  
**RDF（Radial Distribution Function）** および  
**ADF（Angular Distribution Function）** を高精度に計算する Python パッケージです。

本パッケージは以下の特徴を持ちます：

- **周期境界条件（PBC）を完全に考慮した RDF/ADF 計算**
- ASE (`ase.Atoms`) に標準対応
- 中心原子・近接原子を  
  **元素 / ブロック（s,p,d,f）/ 族（G#）/ 周期（P#）/ All**  
  で柔軟にフィルタ可能
- RDF は  
  - 中心原子数で正規化  
  - 4πr² 正規化にも対応  
  - 累積カウント（CN の直接計算）を出力可能
- ADF は  
  - r_min / r_max の距離範囲で近接原子を選択  
  - 中心原子数で正規化  
  - Gaussian smoothing に対応
- 研究用途・論文生成・局所構造解析に最適化

## 2. インストール（Install）

```bash
pip install .
```

アンインストール：

```bash
pip uninstall rdfadf -y
```

## 3. 入力（Inputs）

```python
from ase.io import read

atoms = read("POSCAR")
atoms = read("structure.cif")
atoms = read("output.traj")
```

## 4. ラベル指定

| 種類 | 例 | 意味 |
|------|------|------|
| 元素記号 | "Li", "O", "Co" | 特定の元素 |
| All | "All" | 構造中の全元素 |
| ブロック | "s", "p", "d", "f" | ブロック分類 |
| 族 | "G16" | 16 族 |
| 周期 | "P3" | 第 3 周期 |

## 5. RDF 関数

```
compute_general_rdf(atoms, center_label, neighbor_label, ...)
```

出力：

```
r,
g_raw_norm,
g_smooth,
g_raw2,
g_smooth2,
cumulative_raw,
cumulative_smooth
```

## 6. ADF 関数

```
compute_general_adf(atoms, center_label, neigh1_label, neigh2_label, ...)
```

## 7. RDF 使用例

```python
from rdfadf import compute_general_rdf
from ase.io import read

atoms = read("POSCAR")

r, g_raw, g_smooth, g_raw2, g_smooth2, cn_raw, cn_smooth =     compute_general_rdf(
        atoms,
        center_label="Co",
        neighbor_label="O",
        rcut=6.0,
        bins=200,
        sigma=2.0
    )
```

## 8. ADF 使用例

```python
theta, adf_raw, adf_smooth = compute_general_adf(
    atoms,
    center_label="Co",
    neigh1_label="O",
    neigh2_label="O",
    rcut=4.0,
    r_min=1.0,
    r_max=3.0,
    bins=180,
    sigma=2.0
)
```

## 9. 注意点

- ADF は RDF より計算コストが高い場合あり
- atoms.pbc = True が必要
- 大規模 MD ではフレームごとの処理が必要

## 10. ライセンス

MIT（推奨）

