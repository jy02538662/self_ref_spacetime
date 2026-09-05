# Self-Referential Spacetime — Numerical Closure Tests

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](#)
[![numpy](https://img.shields.io/badge/numpy-%3E%3D1.24-green.svg)](requirements.txt)
[![scipy](https://img.shields.io/badge/scipy-%3E%3D1.10-green.svg)](requirements.txt)

「自指时空」可证伪数值实验合集：厄米耦合矩阵 D（D_ij = D_ji*）如何读出距离、号差、拓扑、空间维数与量子骨架。

> **状态**：研究原型 / **封存期**（约至 2026-11）。数值是可复现计算证据，**不是**已完成物理证明；不声称证明 3+1。  
> **联系**：王超 · 1186306891@qq.com  
> **理论**：vault《量子潮水理论行动指南 v3》《断裂与自指：量子性的自指来源》  
> **公式约定（GitHub 稳妥）**：正文与表格用 Unicode/ASCII；关键公式只用独立 `$$` 块。行内尽量不写 `$...$`，避免渲染错位。

**数据约定**：下文「归档数据」一律取自仓库内 `experiments/*_last_run.json`（复跑后会覆盖）。无 json 的条目标为「文档记录 / 需重跑」。

---

## 理论背景（三分钟）

窄命题：给定厄米 D，能否用同一套数值门控读出几何距离、洛伦兹号差、拓扑荷与量子化骨架——不预设 3+1 坐标。

最小作用量（环靶距离 L_ij；d_ij = 边长 1/|D_ij| 的最短路）：

$$
S = \mathrm{Tr}(D^{2}) + \lambda \sum_{i \lt j} (d_{ij} - L_{ij})^{2}
$$

后续探针：磁通、Tr(D⁴)、Dirac、Hopf / Jones / 纠缠熵。代码回答的是可证伪数值问题，不是「时空公式已推出」。

## 五分钟入门

**1）冒烟**

```bash
cd self_ref_spacetime
pip install -r requirements.txt
python -m tests.test_distance
python -m experiments.exp1_ring_neighbor
```

期望：`exp1_last_run.json` 中 `summary.pass` 相关字段显示约 3/8 种子进环盆地。

**2）核心负结果：号差 ≠ 虚相位**

```bash
python -m experiments.exp2_phase_from_D4
python -m experiments.exp3_directed_signature
```

Exp2 fail（Tr(D⁴) 对相位盲）；Exp3 pass（号差需 Dirac / 有向结构）。

**3）对照论文时优先看**

| 用途 | 入口 |
|------|------|
| 总判词 + 归档数字 | [门控板与数据总表](#gate-board) |
| 逐实验细表 | [实验数据详表](#experimental-data) |
| 论文式子 / 长推导 | vault 笔记或论文 PDF（本 README 不承担讲义） |

---

## 目录

1. [Highlights](#highlights)
2. [门控板与数据总表](#gate-board)
3. [已闭环与钉死](#closed)
4. [实验数据详表](#experimental-data)（完整归档）
5. [物理结论摘要](#physics-notes)
6. [Quick start](#quick-start)
7. [仓库目录](#repository-layout)
8. [范围边界](#scope-boundary)

---

## Highlights

- 厄米 D + 最短路距离 + L-BFGS-B
- Exp1：环近邻可涌现（多峰；归档 3/8）
- Exp2–3：虚相位 ≠ 号差；号差来自有向 / Dirac
- Exp4–5：拓扑荷 → 零模（绕数 / 陈数）
- Exp6a：Hopf 荷数值可读（Q=1 好；归档 Q=2 未收敛到 4）
- 谱维可读；复现几何**不**自发选三维
- 正向：辫子词 → Jones；量子化 = JW 截断 / 有限 N Chebyshev
- 弦网「自指 → SU(2)_k」路线**封存**（已知工具复现）

---

## Gate board

判词以脚本逻辑为准；**数字以 `*_last_run.json` 为准**。

| 门控 | 脚本 | 判词 | 归档关键数字 |
|------|------|------|----------------|
| Exp1 环近邻 | `exp1_ring_neighbor` | **pass** | 3/8 环盆地；best S≈11.903，geo≈0.00313，hit=1.0 |
| Exp2 Tr(D⁴) | `exp2_phase_from_D4` | **fail** | 5 图案 TrD4 全 = 36.0；free 不优于 real |
| Exp3 Dirac 号差 | `exp3_directed_signature` | **pass** | 欧氏 n_neg=0；洛伦兹 n_pos=n_neg=48，不定 |
| Exp4 JR | `exp4_jackiw_rebbi` | **pass** | 均匀 0 零模；畴壁 2 零模，E_min≈9e-18 |
| Exp5 陈数 | `exp5_chern` | **pass** | Chern≈8 → 16 零模 |
| Exp1 v2 磁通 | `exp1_v2_flux` | **pass（负倾向）** | μ=0→10：flux 9.35→1e-7；仍非实矩阵；π-frac≈0.3–0.44 |
| Exp1 v3 曲率 | `exp1_v3_flux_compare` | **pass** | sin→Z₂；wilson→0；cos→quarter_frac=1 |
| Exp1 v4 MCMC | `exp1_v4_mcmc` | **负结果** | β↑ flux↓；β=10 时 π-frac≈0.51；accept≈0.001（未充分热化） |
| 迹反常 | `exp_trace_anomaly` | **边界/负** | 0 维矩阵：g≠0 有斜率，但非时空 RG 迹反常 |
| 量子化自指 | `exp_quantization_selfref` | **pass** | [X,P]≠0；min break=ε²=0.25；坍缩 0.346→1.0 |
| 空间·暖/干净 | `exp_space_3d`（clean/warm） | **pass** | n=27：hit_top6=1.0，geo≈9.8e-4 |
| 空间·冷随机 | `exp_space_3d` | 局部极小 | mean hit≈0.73 |
| 空间·退火 | `exp_space_3d_anneal` | 改善未干净 | mean hit≈0.82，best≈0.90 |
| 空间·冷启动纯 TrD² | `exp_cold_start` | **负结果** | TrD2→~0；谱维塌缩 |
| Exp6a Hopf | `exp6a_hopf_charge` | **部分** | Q_H=(0, 0.99983, 3.604)；Q2 未过关 → pass=false |
| 谱流 | `exp_spectral_flow` | **未达预期** | flow=0（期望 1）；crossings=4；E0≈7e-16 |
| η framing | `exp_eta_framing` | 文档记录 | 无 json；见详表 |
| 谱维对比 | `exp_spectral_dim_compare` | 文档记录 | 无 json；见详表 |
| 维度对比 | `exp_space_dim_compare` | **负结果** | 无 json；文档：暖全 hit=1，冷偏向低维 |
| 弦网 sanity | `exp_string_net_condensation` | **封存** | TEE k=1..5 ≈ 0.693…2.923 |
| Born | `exp_born_deviation` | **pass** | SO(3) 偏≈6e-16；偏好方向偏≈0.50 |
| 纠缠面积律 | `exp_entanglement_area_law` | **pass** | S_loc~0.15 ln L；S_rand~0.38 L |
| 面积律+维 | `exp_entanglement_area_law_dim` | **pass** | 斜率 1.25 / 2.25 |
| 面积=局域 | `exp_area_law_locality` | **pass** | 链面积 vs 随机体积 |
| 裸重连 | `exp_bare_reconnection` | **负结果** | N=30，1e5 步无自组织 |
| 相位联络 | `exp_phase_connection` | **pass（负倾向）** | S_top → −4；平坦 |
| Wilson Q | `exp_wilson_invariant` | **pass** | max_dev=0 |
| SU(2) Wilson | `exp_su2_wilson` | **pass** | skein 成立 |
| 群交换子 L | `exp_group_commutator` | **pass+边界** | SU(2) L_mean≈0.748；U(1) L=0 |
| YBE | `exp_yang_baxter` | **pass** | 误差 ~1e-16 |
| Jones | `exp_jones` | **pass** | 平凡 ≠ Hopf |
| 配对↔交叉 | `exp_pairing_to_braid` | **pass** | skein 一致 |
| 辫子→Jones | `exp_braid_word_to_jones` | **pass** | 21 checks |
| A 从哪来 | `exp_quantization_condition` | **pass（边界）** | 经典 A 自由 |
| Jones–Wenzl | `exp_jones_wenzl` | **pass** | 截断逼单位根 |
| 关系→拓扑 | `exp_relation_to_topology` | **pass（边界）** | 投影幂等；截断=Chebyshev 零 |
| 有限 N | `exp_finite_quantization` | **pass** | Chebyshev 恒等式成立 |
| 纯谱选 3+1 | — | out of scope | — |

---

## Closed

- **已闭环**：量子化自指；欧氏环涌现（Exp1）；号差（Exp3）；零模（Exp4/5）
- **已钉死（归档）**：干净三维近邻 hit=1.0（n=27）；Q_H(Q=1)≈0.99983；面积律 vs 体积律；辫子→Jones 正向链
- **负结果钉子**：虚相位不给号差；经典曲率/热难产非平凡拓扑；裸重连无组织；纯 Tr(D²) 不选维；复现几何偏向低维（文档）
- **归档未过关**：Exp6a 的 Q=2；谱流 net flow≠1——复现时勿把 README 旧口头「全过」当真，以 json 为准

---

## Experimental data

以下数字来自当前仓库归档 json（除非标明「无 json」）。复现时对比同名字段即可。

### Exp1 — 环近邻（`exp1_last_run.json`）

参数：n=6，lam=30，seeds 0–7。

| 项 | 值 |
|----|-----|
| pass | true（判据：环盆地比例 ≥0.25 等） |
| n_ring_basins | **3 / 8**（fraction **0.375**） |
| best_S | **11.903131556922** |
| best_geo | **0.003130630923** |
| best_neighbor_hit_rate | **1.0** |
| best_strength_ratio | **6.23e-6** |

环盆地种子（hit=1）：0、2、4（S≈11.903–11.904）。失败种子 S 量级 10²、hit 0.5–0.83。

```bash
python -m experiments.exp1_ring_neighbor --n 6 --lam 30 --seeds 8
```

### Exp2 — Tr(D⁴) 虚相位（`exp2_last_run.json`）

| 项 | 值 |
|----|-----|
| pass | **false** |
| Protocol A | 全图案 TrD2=12，**TrD4=36.0**，S=48；imag 不优于 all_real |
| Protocol B | 好种子 0,2,4；frac_free_lower_TrD4=**0**；ΔTrD4(free−real)≈+4.5~+4.7 |

解析：等模环 Tr(D⁴)=6n（n=6→36），与边相位无关。

### Exp3 — Dirac 号差（`exp3_last_run.json`）

| 构型 | n_pos | n_neg | n_zero | indefinite |
|------|------:|------:|-------:|:----------:|
| Euclidean | 112 | 0 | 16 | false |
| Lorentzian | 48 | 48 | 32 | **true**（eig 约 ±3.411） |

pass=**true**。

### Exp4 — Jackiw–Rebbi（`exp4_last_run.json`）

n_sites=400。

| 构型 | n_zero | eig_min_abs |
|------|-------:|-------------:|
| 均匀质量 | 0 | ≈1.000 |
| 畴壁 | **2**（位点 200,200） | ≈8.87e-18 |

fermion_doubling=2；pass=**true**。

### Exp5 — 陈数（`exp5_last_run.json`）

Nx=Ny=30，plaquettes=841。

| 磁场 | chern_flux | n_zero |
|------|-----------:|-------:|
| Φ=0 | 0 | 0 |
| Φ≠0 | ≈8.000 | **16** |

pass=**true**。

### Exp1 v2 — 磁通惩罚（`exp1_v2_last_run.json`）

| μ | mean_flux | mean_pi_frac | mean abs(sin φ) |
|---|----------:|-------------:|-------------:|
| 0 | 9.354 | 0.300 | 0.590 |
| 0.1 | 0.00524 | 0.438 | 0.635 |
| 1 | 1.65e-5 | 0.438 | 0.635 |
| 10 | 1.11e-7 | 0.381 | 0.635 |

结论：flux 被压；**phases_flattened_to_real=false**（压到 0/π 的 Z₂，不是全实）。

### Exp1 v3 — 曲率定义（`exp1_v3_last_run.json`）

| kind | mean_pi_frac | mean_quarter_frac | 行为 |
|------|-------------:|------------------:|------|
| sin | 0.4375 | — | 允许 π 通量（Z₂） |
| wilson | **0** | — | 压到零磁通 |
| cos | 0 | **1.0** | 全锁 ±π/2 |

### Exp1 v4 — MCMC（`exp1_v4_last_run.json`）

seeds 0–3。

| β | π_frac 均值 | mean_flux | accept |
|---|------------:|----------:|-------:|
| 0.01 | 0.313 | 9.537 | 0.995 |
| 0.1 | 0.332 | 9.278 | 0.945 |
| 1 | 0.389 | 6.817 | 0.466 |
| 10 | 0.513 | 0.499 | **0.0014** |

`degeneracy_intact=false`（脚本判据：跨 β 是否保持 π_frac~0.5）。物理解读仍是：低温未可靠打破 0/π；β=10 接受率过低，**未充分热化**。

### 迹反常（`exp_trace_anomaly_last_run.json`）

| g | resid_slope_vs_ln λ |
|---|--------------------:|
| 0 | −20.0 |
| 0.3 | −11.33 |
| 1.0 | −7.69 |

analytic g=0 Tr ln M₂ ≈ 4.15888。定位：0 维玩具，不是时空迹反常证明。

### 量子化自指（`exp_quantization_selfref_last_run.json`）

| 子实验 | 结果 |
|--------|------|
| A | commutator_nonzero；matches_2i_sigma2 |
| B | eps=0.5；break_energies=[0, 0.25, 1, 2.25, 4]；min_nonzero=**0.25=ε²** |
| C | overlap 0.346 → **1.0**（坍缩） |

### 空间涌现

**干净 / 暖启动类（`exp_space_3d_clean_last_run.json`）** — n=27（3³），lam=30：

| 项 | 值 |
|----|-----|
| mean_hit_top6 | **1.0** |
| mean_geo | **9.80e-4** |
| mean_degree | 3.0 |

**冷随机（`exp_space_3d_last_run.json`）** — n=8 网格种子：mean_hit≈**0.729**，mean_fp≈0.25。

**退火（`exp_space_3d_anneal_last_run.json`）** — 5 seeds：mean_hit_top6≈**0.817**，best≈**0.904**，n_clean=1。

**冷启动纯 Tr(D²)（`exp_cold_start_last_run.json`）**：TrD2 → ~1e-14；谱维塌缩（负结果）。

**维度对比（`exp_space_dim_compare`）**：无 json。文档记录：2D/3D/4D 暖启动均 hit=1.0；冷启动 excess 偏向低维 → **复现几何不选三维**。

### Exp6a — Hopf（`exp6a_last_run.json`）

| 场 | 归档 Q_H | 期望 | 过关？ |
|----|--------:|------|:------:|
| 平凡 Q=0 | **0.0** | 0 | yes |
| Hopf Q=1 | **0.9998305787** | 1 | yes |
| Q=2 | **3.6043066465** | 4 | **no** |

summary.pass=**false**（卡在 Q=2；divB 质量差）。  
说明：更大盒子（如 N=128,L=10）曾口头报过 Q1≈0.99996、Q2≈4，**当前仓库未归档该次 json**；对外引用请以本表或自行重跑并提交新 json。

```bash
python -m experiments.exp6a_hopf_charge
```

Hopf 积分（文档公式，独立块）：

$$
Q_{H}=\frac{1}{16\pi^{2}}\int A\cdot B\,d^{3}x
$$

### 谱流（`exp_spectral_flow_last_run.json`）

| 项 | 值 |
|----|-----|
| flow | **0**（expected_flow=1） |
| crossings | 4 |
| n_levels | 400 |
| e0_at_lam0 | ≈7.09e-16 |

未达「净谱流=绕数」整数判据；畴壁零模能量接近 0。留 Exp6b。

### η framing（`exp_eta_framing`）— **无 json**

文档记录（需重跑存档）：1D 演示 flow=1 与 η→±1/2，合成 3/2。复现后请写出 `exp_eta_framing_last_run.json`。

### 谱维

**`exp_spectral_dim_last_run.json`**（图拉普拉斯玩具）：

| 图 | measured | expected |
|----|--------:|--------:|
| ring1d | 1.183 | 1 |
| grid2d | 2.396 | 2 |
| grid3d | 3.849 | 3 |

**`exp_spectral_dim_compare`** — 无 json。文档记录（周期超立方热核拟合）：d_s≈2.009 / 3.014 / 4.018。需重跑存档。

### 弦网 sanity（`exp_string_net_condensation_last_run.json`）— 封存

TEE = log D：

| k | TEE |
|--:|----:|
| 1 | 0.693 |
| 2 | 1.386 |
| 3 | 1.979 |
| 4 | 2.485 |
| 5 | 2.923 |

环气体：T:8→0 时密度 0.002→0.5，condensed=true。加权真空 p₀：k=2→0.25；k=3→0.138；k=4→0.083。  
**未**让自指 D 选出 SU(2)；已知 Levin–Wen 工具复现。

### Born（`exp_born_deviation_last_run.json`）

| Part | 结果 |
|------|------|
| A SO(3) | max_dev ≈ **6.11e-16** |
| B 偏好方向 | max_dev ≈ **0.497** |
| C ε 扫描 | 偏差随 ε 近似线性（0→0 … 1→1.656） |

criterion_passes=**true**。

### 纠缠

**面积 vs 体积（`exp_entanglement_area_law_last_run.json`）** N=400：

| L | S_local | S_random |
|--:|--------:|---------:|
| 10 | 0.8465 | 6.814 |
| 40 | 1.0864 | 25.773 |
| 200 | 1.2848 | 77.580 |

斜率：local_vs_logL≈**0.1486**；random_vs_L≈**0.3797**。

**维数（`exp_entanglement_area_law_dim_last_run.json`）**：2D slope≈**1.245**；3D slope≈**2.254**。

**局域性（`exp_area_law_locality_last_run.json`）**：链 → 面积（对数）；随机图 p=0.3 → 体积。

### 裸重连（`exp_bare_reconnection_last_run.json`）

N=30，steps=1e5。闭弦计数/长度来回跳，无「少数大弦主导」趋势。负结果。

### 相位 → Wilson → Jones → 量子化

| 实验 | 归档要点 |
|------|----------|
| `exp_phase_connection` | 3×3，S_top 0.74→**−4**；flux_flat=true |
| `exp_wilson_invariant` | Q_wilson=Q_direct；max_deviation=**0**；invariant=true |
| `exp_su2_wilson` | 迹随重连变；**skein_identity_holds=true** |
| `exp_group_commutator` | SU(2) L_mean≈**0.748**（100% 非零）；U(1) max_L=0 |
| `exp_yang_baxter` | qybe_holds & braid_holds；err~1e-16（k=1..5） |
| `exp_jones` | 平凡 ≠ Hopf；k=1 Jones triv≈(−1.732,0)，Hopf≈(0,1) |
| `exp_pairing_to_braid` | skein 一致；e²−de=0 |
| `exp_braid_word_to_jones` | **checks_passed=true**，n_checks=**21**（含 R-II 等） |
| `exp_quantization_condition` | classical_A_is_free=true；−d(k=1..5)≈1.732…1.950 |
| `exp_jones_wenzl` | hecke 不逼单位根；JW 截断零点随 k 移位 |
| `exp_relation_to_topology` | 谱投影幂等；truncation=Chebyshev 零 |
| `exp_finite_quantization` | chebyshev_identity_holds=true；k=N−1 |

---

## Physics notes

**号差**：厄米 D ⇒ D² 本征值 ≥0；虚相位不改 |λ|²。号差需要有向/非厄米反对称或 Dirac γ（Exp3）。2×2：仅 w=−z（非厄米）使 zw 为负。

**Tr(D⁴)**：等模环上 =6n，相位盲（Exp2）。

**拓扑阶梯**：1D 绕数 / 2D 陈数 → 零模（指标）；3D Hopf = odd CS，不是指标（Exp6a）。

**三维**：暖/干净近邻稳定 ≠ 冷启动自发选三维；维数压力需拓扑层。

**量子化开关**：经典不固定交叉振幅 A；JW 截断 / 有限 N Chebyshev 才逼单位根。

---

## Quick start

```bash
cd self_ref_spacetime
pip install -r requirements.txt

python -m tests.test_distance
python -m experiments.exp1_ring_neighbor
python -m experiments.exp2_phase_from_D4
python -m experiments.exp3_directed_signature
python -m experiments.exp4_jackiw_rebbi
python -m experiments.exp5_chern
python -m experiments.exp1_v2_flux
python -m experiments.exp1_v3_flux_compare
python -m experiments.exp1_v4_mcmc
python -m experiments.exp_trace_anomaly
python -m experiments.exp_quantization_selfref
python -m experiments.exp_space_3d
python -m experiments.exp_space_3d_anneal
python -m experiments.exp_cold_start
python -m experiments.exp_spectral_flow
python -m experiments.exp6a_hopf_charge
python -m experiments.exp_space_dim_compare
python -m experiments.exp_eta_framing
python -m experiments.exp_spectral_dim
python -m experiments.exp_spectral_dim_compare
python -m experiments.exp_string_net_condensation
python -m experiments.exp_born_deviation
python -m experiments.exp_entanglement_area_law
python -m experiments.exp_entanglement_area_law_dim
python -m experiments.exp_area_law_locality
python -m experiments.exp_bare_reconnection
python -m experiments.exp_phase_connection
python -m experiments.exp_wilson_invariant
python -m experiments.exp_su2_wilson
python -m experiments.exp_group_commutator
python -m experiments.exp_yang_baxter
python -m experiments.exp_jones
python -m experiments.exp_pairing_to_braid
python -m experiments.exp_braid_word_to_jones
python -m experiments.exp_quantization_condition
python -m experiments.exp_jones_wenzl
python -m experiments.exp_relation_to_topology
python -m experiments.exp_finite_quantization
```

依赖：`numpy>=1.24`，`scipy>=1.10`。结果写入对应 `experiments/<name>_last_run.json`。

**论文对照建议**：论文表/式 → 本 README 门控行 → 同名 json 字段。缺 json 的脚本（η、谱维对比、维度对比）复现后请提交归档文件。

---

## Repository layout

```
self_ref_spacetime/
  src/
    algebra.py      # 厄米 D：模长 + 相位
    distance.py     # 1/|z| → 最短路
    action.py       # Tr(D²) + λ·geo + λ₄ Tr(D⁴)
    metrics.py      # 诊断
    optimize.py     # L-BFGS-B
    flux.py         # 闭合环磁通
  experiments/      # 门控脚本 + *_last_run.json
  tests/
  requirements.txt
  README.md
```

---

## Scope boundary

| 主张 | 立场 |
|------|------|
| 环近邻可涌现 | ✅ 归档 3/8 |
| 虚相位 ⇒ 洛伦兹号差 | ❌ Exp2 + 代数 |
| Dirac/有向 ⇒ 号差 | ✅ Exp3 |
| 拓扑荷 ⇒ 零模 | ✅ Exp4/5 |
| Hopf Q=1 可读 | ✅ 归档 ≈0.99983 |
| Hopf Q=2→4（本归档） | ❌ 当前 json 未过 |
| D 自发选三维 | ❌ 复现几何负结果 |
| 经典作用量/热/裸重连产拓扑 | ❌ |
| 自指 D → SU(2)_k 弦网 | ⛔ 封存 |
| 证明 3+1 | out of scope |
