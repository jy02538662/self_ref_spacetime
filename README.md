# Self-Referential Spacetime — Numerical Closure Tests

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](#)
[![numpy](https://img.shields.io/badge/numpy-%3E%3D1.24-green.svg)](requirements.txt)
[![scipy](https://img.shields.io/badge/scipy-%3E%3D1.10-green.svg)](requirements.txt)

「自指时空」的可证伪数值实验合集：一个自反性关系算子 $D$（$D_{ij}=D_{ji}^{*}$）如何涌现出距离、号差、拓扑、空间维数与量子骨架。

> **状态**：研究原型 / **封存期**（约至 2026-11，考试复习）。数值结果是可复现的计算证据，**不是**已完成的物理证明；不声称证明 3+1。  
> **联系**：王超 · 1186306891@qq.com  
> **理论路线**：vault《量子潮水理论行动指南 v3（自指主线版）》《断裂与自指：量子性的自指来源》

---

## 理论背景（三分钟）

本仓库检验的是一条窄命题：给定厄米耦合矩阵 $D$，能否用**同一套数值门控**读出几何距离、洛伦兹号差、拓扑荷与量子化骨架——而不预设 3+1 坐标。

最小作用量（环型靶距离 $L_{ij}$）为

$$
S = \mathrm{Tr}(D^{2}) + \lambda \sum_{i \lt j} (d_{ij} - L_{ij})^{2}
$$

其中 $d_{ij}$ 由边长 $1/|D_{ij}|$ 的最短路给出。后续实验在此骨架上加磁通项、$\mathrm{Tr}(D^{4})$、Dirac 结构、Hopf / Jones / 纠缠熵等探针。

代码要回答的，不是「时空公式已推出」，而是一串可证伪的数值问题：环近邻是否自发涌现、虚相位能否给号差、拓扑荷能否读出零模、三维是否被复现几何选中、经典动力学能否生成拓扑、量子化开关在哪里。

## 五分钟入门

新人按下面三条走即可；门控总表与分节在后文。

**1）冒烟（本机）**

```bash
cd self_ref_spacetime
pip install -r requirements.txt
python -m tests.test_distance
python -m experiments.exp1_ring_neighbor
```

看到 Exp1 summary 里 `pass: true`（约 3/8 种子进干净环盆地）即链路通。

**2）核心负结果：号差不来自虚相位**

```bash
python -m experiments.exp2_phase_from_D4
python -m experiments.exp3_directed_signature
```

Exp2 **fail**（Tr(D⁴) 对相位盲）；Exp3 **pass**（号差需 Dirac / 有向结构，不是虚相位）。

**3）模块入口**

| 想看什么 | 先打开 / 先跑 |
|----------|----------------|
| 环近邻涌现 | `experiments/exp1_ring_neighbor.py` |
| 相位 / 磁通 / MCMC | Exp1 v2–v4 · [闭合环磁通](#闭合环磁通项exp1-v2--v3--v4) |
| 号差与拓扑零模 | Exp3–5 · [理论缺口](#理论缺口号差为什么读不出) |
| 空间 / 谱维 / Hopf | Exp6a 等 · [Exp6a](#exp6a-hopf-荷--odd-chernsimons) |
| 纠缠 / 局域性 | [纠缠熵](#纠缠熵面积律-vs-体积律) |
| 辫子 / Jones / 量子化 | [辫子词→Jones](#辫子词--jones-多项式) |
| 设计稿 Exp7 | [`exp7_design.md`](exp7_design.md) |

**阅读顺序建议**：本节 → [门控板](#gate-board) → [已闭环](#已闭环与钉死) → 需要时再进分节细文。

---

## 目录

0. [理论背景（三分钟）](#理论背景三分钟) · [五分钟入门](#五分钟入门)
1. [项目要点](#highlights)
2. [门控板](#gate-board)
3. [已闭环与钉死](#已闭环与钉死)
4. [Pass / Fail 判据（Exp1–5）](#pass-criterion-exp1)
5. [理论缺口：号差](#理论缺口号差为什么读不出)
6. [Exp2 解析定理](#exp2-的解析定理trd4-对相位盲)
7. [闭合环磁通（v2/v3/v4）](#闭合环磁通项exp1-v2--v3--v4)
8. [量子化自指 + 空间涌现](#量子化自指--空间涌现)
9. [Exp6a Hopf](#exp6a-hopf-荷--odd-chernsimons)
10. [谱流 / η / 谱维 / 维度对比](#谱流前置)
11. [弦网 / Born / 纠缠 / 裸重连](#弦网凝聚-sanity-check)
12. [相位联络 → Jones → 量子化](#相位当联络)
13. [快速开始](#quick-start)
14. [仓库目录](#repository-layout)
15. [范围边界](#scope-boundary)

---

## Highlights

- 厄米 $D$ 参数化 + 最短路距离 + 作用量极小化（L-BFGS-B）
- Exp1：环靶距离下近邻图可涌现（多峰景观）
- Exp2–3：钉死「虚相位 ≠ 号差」；号差来自有向 / Dirac 结构
- Exp4–5 / Exp6a：拓扑荷 → 零模；Hopf 荷 = odd Chern–Simons（数值）
- 谱维数可读出 2/3/4；复现几何**不**自发选三维（负结果）
- B 路线正向：辫子词 → Jones；量子化 = Jones–Wenzl 截断 / 有限 $N$ 恒等式
- 诚实封存弦网「自指 → $\mathrm{SU}(2)_{k}$」拼装（已知工具复现，非推进）

## Gate board

| 门控 | 脚本 | 判词 | 一句话 |
|------|------|------|--------|
| Exp1 环近邻 | `exp1_ring_neighbor` | **pass** | 约 3/8 种子进干净环盆地 |
| Exp2 Tr(D⁴) 虚耦合 | `exp2_phase_from_D4` | **fail** | 等模环上相位不降 TrD4 |
| Exp3 Dirac 号差 | `exp3_directed_signature` | **pass** | 欧氏全非负；洛伦兹不定 (−,+) |
| Exp4 Jackiw–Rebbi | `exp4_jackiw_rebbi` | **pass** | 绕数 1 → 2 零模 |
| Exp5 陈数 | `exp5_chern` | **pass** | 陈数 8 → 16 零模 |
| Exp1 v2 磁通 | `exp1_v2_flux` | **pass** | 压到 0/π 的 Z₂ 通量 |
| Exp1 v3 曲率定义 | `exp1_v3_flux_compare` | **pass** | 曲率定义决定磁通零点 |
| Exp1 v4 MCMC | `exp1_v4_mcmc` | **负结果** | 经典热不打破 Z₂ |
| 迹反常 | `exp_trace_anomaly` | **负结果** | 0 维矩阵无量纲嬗变 |
| 量子化自指 | `exp_quantization_selfref` | **pass** | [X,P]≠0、ε²、投影坍缩 |
| 空间·暖启动 | `exp_space_3d --warm-start` | **pass** | hit=1.0，三维近邻稳定 |
| 空间·冷/退火 | `exp_space_3d` / `_anneal` | 局部极小 | hit≈0.62 → 0.82–0.90 |
| 空间·维度对比 | `exp_space_dim_compare` | **负结果** | 复现几何偏向低维 |
| 空间·冷启动 | `exp_cold_start` | **负结果** | 纯 Tr(D²) 不选维 |
| Exp6a Hopf | `exp6a_hopf_charge` | **pass** | Q_H≈0.99996；Q=2→4 |
| 谱流 | `exp_spectral_flow` | **pass** | 畴壁零模 + 穿越 |
| 3/2 = 谱流 + η | `exp_eta_framing` | **pass** | 1+(±1/2) 演示 |
| 谱维对比 | `exp_spectral_dim_compare` | **pass** | d_s≈2.01 / 3.01 / 4.02 |
| 弦网 sanity | `exp_string_net_condensation` | **封存** | 已知工具复现，非推进 |
| Born 偏离 | `exp_born_deviation` | **pass** | SO(3) 协变；外部方向可测 |
| 纠缠面积律 | `exp_entanglement_area_law` | **pass** | 局域 ln L vs 随机体积律 |
| 面积律+维度 | `exp_entanglement_area_law_dim` | **pass** | 指数 ≈ d−1 |
| 面积律=局域 | `exp_area_law_locality` | **pass** | 费米面 vs 随机图 |
| 裸重连 | `exp_bare_reconnection` | **负结果** | 无偏好不组织 |
| 相位联络 | `exp_phase_connection` | **pass（负倾向）** | U(1) 麦克斯韦→平坦 |
| Wilson 不变量 | `exp_wilson_invariant` | **pass** | Q 严格重连守恒 |
| SU(2) Wilson | `exp_su2_wilson` | **pass** | skein / Fierz 成立 |
| 群交换子 L | `exp_group_commutator` | **pass + 边界** | 非对易度量 ≠ 链接数 |
| Yang–Baxter | `exp_yang_baxter` | **pass** | QYBE + braid |
| Jones | `exp_jones` | **pass** | 平凡 vs Hopf 可分 |
| 配对↔交叉 | `exp_pairing_to_braid` | **pass** | skein：交叉=叠加 |
| 辫子词→Jones | `exp_braid_word_to_jones` | **pass** | R-II / R-III 钉死 |
| A 从哪来 | `exp_quantization_condition` | **pass（边界）** | 经典不固定 A |
| Jones–Wenzl | `exp_jones_wenzl` | **pass** | 截断逼出单位根 |
| 关系→拓扑 | `exp_relation_to_topology` | **pass（边界）** | 免费桥 + 付费桥 |
| 有限 N→单位根 | `exp_finite_quantization` | **pass** | Chebyshev 恒等式 |
| 纯谱 / 3+1 选择 | — | out of scope | — |

## 已闭环与钉死

- **已闭环**：量子化自指（断裂+自指→量子性）+ 欧氏涌现（Exp1）+ 号差（Exp3）+ 拓扑零模（Exp4 绕数、Exp5 陈数）
- **已钉死**：三维近邻是 D 的稳定解（暖启动 hit=1.0）；Hopf 荷 ≈ 1（Exp6a）；谱维数可读 2/3/4；总拓扑荷 3/2 的 1D 演示（η framing）
- **负结果钉子**：复现几何不选三维；经典曲率/热/裸重连都不生成非平凡拓扑；虚相位不给号差
- **档 2 维度对比（负结果）**：复现几何不选三维、偏向低维

---

## Pass criterion (Exp1)

- $N=6$ 或 $8$；环距离 $L_{ij}=\min(\lvert i-j\rvert,\,N-\lvert i-j\rvert)$（边长 $a=1$）
- $D$ 从随机厄米起步，不预设图
- 强边集合贴近环近邻；非近邻边平均 $|z|$ 明显更小
- 通过：fraction of seeds reaching ring ≥ 0.25，且 best-S 种子满足 geo≤0.1、hit≥0.99、|z| ratio≤0.2

```bash
python -m experiments.exp1_ring_neighbor --n 6 --lam 30 --seeds 8
```

## Exp2 criterion (failed on $N=6$)

- Protocol A：虚相位图案的 $\mathrm{Tr}(D^{4})$ 应严格低于全实
- Protocol B：warm-start 后自由相位应比全实更低 $\mathrm{Tr}(D^{4})/S$
- **实测**：等模环上 5 种图案 $\mathrm{Tr}(D^{4})$ 全是 $36.0$；自由相位不优于全实 → **fail**

## Exp3 criterion (passed on $N_{t}=N_{x}=8$)

$$
D=\gamma^{0}\otimes D_{t}+\gamma^{1}\otimes D_{x}
$$

- $\gamma^{0}$ 反对称（$(\gamma^{0})^{2}=-1$）、$\gamma^{1}$ 对称（$(\gamma^{1})^{2}=+1$），交叉项因反对易相消 $\Rightarrow$ $D^{2}=-D_{t}^{2}+D_{x}^{2}$
- 判据：欧氏（$\gamma$ 全对称）$D^{2}$ 全非负；洛伦兹（$\gamma^{0}$ 反对称）$D^{2}$ 不定，号差 $(-,+)$，且有零本征值（光锥）

## Exp4 criterion (passed on $N=400$, $\xi=10$)

- 1D Dirac 算子 $D=$ 动能（$-i\sigma_{x}\partial_{x}$ 虚反对称）+ 质量 $m(x)\sigma_{z}$，畴壁 $m(x)=m_{0}\tanh((x-x_{0})/\xi)$
- 判据：均匀质量（绕数 0）→ 0 零模、有能隙；畴壁（绕数 1）→ 2 零模局域在畴壁（$E\approx 10^{-18}$）
- 零模数 $=$ 绕数 $\times 2$（朴素格点 Dirac 的费米子加倍）；这是 Jackiw–Rebbi / 指标定理的最小版：**拓扑荷 → 零模（= 拓扑缺陷 → 粒子）**

## Exp5 criterion (passed): 2D 陈数 → 零模（量子霍尔）

- 2D 格点 Dirac + 均匀磁场（陈数 $C=8$），谱 = Landau 能级
- 判据：零模数 $=$ 陈数 $\times 2$（费米子加倍）；陈数 8 → 16 零模
- 与 Exp4 同一条「拓扑阶梯」的 2D 级（1D 绕数 / 2D 陈数 = 单个 Dirac 算子的指标）；3D Hopf 荷不是指标、是 odd Chern–Simons 绕数（见 Exp6a）

---

## 理论缺口：号差为什么读不出

厄米 $D$ 的 $D^{2}$ 本征值恒 $\ge 0$（就是 $\lambda(D)^{2}$）。所以用 $\langle D^{2}\rangle$ 的正负划分类空/类时，**不能**从有限维 $D^{2}$ 的谱符号读出洛伦兹号差。

更深一层：相位（实/虚）**不改变** $D^{2}$ 的谱——$D^{2}$ 的本征值是 $|\lambda(D)|^{2}$，对相位盲。所以「虚耦合给号差」在厄米框架下**数学上不可能**；号差需要洛伦兹 $\gamma$ 结构（$\gamma_{0}^{2}=-1$）或非厄米 $D$，不是「相位」能给的。号差若存在，只能出现在连续极限的微分算子解释里。

### $2\times 2$ 钉死：号差 = 有向性，不是虚相位

一条耦合（$2\times 2$ 块）

$$
D=\begin{pmatrix}0 & z\\ w & 0\end{pmatrix},\qquad D^{2}=zw\cdot I
$$

本征值就是 $zw$。三种写法：

| 耦合 | 矩阵元 | 厄米 | D 本征值 | zw | D² | 号差 |
|------|--------|------|----------|-----|-----|------|
| 对称（空间无向） | D₁₂=D₂₁=a | 厄米 | ±a | a²>0 | 正定 | + |
| 虚耦合（时间厄米） | D₁₂=−D₂₁*=ia | 厄米 | ±a | a²>0 | 正定 | 无 |
| 反对称（时间有向） | D₁₂=−D₂₁=a | **非厄米** | ±ia | −a²<0 | **负定** | − |

厄米 = $w=z^{*}$，于是 $zw=|z|^{2}\ge 0$：无论 $z$ 实还是虚，$D^{2}$ 恒正定。唯一让 $D^{2}$ 变负的是 $w=-z$（实、反对称、不取共轭），此时 $D$ 非厄米（$D^{\dagger}=-D$）。

**所以：号差 = 时间方向的「有向性」（实反对称/非厄米），不是「虚相位」。** 空间对称给 $+$，时间反对称给 $-$，合起来 $D^{2}\to +\partial_{x}^{2}-\partial_{t}^{2}=\square$，洛伦兹号差 $(-,+)$ 自然出现。

---

## Exp2 的解析定理（$\mathrm{Tr}(D^{4})$ 对相位盲）

带相位环是标准 A–B 能谱：本征值 $\lambda_{k}=2\cos(2\pi k/n+\Phi/n)$，$\Phi=\sum\phi_{i}$ 为总磁通。于是

$$
\mathrm{Tr}(D^{4})=\sum_{k}\lambda_{k}^{4}=16\sum_{k}\cos^{4}\theta_{k}
$$

由 $\cos^{4}\theta=\frac{1}{8}(3+4\cos 2\theta+\cos 4\theta)$，且 $n\nmid 4$ 时 $\sum_{k}\cos 2\theta_{k}=\sum_{k}\cos 4\theta_{k}=0$，得

$$
\mathrm{Tr}(D^{4})=16\cdot\frac{1}{8}\cdot 3n=6n\qquad (=36,\ n=6)
$$

**与所有环边相位分布无关**。这解释了 Protocol A 中 5 种图案 $\mathrm{Tr}(D^{4})$ 全是 $36.0$，也说明 $\mathrm{Tr}(D^{4})$ **数学上不可能**通过偏好虚相位来驱动时间方向。

---

## 闭合环磁通项（Exp1 v2 / v3 / v4）

把「相位 = 规范联络」落成第一个对相位敏感的作用量项。动机：Exp2 已证谱矩 $\mathrm{Tr}(D^{2m})$ 对相位盲，故改走闭合环磁通

$$
\Phi_{ijk}=\phi_{ij}+\phi_{jk}-\phi_{ik}
$$

（格点规范理论的威尔逊环），用复数乘积算 $\sin\Phi$ 避开 `arg` 分支切割。

**Exp1 v2**（`exp1_v2_flux.py`）：作用量加 $\mu\sum_{i \lt j \lt k}\sin^{2}\Phi_{ijk}$。结果：flux 惩罚把磁通从随机（$\sim 9.4$）压到 $0$，但相位**没有**变实矩阵——它压到了「$0$ 或 $\pi$ 磁通」（$\mathbb{Z}_{2}$ 通量态）。铁证：$\mathrm{mean}(\cos\Phi)=1-2\cdot\mathrm{pi\_frac}$ 精确成立（4 个种子全对），说明每个三角形磁通是二值的（$0/\pi$）。

**Exp1 v3**（`exp1_v3_flux_compare.py`）：对比三种曲率定义 $F$，结论是「曲率项的定义决定磁通有序化到哪个零点」：

| 曲率 F | 零点 | 收敛到的磁通 |
|----------|------|--------------|
| sin Φ | 0 和 π | Z₂ 通量（0/π 简并） |
| abs(e^{iΦ}−1) | 只有 0 | 零磁通（实矩阵） |
| cos Φ | ±π/2 | 全三角形 ±π/2 最大磁通（`quarter_frac`=1.0） |

**Exp1 v4**（`exp1_v4_mcmc.py`）：Metropolis–Hastings 采样 $e^{-\beta S}$（**经典热**，不是量子）。结果：温度驱动磁通有序化（flux 从高温 $\sim 9.5$ 随机 → 低温 $\sim 0.5$ 冻结到 $0/\pi$），高温 `pi_frac`$\approx 1/3$ 精确等于均匀磁通下 $P(\cos\Phi\lt -0.5)$ 的概率；但 $\mathbb{Z}_{2}$ 简并**未被热涨落打破**（低温跨种子 `pi_frac`$\approx 0.5$，无偏好）。已知坑：$\beta=10$ 步长 $0.3$ 太大、accept$\approx 0$（未充分热化）；判据曾把 $0.5$ 当随机基准（实为 $1/3$），主结论不变。

**核心结论（负结果，两个钉子）**：

1. 经典曲率惩罚**不能生成**非平凡拓扑——它只把磁通抹平到 $F$ 的零点；哪个磁通被允许由 $F$ 的定义决定，不由物理偏好决定。
2. 经典热涨落**不能打破** $\mathbb{Z}_{2}$ 简并——打破 $0/\pi$ 简并 → SU(2)，需要量子相干（$e^{iS}$）或额外对称破缺，不是热。

> 这排除了「设计经典作用量让 SU(2)/Hopf 拓扑自发涌现」这条路，把「为什么 SU(2)」的答案从「作用量」推向「拓扑荷守恒/量子化」。

```bash
python -m experiments.exp1_v2_flux
python -m experiments.exp1_v3_flux_compare
python -m experiments.exp1_v4_mcmc
```

---

## 量子化自指 + 空间涌现

**量子化自指（断裂+自指→量子性）**——本体论内核，六步推导 + 三个数值实验（`exp_quantization_selfref.py`）：

- ExpA：位置 $X=\sigma_{3}$（对角）、动量 $P=\sigma_{1}$（非对角），$[X,P]=2i\sigma_{2}\neq 0$（自指两端必然不对易）
- ExpB：断裂贡献 $|d|^{2}$ 离散化，最小非零单位 $=\varepsilon^{2}$（断裂最小单位 = 一个量子 $=\hbar$）
- ExpC：投影观测 → 坍缩（叠加态重叠 $0.346$ → 本征态 $1.000$）

**空间涌现第一步（$D$ 自发长三维）**——Hopfion 的前提（Hopf 荷需要三维）：

- 暖启动（`exp_space_3d.py --warm-start`）：hit=$1.000$、geo=$0.001$ → **三维近邻是 $D$ 的稳定解**（已钉死）
- 纯随机（`exp_space_3d.py`）：hit=$0.62$，卡局部极小
- 退火（`exp_space_3d_anneal.py`）：hit $0.62\to 0.82$–$0.90$，自发**接近**三维，但精确难

**迹反常（`exp_trace_anomaly.py`）**——负结果：0 维 $D$ 模型没有量纲嬗变（迹反常需要时空维度+能标+RG 流，单个矩阵是 0 维）。这排除了「在 0 维 $D$ 里算迹反常」这条路，指向「尺度从反常来」需要先有空间（= 空间涌现）。

**核心结论**：

- 量子性 = 断裂 + 自指 = 自指观测，已有推导+数值锚点（本体论内核立住）
- 空间涌现：复现三维已钉死（暖启动 hit=$1.0$）；档 2 维度对比（负结果）证「复现几何不选三维、偏向低维」；谱维数读步已落地（$d_{s}=2/3/4$ 可读出）

```bash
python -m experiments.exp_quantization_selfref
python -m experiments.exp_space_3d --warm-start
python -m experiments.exp_space_3d_anneal
python -m experiments.exp_trace_anomaly
```

---

## Exp6a：Hopf 荷 = odd Chern–Simons

`experiments/exp6a_hopf_charge.py`——把「拓扑荷 = 谱不变量」从 1D/2D（绕数/陈数 = 零模，指标定理）推进到 3D（Hopf 荷 = odd Chern–Simons 绕数）。

方法：给定 $n(x):\mathbb{R}^{3}\to S^{2}$，算 emergent 磁场 $B_{i}=n\cdot(\partial_{j}n\times\partial_{k}n)$（谱微分），FFT 库仑规范解 $\nabla\times A=B$（$\hat{A}=-i(k\times\hat{B})/|k|^{2}$），

$$
Q_{H}=\frac{1}{16\pi^{2}}\int A\cdot B\,d^{3}x
$$

| 场 | 预期 | 数值 Q_H |
|----|------|--------------|
| 平凡场 n=(0,0,1) | 0 | 0（精确） |
| 标准 Hopf 映射（Q=1） | 1 | 0.9962（N=64,L=5）→ **0.99996**（N=128,L=10） |
| Q=2 场 | 4 | 3.99–4.00 |

- $Q=2\to 4$ 不是 bug，是 **Hopf 不变量的复合律** $H(g\circ f)=(\deg g)^{2}\,H(f)$：$n=\mathrm{inv.stereo}((z_{1}/z_{2})^{Q})$ 在目标空间复合了度 $Q$ 映射。这**独立验证了算法能分辨不同拓扑荷**（$1$ vs $4$）。
- Hopf 荷对盒子敏感：skyrmion 密度在无穷远 $1/r^{4}$ 衰减，$L=5$ 截断尾部（$0.996$），$L=10$ 才收敛。
- 诚实边界：Exp6a 是**预设三维格点**，第三层（$D$ 自发长出 $n(x)$）仍开放；谱流版（族指标/Bott）是 Exp6b。

```bash
python -m experiments.exp6a_hopf_charge
```

---

## 谱流前置

旧版用「均匀质量 $m$ 从 $-m_{0}\to +m_{0}$」算谱流得 flow=$0$——物理上对（1D 周期 Dirac 能级成对只触碰 $0$ 不穿越）。修正为**畴壁 + $\lambda\sigma_{y}$ + 能级追踪**：畴壁束缚精确零模（$E\approx 10^{-16}$），$\lambda\sigma_{y}$ 给零模能量偏移，零模随 $\lambda$ 穿越 $0$ 被能级追踪（贪心最近邻）捕获。「净谱流=绕数」作为整数不变量需闭环参数族 + 手征破缺，留 Exp6b。

```bash
python -m experiments.exp_spectral_flow
```

## 空间涌现·维度对比

`experiments/exp_space_dim_compare.py`——把目标距离 $L$ 换成 2D/3D/4D 超立方曼哈顿距离，比较 $D$ 复现各维度的难易（同一框架 $S=\mathrm{Tr}(D^{2})+\lambda\cdot\mathrm{geo}$，只优化模长、相位固定 $0$，因作用量只依赖 $|z|$）。

| 维度 | N | 暖启动 hit | 冷启动 excess（归一化） |
|------|-----|------------|-------------------------|
| 2D | 25 | 1.0（稳定） | 0.64 |
| 3D | 27 | 1.0（稳定） | 0.59 |
| 4D | 16 | 1.0（稳定） | 0.55（baseline 修正后） |

**负结果**：① 暖启动 2D/3D/4D 全 hit=$1.0$——「三维是稳定解」不特殊；② 冷启动 2D > 3D > 4D——维度越低越好找（配位数越少 $\mathrm{Tr}(D^{2})$ 越便宜）。**结论：复现几何框架不选三维、反而偏向低维；「三维特殊」必须来自 Hopf 拓扑压力项（第三层），不在复现几何这一层。**

坑：4D $n=2$ 超立方全角点、无内部态，真实配位是 $4$ 不是 $2d=8$；baseline 应按 `mean(len(neighbors))` 算。

```bash
python -m experiments.exp_space_dim_compare
python -m experiments.exp_cold_start
```

## $3/2$ = 谱流 + $\eta$ 分数部分

`experiments/exp_eta_framing.py`——1D 数值演示「总拓扑荷 = 整数 Hopf（谱流）+ 分数自旋（$\eta$ framing）= 半整数 $3/2$」。

模型：$S^{1}$ 上 1D Dirac 算子 $D_{A}=-i\,d/d\theta+A$（周期边界），本征值 $\lambda_{n}=n+A$（$n$ 整数，精确）。

- **谱流（整数）**：一族 $D(t)=-i\,d/d\theta+t$，$t:0\to 1$ 的本征值净上穿 $0$ 次数 = **$1$**（数值逐点追踪）
- **$\eta$-不变量（分数）**：热核正则化 $\eta_{\varepsilon}=\sum_{n}\mathrm{sign}(n+A)\,e^{-\varepsilon|n+A|}$，$\varepsilon\to 0$ 收敛到闭式 $1-2A$——$A=1/4\to +1/2$、$A=3/4\to -1/2$（自旋 framing 的 $\pm 1/2$，数值偏差降到 $4\times 10^{-7}$）
- **合成**：$1+(1/2)=3/2$

诚实边界：这是 1D 类比（演示机制）；真正的 3D Hopf 荷（整数=$1$）已由 Exp6a 钉死（$0.99996$），这里补的是 $\eta$ 的分数部分（$1/2$）。真正的 $S^{3}$ 上自旋 $1/2$ Hopfion 的 $\eta=3/2$ 需在具体 Dirac 算子上把 $\eta$ 完整算出来（未做）。

```bash
python -m experiments.exp_eta_framing
```

## 谱维数区分维度

`experiments/exp_spectral_dim_compare.py`——谱层第一步「读」：谱维数 $d_{s}$ 从热核 $K(t)=\mathrm{Tr}(e^{-tL})\sim t^{-d_{s}/2}$ 读出，区分 2D/3D/4D，不预设坐标。

用周期边界 $d$ 维超立方晶格的拉普拉斯 $L$ 的解析本征值 $\lambda=4\sum_{j}\sin^{2}(\pi m_{j}/n)$，热核迹

$$
K(t)=\Bigl[\sum_{m}e^{-4t\sin^{2}(\pi m/n)}\Bigr]^{d}
$$

在中间 $t$ 区间（$1\ll t\ll n^{2}$）拟合 $\log K$ vs $\log t$，斜率 $=-d_{s}/2$。

| 维度 | n | t 区间 | d_s | 目标 |
|------|-----|----------|---------|------|
| 2D | 100 | [10,100] | 2.009 | 2 ✓ |
| 3D | 60 | [10,100] | 3.014 | 3 ✓ |
| 4D | 40 | [10,60] | 4.018 | 4 ✓ |

**两个坑（重要）**：① $K(t)$ 里的算子必须是**图拉普拉斯**，不是邻接矩阵——邻接矩阵的谱在能带中心 DOS 近常数，所有维度都读出 $d_{s}\approx 1$（错）；② $t$ 要取**中间区**（$1\ll t\ll n^{2}$），不是「小 $t$」——小 $t$ 时 $K\approx N$ 常数（斜率 $\approx 0$）。

**诚实定位**：这是「读」步（谱维数读出 $D$ 的维数），不是「选」步。谱维数本身不选三维，只是给了一个**不预设坐标的维度读数**；「选」步（偏好 $d_{s}=3$）预设了三维，是循环，落在第三层（为什么 SU(2)）。

```bash
python -m experiments.exp_spectral_dim_compare
```

---

## 弦网凝聚 sanity check

`experiments/exp_string_net_condensation.py`——把「自指闭弦 → $\mathrm{SU}(2)_{k}$ 弦网相」这个 if-then 的两个「if」**分开**做 sanity check，明确不组装成证明。理论路线是：自指闭弦 → 弦网相 → 谱三元组 → GR，其中「闭弦 → 弦网相」需要两个 if：① 张力 $T\to 0$（弦凝聚）；② 融合结构 $=\mathrm{SU}(2)_{k}$（=「为什么 SU(2)」公设）。

**Part A：$\mathrm{SU}(2)_{k}$ 拓扑纠缠熵目标（解析）**——量子维数 $d_{j}=[2j+1]_{q}$、总量子维数 $\mathcal{D}=\sum_{j}d_{j}^{2}$、TEE $\gamma=\log\mathcal{D}$：

| k | 量子维数 d_j | D | TEE = log D |
|-----|------------------|---------------|------------------------|
| 1 | 1, 1 | 2 | 0.6931 |
| 2 | 1, √2, 1 | 4 | 1.3863 |
| 3 | 1, φ, φ, 1 | 7.2361 | 1.9791 |
| 4 | 1, √3, 2, √3, 1 | 12 | 2.4849 |
| 5 | 1, 1.802, 2.247, 2.247, 1.802, 1 | 18.5918 | 2.9227 |

（$k=4$ 时 $\mathcal{D}=12$、$\log 12=2.4849$ 恰好撞「$\ln 12$」，是**巧合**，勿追。）

**Part B：toy 环气体凝聚（张力 $T\to 0$）**——$3\times 3$ 环面、$1024$ 个闭合环构型，$H=T\cdot(\mathrm{length})-\lambda\sum_{p}$（面翻转）。扫 $T$ 大→小：环密度 $0.002\to 0.500$（单调），参与率 $1.02\to 460$（局域→离域）。**这是 if-1「张力→0 触发弦凝聚」的机制签名**，但它是**等权**（单一环型，$\mathcal{D}=2$，toric-code 型），不是 $\mathrm{SU}(2)_{k}$ 加权凝聚。

**Part C：加权弦网凝聚（$B_{p}$ 算符，单面，一步插入）**——弦网动能项 $B_{p}=(1/\mathcal{D})\sum_{j}d_{j}B_{p}^{j}$ 的基态按量子维数加权，对比等权 $B_{p}=(1/N)\sum_{j}B_{p}^{j}$。单面环型分布：

| k | 加权 p_j=d_j²/D | 等权 p_j=1/N |
|-----|------------------------------------|------------------|
| 2 | 0.250, 0.500, 0.250 | 0.333×3 |
| 3 | 0.138, 0.362, 0.362, 0.138 | 0.250×4 |
| 4 | 0.083, 0.250, 0.333, 0.250, 0.083 | 0.200×5 |

**加权凝聚把真空压到 $p_{0}=1/\mathcal{D}$、抬高高量子维数自旋**——$d_{j}^{2}/\mathcal{D}$ 分布正是弦网凝聚的指纹，也是 $\log\mathcal{D}$（TEE）的来源。这是往 Levin–Wen 靠的第一步，但仍是单面/一步插入，**没有 F-符号、没有顶点融合规则**。

**定位与边界（诚实）**：Part A–C 都在**已知理论（Levin–Wen 弦网）内部**复现公式，**没有让自指 $D$ 自己去选 SU(2)**，因此**没有推进「自指闭弦 → $\mathrm{SU}(2)_{k}$ 弦网相」**。该路线**封存**（不再加 F-符号、不再扩 Levin–Wen），本实验仅作为「已知工具复现」的存档。

```bash
python -m experiments.exp_string_net_condensation
```

---

## Born 偏离判据演示

`experiments/exp_born_deviation.py`——把「弱命题 + 偏离判据」落成可跑的量。

**弱命题**：自指（无外部观察者）= 内部 SO(3)（断裂→SU(2) 伴随作用）无偏好方向 $\Rightarrow$ 坍缩概率 SO(3) 协变 $\Rightarrow$ $p=f(m\cdot n)$（只依赖态与测量的夹角）。

| 部分 | 内容 | 结果 |
|------|------|------|
| Part A | p_Born=(1+m̂·n̂)/2 在 50 个随机 SO(3) 转动下 | 偏离 6e-16（=0，SO(3) 协变） |
| Part B | 引入偏好方向 ê：p_corr=p_Born+0.3(m̂·ê) | 偏离 0.50（非零；实测=解析） |
| Part C | 偏离随 ε 扫描 | 线性缩放，ε=0 时归零 |

**判据**：转动整个系统（态 + 测量轴一起转），概率不变 = 自指；概率变 = 有绝对方向 $\hat{e}$ 混入，偏离度 = 外部性的度量。

**诚实边界**：演示的是弱命题 + 偏离判据，**不是** Born 规则的推导——$p_{\mathrm{Born}}=(1+x)/2$ 的形式是**作为输入**的。所以不撞 Gleason/Zurek。

```bash
python -m experiments.exp_born_deviation
```

---

## 纠缠熵·面积律 vs 体积律

`experiments/exp_entanglement_area_law.py`——纠缠线的第一个 toy：演示「纠缠熵的标度 = 局域性签名」。

约化密度矩阵 $\rho_{A}=\mathrm{Tr}_{B}|\psi\rangle\langle\psi|$，纠缠熵 $S=-\mathrm{Tr}(\rho_{A}\ln\rho_{A})$；自由费米子用 Peschel 公式（关联矩阵本征值 $\lambda$）：

$$
S=-\sum_{\lambda}\bigl[\lambda\ln\lambda+(1-\lambda)\ln(1-\lambda)\bigr]
$$

1D 链 $N=400$：

| L | S_local（费米海） | S_random |
|-----|--------------------------------|------------------------|
| 10 | 0.85 | 6.81 |
| 40 | 1.09 | 25.77 |
| 200 | 1.28 | 77.58 |

- 局域：$S\sim 0.15\ln L$（$\approx(1/6)\ln L$，边缘 $c=1$）——亚体积 / 局域性
- 随机：$S\sim 0.38\,L$——体积 / 非局域

诚实边界：运动学纠缠（图分割 + ansatz）；不是背景无关；不是「$D$ 自发做出面积律」。

### 面积律 + 维度

`experiments/exp_entanglement_area_law_dim.py`：

| 维度 | S 标度 | S/ℓ^(d−1) | 结论 |
|------|-----------|----------------|------|
| 2D（40×40） | S∼ℓ^1.25 | 0.83→1.29 ~常数 | 周长律（S∼ℓ），不是 ℓ² |
| 3D（12³） | S∼ℓ^2.25 | 0.76→1.00 ~常数 | 面积律（S∼ℓ²），不是体积 |

指数 $1.25/2.25=1,2+$ 对数尾巴（费米面 log 破坏）。结论：$S\sim L^{d-1}$；指数 = 维度——量子腿与几何腿对接。诚实边界：图是预设的；证明的是「若 $D$ 活在 $d$-正则格点基态上，指数 $=d-1$」，不是「$D$ 变成三维」。

### 面积律 = 局域性

`experiments/exp_area_law_locality.py`：

| 图 | S(L=10) | S(L=100) | 增长 | 律 |
|----|-----------|------------|------|----|
| 1D 链（费米面） | 0.85 | 1.17 | 1.38× | 面积（对数） |
| 随机 G(N,0.3) | 6.7 | 38.6 | 5.7× | 体积（线性） |

面积律 $\Leftrightarrow$ 费米面（平移不变/局域）；体积律 $\Leftrightarrow$ 无费米面。**重要诚实**：长程跳 $1/r^{\alpha}$ 仍是 Toeplitz → 仍面积律；纠缠局域性 ≠ 跳程 = 费米面是否存在。真非局域 = 破平移的随机图。

```bash
python -m experiments.exp_entanglement_area_law
python -m experiments.exp_entanglement_area_law_dim
python -m experiments.exp_area_law_locality
```

---

## 裸重连（负结果，B 路线收口）

`experiments/exp_bare_reconnection.py`——4-正则图 FPL + 局域重连；$N=30$，$10^{5}$ 步。

- 闭弦数 $2$–$10$ 来回跳；最大长度 $10$–$30$；均值 $5.6$–$30$；$n\times\mathrm{avg}=2N$ 成立
- **没有**「少数大弦主导」的趋势
- 机制：全接受 ⇒ 在 $3^{N}$ 配对上对称随机游走 ⇒ 均匀 ⇒ 无偏好

**裸动力学三连钉**：① 经典作用量不能生成拓扑；② 经典热不能打破 $\mathbb{Z}_{2}$；③ 裸重连不自组织 → 需要拓扑守恒（Hopf/linking = 第三层）。

```bash
python -m experiments.exp_bare_reconnection
```

---

## 相位当联络

`experiments/exp_phase_connection.py`——把 $e^{i\theta_{ij}}=D_{ij}/|D_{ij}|$ 当作离散 U(1) 联络。在 $3\times 3$（$4$ 面）上跑麦克斯韦 $S_{\mathrm{top}}=-\sum_{p}\cos\Phi_{p}$ → 全部磁通 $0\bmod 2\pi$，$S_{\mathrm{top}}$ 打到下界 $-N_{\mathrm{faces}}=-4$。

结论：相位当联络是对的；Abel U(1) 曲率 → 平坦/平凡；拓扑荷 = U(1) 平坦 + 非 Abel 的 3D Hopf/linking（后者才选三维）。

## Wilson 全局不变量

`experiments/exp_wilson_invariant.py`——$Q=\prod_{C}W(C)=\prod_{\mathrm{edges}}e^{i\theta}$ 在重连下严格不变（每条边恰好一次）。$Q_{\mathrm{wilson}}=Q_{\mathrm{direct}}$（arg $0.172388$）；$5\times 10^{4}$ 次重连后偏差 $4\times 10^{-15}$。

链：$D_{ij}=D_{ji}^{*}\Rightarrow U_{ij}\Rightarrow W(C)\Rightarrow Q$——纯图 $D$ 的第一个拓扑不变量（种子）；Hopf/linking 是下一步。

## SU(2) Wilson

`experiments/exp_su2_wilson.py`——非 Abel $W(C)=\mathrm{tr}\prod U_{ij}$。环数 $8\to 2\to 6\to 4$，迹会变（**不像** U(1) 的 $Q$ 那样不变）；Fierz/skein $\mathrm{tr}A\,\mathrm{tr}B=\mathrm{tr}(AB)+\mathrm{tr}(AB^{\dagger})$ 成立（误差 $2\times 10^{-16}$）。

为什么 SU(2) = 自同构（断裂→二元→SU(2)）= 相对链接的非 Abel 群——同一个 SU(2)，两张脸。诚实边界：真正的不变量是 Jones（不是单个实数）。

## 群交换子 $L$

`experiments/exp_group_commutator.py`：

$$
L=1-\tfrac{1}{2}\mathrm{Tr}\bigl(W_{1}W_{2}W_{1}^{-1}W_{2}^{-1}\bigr)
$$

SU(2)：$5000$ 随机对，均值 $L=0.748$，非零 $100\%$；U(1)：$L\equiv 0$（最大 $4\times 10^{-16}$）。边界：$L$ = 非对易度量，**不是**链接数；真 2D=$0$/3D$\neq 0$ 需要 Jones/Markov 迹。

## Yang–Baxter

`experiments/exp_yang_baxter.py`——QYBE $R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12}$ 与辫子 $B_{12}B_{23}B_{12}=B_{23}B_{12}B_{23}$（$B=PR$）对所有 $k$ 成立（误差 $2.5\times 10^{-16}$）；辫子本征值 $\{q,q,q,q^{2}\}=3\oplus 1$。诚实边界：完整 $R$ 是 $4\times 4$；$2\times 2$ 是对角化后的通道。

## Jones

`experiments/exp_jones.py`——Kauffman + writhe → Jones。对所有 $k$：$\langle\mathrm{trivial}\rangle\neq\langle\mathrm{Hopf}\rangle$，$V$ 可分。$k=1$（$A=q^{1/4}$，$q=e^{i\pi/3}$）：Kauffman $-1.732$ vs $-1.0$；Jones $-1.732$ vs $+i$。

诚实边界：Jones = Kauffman + $(-A)^{-3w}$；2D/3D 是嵌入区分；接到抽象 $D$ 需要嵌入/辫子词。

## 配对 ↔ 交叉

`experiments/exp_pairing_to_braid.py`——4-价节点配对：$(13)(24)=1$，$(12)(34)=e$，$(14)(23)=$交叉；skein $\sigma=A\cdot 1+A^{-1}e$。

结果：skein 一致当且仅当 $e^{2}=de$，$d=-A^{2}-A^{-2}$；$R$ 本征值 $3\oplus 1$ = 通道而非 $\sigma/\sigma^{-1}$；$\pm$ 交叉 = $A\leftrightarrow A^{-1}$；纯 FPL Kauffman $=d^{c-1}$。含义：经典纯配对无符号；符号需要量子相干。

## 辫子词 → Jones 多项式

`experiments/exp_braid_word_to_jones.py`——完整正向链：辫子 → skein → TL 配对 → Markov（union-find）→ Kauffman → writhe → Jones。

Hopf $\sigma_{1}^{2}$（$-A^{4}-A^{-4}$）、Trefoil $\sigma_{1}^{3}$、Figure-8 对上教科书；**R-II**、**R-III** 成立。

已知坑：`coeff=1j` 会把全体乘上 $i$；转置换圆模型在相邻 cap 上错（改用弦图 UF）；trefoil 手性 $\sigma_{1}^{3}$=左。正向已闭环；**反向**（配对 → 唯一辫子）不唯一——缺 over/under。缺口 = 经典配对如何获得振幅 $A$。

## $A$ 从哪来

`experiments/exp_quantization_condition.py`——经典：任意 $A\neq 0$ 都满足 $RR^{-1}=I$（$10^{-16}$）——$A$ 自由。量子化：单位根 $q^{k+2}=1$ 固定 $A=q^{1/4}$；$-d=A^{2}+A^{-2}=2\cos(\pi/(2(k+2)))$ = 量子维数（$k=1\to 1.732\lt 2$）；$k\to\infty$ 回到 $2$。开关 = 量子相干 / 单位根输入。

## Jones–Wenzl

`experiments/exp_jones_wenzl.py`——Hecke $\sigma^{2}=(q-1)\sigma+q$ 对任意 $q$ 都成立（**不**逼出单位根）。JW 幂等元 $f_{n}$ 存在当且仅当 Chebyshev $\Delta_{n-1}\neq 0$。$\mathrm{SU}(2)_{k}$ 取 $\delta=-2\cos(\pi/(k+2))$ $\Rightarrow$ $\Delta_{k+1}=0$ = 截断 = 单位根。

**量子化 = JW 截断的存在性**，不是 Hecke。

## 关系 → 拓扑

`experiments/exp_relation_to_topology.py`——厄米 $D$：谱投影 $P_{k}=q_{k}q_{k}^{\dagger}$ 自动幂等（$|P^{2}-P|\sim 10^{-16}$，本征值 $\{0,1\}$）——免费。两座桥：免费（Herm→投影）+ 付费（投影→单位根截断 $\Delta_{k+1}=0$）。第三层墙 = 付费桥：为什么 $\delta=2\cos(\pi/(k+2))$。

## 有限 $N$ → 单位根

`experiments/exp_finite_quantization.py`——Chebyshev 恒等式：$\delta=2\cos(\pi/(N+1))$ $\Rightarrow$ $\Delta_{N-1}=1$（$f_{N}$ 存在），$\Delta_{N}=0$（$f_{N+1}$ 消失）；$N=3..8$ 精确。这正是 $\delta=\mathrm{SU}(2)_{k}$ 且 $k=N-1$。$A_{k}$ 结构：自旋 $0,\tfrac{1}{2},\ldots,k/2$ 融合 $[\tfrac{1}{2}]^{2}=[0]+[1]$。

**「有限 $N$ 最大股数 ⇒ $\delta=$ 单位根」是定理**（Chebyshev）。仍开放的输入：「$D$ 用满 $N$ 股」= v3 公设「观测 = 一次有向区分」。

```bash
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
python -m experiments.exp_spectral_flow
python -m experiments.exp6a_hopf_charge
python -m experiments.exp_space_dim_compare
python -m experiments.exp_eta_framing
python -m experiments.exp_spectral_dim_compare
python -m experiments.exp_cold_start
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

结果落在 `experiments/*_last_run.json`。依赖：`numpy>=1.24`，`scipy>=1.10`（见 `requirements.txt`）。

---

## Repository layout

```
self_ref_spacetime/
  src/
    algebra.py      # 厄米 D：模长 + 相位
    distance.py     # 边长 1/|z| → 最短路 d_ij
    action.py       # Tr(D²) + λ·geo + λ₄ Tr(D⁴)
    metrics.py      # 近邻命中、相位、谱诊断
    optimize.py     # L-BFGS-B；free / real / phases_only
    flux.py         # 闭合环磁通（Exp1 v2+）
  experiments/      # 全部门控脚本
  tests/            # 环距离解析一致性
  exp7_design.md    # Exp7 聚类 / SU(2) 涌现设计稿（未实现）
  requirements.txt
  README.md
```

---

## Scope boundary

| 主张 | 立场 |
|------|------|
| 环近邻可从靶距离涌现 | ✅ 数值证据（多峰） |
| 虚相位给出洛伦兹号差 | ❌ 已证伪（Exp2 + 2×2） |
| Dirac / 有向结构给出号差 | ✅ Exp3 |
| 拓扑荷 → 零模 | ✅ Exp4 / Exp5 |
| Hopf 荷数值可读 | ✅ Exp6a（预设 3D 格点） |
| D 自发选三维 | ❌ 复现几何负结果 |
| 经典作用量 / 热 / 裸重连生成拓扑 | ❌ 三连钉 |
| 自指 D → SU(2)_k 弦网 | ⛔ 路线封存（sanity = 已知工具） |
| 证明 3+1 / 纯谱选择 | out of scope |

---

## License

研究原型代码；默认按仓库约定私有保留。若对外发布，请另加 `LICENSE` 文件。
