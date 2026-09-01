# Self-Referential Spacetime — Numerical Closure Tests

把「自指时空」里可证伪的两步做成数值实验。挂在 `AI_Analysis/self_ref_spacetime`，与业务代码隔离。

**v0.1 命题（环）**：固定环型目标距离 \(L_{ij}\) 时，
\(S=\mathrm{Tr}(D^2)+\lambda\sum_{i<j}(d_{ij}-L_{ij})^2\)
的数值极小解，是否自发变成近邻环图。

不声称证明 3+1。

## Setup

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
```

## Status

| Item | Result |
|------|--------|
| Exp1 环 + 近邻涌现 | **pass**（多峰；约 3/8 种子进干净环盆地） |
| Exp2 \(\mathrm{Tr}(D^4)\) 虚耦合 | **fail**（等模环上相位不降 TrD4；自由相位不优于全实） |
| Exp3 Dirac 结构 → 洛伦兹号差 | **pass**（欧氏 D² 全非负；洛伦兹 D² 正负 48/48/32，不定 \(-+\)） |
| Exp4 拓扑荷 → 零模（Jackiw-Rebbi） | **pass**（绕数 0 无零模；绕数 1 → 2 零模在畴壁，E≈1e-18） |
| Exp5 二维陈数 → 零模（量子霍尔） | **pass**（陈数 8 → 16 零模 = 陈数 × 2 费米子加倍） |
| Exp1 v2 闭合环磁通项 | **pass**（flux 惩罚把磁通从随机压到 0/π Z₂ 通量，0/π 简并） |
| Exp1 v3 曲率项定义对比 | **pass**（sin→Z₂、wilson→零磁通、cos→±π/2；曲率定义决定磁通零点） |
| Exp1 v4 热采样（MCMC） | **负结果**（温度驱动磁通有序化，但经典热不打破 Z₂ 简并） |
| 迹反常（exp_trace_anomaly） | **负结果**（0 维 D 模型无量纲嬗变——迹反常需时空维度+能标，单个矩阵是 0 维） |
| 量子化自指三实验（exp_quantization_selfref） | **pass**（[X,D]=2iσ₂ 不对易、断裂最小单位=ε²、投影坍缩） |
| 空间涌现·暖启动（exp_space_3d --warm-start） | **pass**（hit=1.000、geo=0.001，三维近邻是稳定解） |
| 空间涌现·纯随机（exp_space_3d） | 卡局部极小（hit 0.62） |
| 空间涌现·退火（exp_space_3d_anneal） | 自发接近三维（hit 0.82–0.90），精确难 |
| 纯谱 / 3+1 选择 | out of scope |

## Pass criterion (Exp1)

- \(N=6\) 或 \(8\)；\(L_{ij}=\min(|i-j|,N-|i-j|)\)（边长 \(a=1\)）
- \(D\) 从随机厄米起步，不预设图
- 强边集合贴近环近邻；非近邻边平均 \(|z|\) 明显更小

## Exp2 criterion (failed on N=6)

- Protocol A：虚相位图案的 \(\mathrm{Tr}(D^4)\) 应严格低于全实
- Protocol B：warm-start 后自由相位应比全实更低 \(\mathrm{Tr}(D^4)/S\)

## Exp3 criterion (passed on Nt=Nx=8)

- \(D=\gamma^0\otimes D_t+\gamma^1\otimes D_x\)，\(\gamma^0\) 反对称（\((\gamma^0)^2=-1\)）、\(\gamma^1\) 对称（\((\gamma^1)^2=+1\)），交叉项因反对易相消 \(\Rightarrow D^2=-D_t^2+D_x^2\)
- 判据：欧氏（γ 全对称）\(D^2\) 全非负；洛伦兹（γ⁰ 反对称）\(D^2\) 正负对称（不定，号差 \(-+\)），且有零本征值（光锥）

## Exp4 criterion (passed on N=400, xi=10)

- 1D Dirac 算子 \(D=\) 动能（\(-i\sigma_x\partial_x\) 虚反对称）+ 质量 \(m(x)\sigma_z\)，畴壁 \(m(x)=m_0\tanh((x-x_0)/\xi)\)
- 判据：均匀质量（绕数 0）→ 0 零模、有能隙；畴壁（绕数 1）→ 2 零模局域在畴壁（E≈1e-18）
- 零模数 = 绕数 × 2（朴素格点 Dirac 的费米子加倍）；这是 Jackiw-Rebbi / 指标定理的最小版：**拓扑荷 → 零模（= 拓扑缺陷 → 粒子）**

## 理论缺口（号差为什么读不出）

厄米 \(D\) 的 \(D^2\) 本征值恒 \(\ge 0\)（就是 \(\lambda(D)^2\)）。所以用 \(\langle D^2\rangle\) 的正负划分类空/类时，**不能**从有限维 \(D^2\) 的谱符号读出洛伦兹号差。

更深一层：相位（实/虚）**不改变** \(D^2\) 的谱——\(D^2\) 的本征值是 \(|\lambda(D)|^2\)，对相位盲。所以"虚耦合给号差"在厄米框架下**数学上不可能**；号差需要洛伦兹 \(\gamma\) 结构（\(\gamma_0^2=-1\)）或非厄米 \(D\)，不是"相位"能给的。号差若存在，只能出现在连续极限的微分算子解释里。

### 2×2 钉死：号差 = 有向性，不是虚相位

一条耦合 \(D=\begin{pmatrix}0&z\\w&0\end{pmatrix}\)，平方 \(D^2=zw\cdot I\)，本征值就是 \(zw\)。三种写法：

| 耦合 | \(D\) | 厄米 | \(D\) 本征值 | \(zw\) | \(D^2\) | 号差 |
|---|---|---|---|---|---|---|
| 对称（空间无向） | \(\begin{pmatrix}0&a\\a&0\end{pmatrix}\) | 厄米 | \(\pm a\) | \(a^2>0\) | 正定 | \(+\) |
| 虚耦合（时间厄米） | \(\begin{pmatrix}0&ia\\-ia&0\end{pmatrix}\) | 厄米 | \(\pm a\) | \(a^2>0\) | 正定 | 无 |
| 反对称（时间有向） | \(\begin{pmatrix}0&a\\-a&0\end{pmatrix}\) | **非厄米** | \(\pm ia\) | \(-a^2<0\) | **负定** | \(-\) |

厄米 \(=\) \(w=\overline z\)，于是 \(zw=|z|^2\ge0\)：无论 \(z\) 实还是虚，\(D^2\) 恒正定。唯一让 \(D^2\) 变负的是 \(w=-z\)（实、反对称、不取共轭），此时 \(D\) 非厄米（\(D^\dagger=-D\)）。

**所以：号差 = 时间方向的"有向性"（实反对称/非厄米），不是"虚相位"。** 空间对称给 \(+\)，时间反对称给 \(-\)，合起来 \(D^2\to+\partial_x^2-\partial_t^2=\square\)，洛伦兹号差 \((-,+)\) 自然出现。

## Exp2 的解析定理（Tr(D⁴) 对相位盲）

带相位环是标准 A-B 能谱：本征值 \(\lambda_k=2\cos(2\pi k/n+\Phi/n)\)，\(\Phi=\sum\phi_i\) 为总磁通。于是

$$\mathrm{Tr}(D^4)=\sum_k\lambda_k^4=16\sum_k\cos^4\theta_k$$

由 \(\cos^4\theta=\frac18(3+4\cos2\theta+\cos4\theta)\)，且 \(n\nmid4\) 时 \(\sum_k\cos2\theta_k=\sum_k\cos4\theta_k=0\)，得

$$\mathrm{Tr}(D^4)=16\cdot\tfrac18\cdot 3n=6n\quad(=36,\ n=6)$$

**与所有环边相位分布无关**。这解释了 Protocol A 中 5 种图案 Tr(D⁴) 全是 36.0，也说明 Tr(D⁴) **数学上不可能**通过偏好虚相位来驱动时间方向。

## 闭合环磁通项（Exp1 v2 / v3 / v4，2026-09-01）

把「相位 = 规范联络」落成第一个对相位敏感的作用量项。动机：Exp2 已证谱矩 \(\mathrm{Tr}(D^{2m})\) 对相位盲，故改走闭合环磁通 \(\Phi_{ijk}=\phi_{ij}+\phi_{jk}-\phi_{ik}\)（格点规范理论的威尔逊环），用复数乘积算 \(\sin\Phi\) 避开 arg 分支切割。

**Exp1 v2**（`exp1_v2_flux.py`）：作用量加 \(\mu\sum_{i<j<k}\sin^2\Phi_{ijk}\)。结果：flux 惩罚把磁通从随机（~9.4）压到 0，但相位**没有**变实矩阵——它压到了「0 或 π 磁通」（Z₂ 通量态）。铁证：\(\text{mean\_cos}\Phi = 1-2\times\text{pi\_frac}\) 精确成立（4 个种子全对），说明每个三角形磁通是二值的（0/π）。

**Exp1 v3**（`exp1_v3_flux_compare.py`）：对比三种曲率定义 \(F\)，结论是「曲率项的定义决定磁通有序化到哪个零点」：

| 曲率 F | 零点 | 收敛到的磁通 |
|---|---|---|
| \(\sin\Phi\) | 0 和 π | Z₂ 通量（0/π 简并） |
| \(|e^{i\Phi}-1|\) | 只有 0 | 零磁通（实矩阵） |
| \(\cos\Phi\) | ±π/2 | 全三角形 ±π/2 最大磁通（quarter_frac=1.0） |

**Exp1 v4**（`exp1_v4_mcmc.py`）：Metropolis-Hastings 采样 \(e^{-\beta S}\)（**经典热**，不是量子）。结果：温度驱动磁通有序化（flux 从高温 ~9.5 随机 → 低温 ~0.5 冻结到 0/π），高温 pi_frac≈1/3 精确等于均匀磁通下 \(\cos\Phi<-0.5\) 的概率；但 Z₂ 简并**未被热涨落打破**（低温跨种子 pi_frac≈0.5，无偏好）。已知坑：β=10 步长 0.3 太大、accept≈0（未充分热化）；判据曾把 0.5 当随机基准（实为 1/3），结论不受影响。

**核心结论（负结果，两个钉子）**：
1. 经典曲率惩罚**不能生成**非平凡拓扑——它只把磁通抹平到 F 的零点；哪个磁通被允许由 F 的定义决定，不由物理偏好决定。
2. 经典热涨落**不能打破** Z₂ 简并——打破 0/π 简并 → SU(2)，需要量子相干（\(e^{iS}\)）或额外对称破缺，不是热。

> 这排除了「设计经典作用量让 SU(2)/Hopf 拓扑自发涌现」这条路，把「为什么 SU(2)」的答案从「作用量」推向「拓扑荷守恒/量子化」（= B4 迹反常 + 弦网凝聚）。已同步记入 vault 的《涌现时空模型：后续可证明与扩展方向》（排除的错路）与《字典观测点》（记录表）。

## 量子化自指 + 空间涌现（exp_quantization_selfref / exp_space_3d / exp_space_3d_anneal，2026-09-01）

**量子化自指（断裂+自指→量子性）**——本体论内核，六步推导 + 三个数值实验（`exp_quantization_selfref.py`）：
- ExpA：位置 $X=\sigma_3$（对角）、动量 $P=\sigma_1$（非对角），$[X,P]=2i\sigma_2\ne0$（自指两端必然不对易）；
- ExpB：断裂贡献 $|\mathbf d|^2$ 离散化，最小非零单位 $=\varepsilon^2$（断裂最小单位 = 一个量子 = $\hbar$）；
- ExpC：投影观测 → 坍缩（叠加态重叠 0.346 → 本征态 1.000）。

**空间涌现第一步（D 自发长三维）**——Hopfion 的前提（Hopf 荷需要三维）：
- 暖启动（`exp_space_3d.py --warm-start`）：hit=1.000、geo=0.001 → **三维近邻是 D 的稳定解**（已钉死）；
- 纯随机（`exp_space_3d.py`）：hit=0.62，卡局部极小；
- 退火（`exp_space_3d_anneal.py`）：hit 0.62→0.82–0.90，自发**接近**三维，但精确难。

**迹反常（`exp_trace_anomaly.py`）**——负结果：0 维 D 模型没有量纲嬗变（迹反常需要时空维度+能标+RG 流，单个矩阵是 0 维）。这排除了"在 0 维 D 里算迹反常"这条路，指向"尺度从反常来"需要先有空间（= 空间涌现）。

**核心结论 + 下一步**：
- 量子性 = 断裂 + 自指 = 自指观测，已有推导+数值锚点（本体论内核立住）；
- 空间涌现：复现三维已钉死（暖启动 hit=1.0）；自发选维缺"维数压力项"——极可能就是 Hopf 拓扑本身（只有三维能承载 Hopf 荷，所以 D 自动落向三维）。下一步：Exp6a（Hopf 荷=谱流）。

## 下一步

- **已闭环**：量子化自指（断裂+自指→量子性）+ 欧氏涌现（Exp1）+ 号差（Exp3）+ 拓扑零模（Exp4 绕数、Exp5 陈数）
- **已钉死**：三维近邻是 D 的稳定解（暖启动 hit=1.0）
- **下一站**：Exp6a（Hopf 荷 = 谱流，odd Chern-Simons 绕数，不是 η）→ 维数压力项（偏好 Hopf 荷，让 D 自发掉三维）
- **远景**：玻恩规则、真统一

## Layout

```
self_ref_spacetime/
  src/           # D 参数化、最短路、作用量、指标、闭合环磁通
  experiments/   # Exp1 环近邻；Exp2 Tr(D^4)/相位；Exp3 Dirac 号差；Exp4 拓扑零模；Exp5 陈数；Exp1 v2/v3/v4 磁通项；exp_trace_anomaly 迹反常；exp_quantization_selfref 量子化自指；exp_space_3d(+anneal) 空间涌现
  tests/         # 已知环 D → d 与解析一致
```
