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
```

## Status

| Item | Result |
|------|--------|
| Exp1 环 + 近邻涌现 | **pass**（多峰；约 3/8 种子进干净环盆地） |
| Exp2 \(\mathrm{Tr}(D^4)\) 虚耦合 | **fail**（等模环上相位不降 TrD4；自由相位不优于全实） |
| Exp3 Dirac 结构 → 洛伦兹号差 | **pass**（欧氏 D² 全非负；洛伦兹 D² 正负 48/48/32，不定 \(-+\)） |
| Exp4 拓扑荷 → 零模（Jackiw-Rebbi） | **pass**（绕数 0 无零模；绕数 1 → 2 零模在畴壁，E≈1e-18） |
| Exp5 二维陈数 → 零模（量子霍尔） | **pass**（陈数 8 → 16 零模 = 陈数 × 2 费米子加倍） |
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

## 下一步

- **已闭环**：欧氏涌现（Exp1）+ 号差（Exp3）+ 拓扑零模（Exp4 绕数、Exp5 陈数）——同一个 \(D\) 生成 距离 + 号差 + 拓扑
- **下一站**：三维 Hopf 荷 = 谱不变量（η 不变量 / 非交换几何指标，开放问题，需先厘清方向）
- **远景**：Hopfion 拓扑（量子潮水核心）

## Layout

```
self_ref_spacetime/
  src/           # D 参数化、最短路、作用量、指标
  experiments/   # Exp1 环近邻；Exp2 Tr(D^4)/相位；Exp3 Dirac 号差；Exp4 拓扑零模；Exp5 陈数
  tests/         # 已知环 D → d 与解析一致
```
