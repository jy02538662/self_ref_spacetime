# Self-Referential Spacetime — Numerical Closure Tests

「自指时空」的可证伪数值实验合集：一个自反性关系算子 $D$（$D_{ij}=D_{ji}^*$）如何涌现出距离、号差、拓扑、空间维数与量子骨架。

> **归档状态（2026-09-03）**：项目处于**封存期**（考试复习，约至 2026-11），不再新增实验，仅保留可复现的记录。理论完整路线见 vault 笔记 [[量子潮水理论行动指南 v3（自指主线版）]] 与 [[断裂与自指：量子性的自指来源（严格推导+开放问题+实验任务）]]。

**v0.1 命题（环）**：固定环型目标距离 $L_{ij}$ 时，

$$
S=\mathrm{Tr}(D^2)+\lambda\sum_{i<j}(d_{ij}-L_{ij})^2
$$

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
python -m experiments.exp_spectral_flow
python -m experiments.exp6a_hopf_charge
python -m experiments.exp_space_dim_compare
python -m experiments.exp_eta_framing
python -m experiments.exp_spectral_dim_compare
python -m experiments.exp_cold_start
python -m experiments.exp_spectral_dim
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

## Status

| Item | Result |
|------|--------|
| Exp1 环 + 近邻涌现 | **pass**（多峰；约 3/8 种子进干净环盆地） |
| Exp2 $\mathrm{Tr}(D^4)$ 虚耦合 | **fail**（等模环上相位不降 TrD4；自由相位不优于全实） |
| Exp3 Dirac 结构 → 洛伦兹号差 | **pass**（欧氏 $D^2$ 全非负；洛伦兹 $D^2$ 正负 48/48/32，不定 $(-,+)$） |
| Exp4 拓扑荷 → 零模（Jackiw-Rebbi） | **pass**（绕数 0 无零模；绕数 1 → 2 零模在畴壁，E≈1e-18） |
| Exp5 二维陈数 → 零模（量子霍尔） | **pass**（陈数 8 → 16 零模 = 陈数 × 2 费米子加倍） |
| Exp1 v2 闭合环磁通项 | **pass**（flux 惩罚把磁通从随机压到 $0$/$\pi$ Z₂ 通量，$0$/$\pi$ 简并） |
| Exp1 v3 曲率项定义对比 | **pass**（sin→Z₂、wilson→零磁通、cos→$\pm \frac{\pi}{2}$；曲率定义决定磁通零点） |
| Exp1 v4 热采样（MCMC） | **负结果**（温度驱动磁通有序化，但经典热不打破 Z₂ 简并） |
| 迹反常（exp_trace_anomaly） | **负结果**（0 维 D 模型无量纲嬗变——迹反常需时空维度+能标，单个矩阵是 0 维） |
| 量子化自指三实验（exp_quantization_selfref） | **pass**（$[X,D]=2i\sigma_2$ 不对易、断裂最小单位=$\varepsilon^2$、投影坍缩） |
| 空间涌现·暖启动（exp_space_3d --warm-start） | **pass**（hit=1.000、geo=0.001，三维近邻是稳定解） |
| 空间涌现·纯随机（exp_space_3d） | 卡局部极小（hit 0.62） |
| 空间涌现·退火（exp_space_3d_anneal） | 自发接近三维（hit 0.82–0.90），精确难 |
| 谱流 = 绕数（exp_spectral_flow） | **pass**（畴壁零模精确存在 E≈1e-16，$\lambda\sigma_y$ 驱动下穿越 0，能级追踪捕获；净谱流=绕数的整数不变量验证留 Exp6b） |
| **Exp6a Hopf 荷 = odd Chern-Simons**（exp6a_hopf_charge） | **pass**（标准 Hopf 映射 → 0.99996；平凡场 → 0；Q=2 → 4=$Q^2$ 复合律；FFT 解矢势） |
| 空间涌现·维度对比（exp_space_dim_compare） | **负结果**（暖启动 2D/3D/4D 全 hit=1.0——三维不特殊；冷启动 2D>3D>4D——复现几何框架偏向低维） |
| **3/2 = 谱流 + $\eta$**（exp_eta_framing） | **pass**（1D 谱流=1、热核 $\eta \to \pm \frac{1}{2}$、合成 3/2 = 整数 Hopf + 分数自旋） |
| 谱维数区分维度（exp_spectral_dim_compare） | **pass**（热核 $K(t)=\mathrm{Tr}(e^{-tL})$ 读出 $d_s$：2D→2.01、3D→3.01、4D→4.02，不预设坐标） |
| 谱维数·早期版（exp_spectral_dim） | 谱维数读出的早期实现，已被 exp_spectral_dim_compare 取代（后者更干净、可解析） |
| 空间涌现·冷启动（exp_cold_start） | **负结果**（随机 D + 纯 $\mathrm{Tr}(D^2)$ 不自发选维、只塌缩——暴露「自发选维」需要一个维度无关的维数压力项 = 第三层） |
| **弦网凝聚 sanity check**（exp_string_net_condensation） | **已知工具复现（非推进，封存）**（Part A：SU(2)_k 的 TEE 目标 log D；Part B：张力→0 环密度 0→0.5、参与率 1→460 = 凝聚机制；Part C：加权 B_p → $d_j^2$/D 分布，真空被压到 1/D） |
| **Born 偏离判据演示**（exp_born_deviation） | **pass**（自指：p_Born SO(3) 协变，偏离 6e-16；引入偏好方向 e 后偏离 0.50；偏离随 e 强度线性缩放、e=0 归零） |
| **纠缠熵·面积律 vs 体积律**（exp_entanglement_area_law） | **pass**（1D 链：局域基态 S~0.15 ln L 次体积、随机态 S~0.38 L 体积律——演示纠缠熵标度 = 局域性签名） |
| **纠缠熵·面积律 + 维度**（exp_entanglement_area_law_dim） | **pass**（2D：S~l^1.25 周长律非面积；3D：S~l^2.25 面积律非体积——指数 = 维度） |
| **面积律 = 局域性判据**（exp_area_law_locality） | **pass**（1D 链有费米面→面积律 S~0.85→1.17；随机图无费米面→体积律 S~6.7→38.6；判据 = 费米面非跳跃范围） |
| **裸重连动力学**（exp_bare_reconnection） | **负结果**（4-正则图上满填充圈 + 局部重连，全部接受→均匀游走→弦统计乱跳、不组织；印证「裸动力学无偏好不产生结构」） |
| **相位当联络**（exp_phase_connection） | **pass（负倾向）**（相位=U(1) 联络，Wilson 环/面通量是拓扑对象；但麦克斯韦项 $-\sum_p\cos\Phi_p$ 极小化到平坦/平庸——三维靠非阿贝尔 Hopf/链接） |
| **Wilson 全局不变量**（exp_wilson_invariant） | **pass**（$Q=\prod_C W(C)=\prod_{\text{edges}} e^{i\theta}$ 严格重连不变量，Q_wilson=Q_direct、5 万步偏差 4e-15——拓扑荷的种子） |
| **SU(2) Wilson 环**（exp_su2_wilson） | **pass**（U(1) 升级 SU(2) 后单条 trace 会变、Fierz/skein 关系成立 2e-16——U(1) 种子 → SU(2) 相对链接，为什么 SU(2) 的第二张脸） |
| **群交换子 L**（exp_group_commutator） | **pass + 边界**（$L=1-\frac12\mathrm{Tr}([W_1,W_2])$：U(1) 恒 0、SU(2) 非 0（均值 0.748）——但 L 是「非对易」度量不是「链接数」，2D/3D 都非零，链接数需 Jones 多项式） |
| **辫子 R 矩阵**（exp_yang_baxter） | **pass**（SU(2) R 矩阵满足 QYBE + braid 关系（2.5e-16），本征值 $3\oplus 1$——交叉 → 辫群的第一块砖，Yang–Baxter = 「不分先后」） |
| **Jones 多项式**（exp_jones） | **pass**（Kauffman + writhe 区分平凡链接(2D) vs Hopf 链接(3D)：V(trivial)=-1.73 vs V(Hopf)=i——一个可算量区分 2D/3D） |
| **配对 ↔ 辫子交叉**（exp_pairing_to_braid） | **pass**（skein：交叉=配对{1,e}的量子叠加 $\sigma=A\cdot 1+A^{-1}e$，正负=系数 $A\leftrightarrow A^{-1}$ 互换；R 本征值 $3\oplus 1$=通道非正负；纯 FPL 构型 $Kauffman=d^{c-1}$ 只依赖闭合弦数——补上「配对→交叉」这一概念缺口的精确版） |
| **辫子词 → Jones**（exp_braid_word_to_jones） | **pass**（正向闭环：辫子词→skein→配对→Markov迹(数圈)→Kauffman→Jones；Hopf/Trefoil/Figure-8 对教科书值，R-II/R-III 不变性成立——B 路线第二~四步的正向链） |
| **A 从哪来（量子化开关）**（exp_quantization_condition） | **pass（边界澄清）**（经典自洽不固定 $A$——任意 $A\neq 0$ 都满足 $R\cdot R^{-1}=I$，$A$ 是自由参数；量子化=单位根 $q^{k+2}=1$ 才固定 $A$，此时 $-d=A^2+A^{-2}$=量子维度（$k=1\to 1.732<2$），经典极限 $k\to\infty$ 恢复 2——第三层唯一剩下的开关 = 量子相干，非经典自洽能推出） |
| **量子化真正来自 Jones-Wenzl 截断**（exp_jones_wenzl） | **pass（机制澄清）**（Hecke 关系对任意 $q$ 自洽、不逼出单位根；逼出单位根的是 Jones-Wenzl 幂等元存在性——Chebyshev 量子维度 $\Delta_n$ 在 $\delta=-2\cos\!\left(\frac{\pi}{k+2}\right)$ 时 $\Delta_{k+1}=0$，$f_{k+2}$ 消失 = SU(2)_k 截断。「绕自己一圈不多不少」（$f_n^2=f_n$）才是「整个一不能看出破绽」的精确数学形式） |
| **关系自指 → 拓扑自指（桥的性质）**（exp_relation_to_topology） | **pass（边界澄清）**（厄米 $D_{ij}=D_{ji}^*$ 的谱投影 $P_k$ 自动幂等 $P_k^2=P_k$——"关系→拓扑"的免费桥（谱定理），但只给 $\{0,1\}$ 平凡投影，不给量子化；量子化=单位根需额外截断 $\Delta_{k+1}=0$（Chebyshev 零点），厄米性推不出——第三层墙 = 第二个付费桥） |
| **有限 N → 单位根（量子化桥）**（exp_finite_quantization） | **pass（数学恒等式）**（N 条最大弦 ⟹ $\delta=2\cos\!\left(\frac{\pi}{N+1}\right)$ 使 $\Delta_{N-1}=1$、$\Delta_N=0$（Chebyshev 零点恒等式，N=3..8 全部精确）⟹ 单位根 ⟹ level $k=N-1$；A_k 结构（N 个自旋 0..k/2）融合规则验证。「有限→量子化」是**证明的数学恒等式**；唯一剩下的开放输入 =「D 用满 N 条弦」= 区分所有节点 = v3 公设「观察=一次有向区分」） |
| 纯谱 / 3+1 选择 | out of scope |

## Pass criterion (Exp1)

- $N=6$ 或 $8$；$L_{ij}=\min(|i-j|,N-|i-j|)$（边长 $a=1$）
- $D$ 从随机厄米起步，不预设图
- 强边集合贴近环近邻；非近邻边平均 $|z|$ 明显更小

## Exp2 criterion (failed on N=6)

- Protocol A：虚相位图案的 $\mathrm{Tr}(D^4)$ 应严格低于全实
- Protocol B：warm-start 后自由相位应比全实更低 $\mathrm{Tr}(D^4)/S$

## Exp3 criterion (passed on Nt=Nx=8)

- $D=\gamma^0\otimes D_t+\gamma^1\otimes D_x$，$\gamma^0$ 反对称（$(\gamma^0)^2=-1$）、$\gamma^1$ 对称（$(\gamma^1)^2=+1$），交叉项因反对易相消 $\Rightarrow D^2=-D_t^2+D_x^2$
- 判据：欧氏（$\gamma$ 全对称）$D^2$ 全非负；洛伦兹（$\gamma^0$ 反对称）$D^2$ 正负对称（不定，号差 $(-,+)$），且有零本征值（光锥）

## Exp4 criterion (passed on N=400, xi=10)

- 1D Dirac 算子 $D=$ 动能（$-i\sigma_x\partial_x$ 虚反对称）+ 质量 $m(x)\sigma_z$，畴壁 $m(x)=m_0\tanh((x-x_0)/\xi)$
- 判据：均匀质量（绕数 0）→ 0 零模、有能隙；畴壁（绕数 1）→ 2 零模局域在畴壁（E≈1e-18）
- 零模数 = 绕数 × 2（朴素格点 Dirac 的费米子加倍）；这是 Jackiw-Rebbi / 指标定理的最小版：**拓扑荷 → 零模（= 拓扑缺陷 → 粒子）**

## Exp5 criterion (passed): 2D 陈数 → 零模（量子霍尔）

- 2D 格点 Dirac + 均匀磁场（陈数 C=8），谱 = Landau 能级
- 判据：零模数 = 陈数 × 2（费米子加倍）；陈数 8 → 16 零模
- 与 Exp4 同一条「拓扑阶梯」的 2D 级（1D 绕数 / 2D 陈数 = 单个 Dirac 算子的指标，零模=荷）；3D Hopf 荷不是指标、是 odd Chern-Simons 绕数（见 Exp6a）

## 理论缺口（号差为什么读不出）

厄米 $D$ 的 $D^2$ 本征值恒 $\ge 0$（就是 $\lambda(D)^2$）。所以用 $\langle D^2\rangle$ 的正负划分类空/类时，**不能**从有限维 $D^2$ 的谱符号读出洛伦兹号差。

更深一层：相位（实/虚）**不改变** $D^2$ 的谱——$D^2$ 的本征值是 $|\lambda(D)|^2$，对相位盲。所以"虚耦合给号差"在厄米框架下**数学上不可能**；号差需要洛伦兹 $\gamma$ 结构（$\gamma_0^2=-1$）或非厄米 $D$，不是"相位"能给的。号差若存在，只能出现在连续极限的微分算子解释里。

### $2\times 2$ 钉死：号差 = 有向性，不是虚相位

一条耦合（$2\times 2$ 块）

$$D=\begin{pmatrix}0&z\\ w&0\end{pmatrix}$$

平方 $D^2=zw\cdot I$，本征值就是 $zw$。三种写法：

| 耦合 | 矩阵元 | 厄米 | $D$ 本征值 | $zw$ | $D^2$ | 号差 |
|---|---|---|---|---|---|---|
| 对称（空间无向） | $D_{12}=D_{21}=a$ | 厄米 | $\pm a$ | $a^2>0$ | 正定 | $+$ |
| 虚耦合（时间厄米） | $D_{12}=-D_{21}^*=ia$ | 厄米 | $\pm a$ | $a^2>0$ | 正定 | 无 |
| 反对称（时间有向） | $D_{12}=-D_{21}=a$ | **非厄米** | $\pm ia$ | $-a^2<0$ | **负定** | $-$ |

厄米 $=$ $w=\overline z$，于是 $zw=|z|^2\ge0$：无论 $z$ 实还是虚，$D^2$ 恒正定。唯一让 $D^2$ 变负的是 $w=-z$（实、反对称、不取共轭），此时 $D$ 非厄米（$D^\dagger=-D$）。

**所以：号差 = 时间方向的"有向性"（实反对称/非厄米），不是"虚相位"。** 空间对称给 $+$，时间反对称给 $-$，合起来 $D^2 \to +\partial_x^2-\partial_t^2=\square$，洛伦兹号差 $(-,+)$ 自然出现。

## Exp2 的解析定理（$\mathrm{Tr}(D^4)$ 对相位盲）

带相位环是标准 A-B 能谱：本征值 $\lambda_k=2\cos(2\pi k/n+\Phi/n)$，$\Phi=\sum\phi_i$ 为总磁通。于是

$$\mathrm{Tr}(D^4)=\sum_k\lambda_k^4=16\sum_k\cos^4\theta_k$$

由 $\cos^4\theta=\frac18(3+4\cos2\theta+\cos4\theta)$，且 $n\nmid4$ 时 $\sum_k\cos2\theta_k=\sum_k\cos4\theta_k=0$，得

$$\mathrm{Tr}(D^4)=16\cdot\tfrac18\cdot 3n=6n\quad(=36,\ n=6)$$

**与所有环边相位分布无关**。这解释了 Protocol A 中 5 种图案 $\mathrm{Tr}(D^4)$ 全是 36.0，也说明 $\mathrm{Tr}(D^4)$ **数学上不可能**通过偏好虚相位来驱动时间方向。

## 闭合环磁通项（Exp1 v2 / v3 / v4，2026-09-01）

把「相位 = 规范联络」落成第一个对相位敏感的作用量项。动机：Exp2 已证谱矩 $\mathrm{Tr}(D^{2m})$ 对相位盲，故改走闭合环磁通 $\Phi_{ijk}=\phi_{ij}+\phi_{jk}-\phi_{ik}$（格点规范理论的威尔逊环），用复数乘积算 $\sin\Phi$ 避开 arg 分支切割。

**Exp1 v2**（`exp1_v2_flux.py`）：作用量加 $\mu\sum_{i<j<k}\sin^2\Phi_{ijk}$。结果：flux 惩罚把磁通从随机（~9.4）压到 0，但相位**没有**变实矩阵——它压到了「$0$ 或 $\pi$ 磁通」（Z₂ 通量态）。铁证：$\overline{\cos\Phi} = 1 - 2\,\mathrm{pi\_frac}$ 精确成立（4 个种子全对），说明每个三角形磁通是二值的（$0$/$\pi$）。

**Exp1 v3**（`exp1_v3_flux_compare.py`）：对比三种曲率定义 $F$，结论是「曲率项的定义决定磁通有序化到哪个零点」：

| 曲率 F | 零点 | 收敛到的磁通 |
|---|---|---|
| $\sin\Phi$ | $0$ 和 $\pi$ | Z₂ 通量（$0$/$\pi$ 简并） |
| $|e^{i\Phi}-1|$ | 只有 $0$ | 零磁通（实矩阵） |
| $\cos\Phi$ | $\pm \frac{\pi}{2}$ | 全三角形 $\pm \frac{\pi}{2}$ 最大磁通（quarter_frac=1.0） |

**Exp1 v4**（`exp1_v4_mcmc.py`）：Metropolis-Hastings 采样 $e^{-\beta S}$（**经典热**，不是量子）。结果：温度驱动磁通有序化（flux 从高温 ~9.5 随机 → 低温 ~0.5 冻结到 $0$/$\pi$），高温 pi_frac≈1/3 精确等于均匀磁通下 $\cos\Phi<-0.5$ 的概率；但 Z₂ 简并**未被热涨落打破**（低温跨种子 pi_frac≈0.5，无偏好）。已知坑：$\beta=10$ 步长 0.3 太大、accept≈0（未充分热化）；判据曾把 0.5 当随机基准（实为 1/3），主结论不变。

**核心结论（负结果，两个钉子）**：
1. 经典曲率惩罚**不能生成**非平凡拓扑——它只把磁通抹平到 F 的零点；哪个磁通被允许由 F 的定义决定，不由物理偏好决定。
2. 经典热涨落**不能打破** Z₂ 简并——打破 $0$/$\pi$ 简并 → SU(2)，需要量子相干（$e^{iS}$）或额外对称破缺，不是热。

> 这排除了「设计经典作用量让 SU(2)/Hopf 拓扑自发涌现」这条路，把「为什么 SU(2)」的答案从「作用量」推向「拓扑荷守恒/量子化」（= B4 迹反常 + 弦网凝聚）。已同步记入 vault 的《涌现时空模型：后续可证明与扩展方向》（排除的错路）与《字典观测点》（记录表）。

## 量子化自指 + 空间涌现（exp_quantization_selfref / exp_space_3d / exp_space_3d_anneal，2026-09-01）

**量子化自指（断裂+自指→量子性）**——本体论内核，六步推导 + 三个数值实验（`exp_quantization_selfref.py`）：
- ExpA：位置 $X=\sigma_3$（对角）、动量 $P=\sigma_1$（非对角），$[X,P]=2i\sigma_2 \ne 0$（自指两端必然不对易）；
- ExpB：断裂贡献 $|\mathbf d|^2$ 离散化，最小非零单位 $=\varepsilon^2$（断裂最小单位 = 一个量子 = $\hbar$）；
- ExpC：投影观测 → 坍缩（叠加态重叠 0.346 → 本征态 1.000）。

**空间涌现第一步（D 自发长三维）**——Hopfion 的前提（Hopf 荷需要三维）：
- 暖启动（`exp_space_3d.py --warm-start`）：hit=1.000、geo=0.001 → **三维近邻是 D 的稳定解**（已钉死）；
- 纯随机（`exp_space_3d.py`）：hit=0.62，卡局部极小；
- 退火（`exp_space_3d_anneal.py`）：hit 0.62→0.82–0.90，自发**接近**三维，但精确难。

**迹反常（`exp_trace_anomaly.py`）**——负结果：0 维 D 模型没有量纲嬗变（迹反常需要时空维度+能标+RG 流，单个矩阵是 0 维）。这排除了"在 0 维 D 里算迹反常"这条路，指向"尺度从反常来"需要先有空间（= 空间涌现）。

**核心结论**：
- 量子性 = 断裂 + 自指 = 自指观测，已有推导+数值锚点（本体论内核立住）；
- 空间涌现：复现三维已钉死（暖启动 hit=1.0）；档 2 维度对比（负结果）证「复现几何不选三维、偏向低维」；谱维数读步已落地（$d_s$=2/3/4 可读出）。

## 实验状态总结

- **已闭环**：量子化自指（断裂+自指→量子性）+ 欧氏涌现（Exp1）+ 号差（Exp3）+ 拓扑零模（Exp4 绕数、Exp5 陈数）
- **已钉死**：三维近邻是 D 的稳定解（暖启动 hit=1.0）；Hopf 荷 = odd Chern-Simons（Exp6a，0.99996）；谱维数读出（exp_spectral_dim_compare，$d_s$=2/3/4）；总拓扑荷 3/2 的 1D 演示（exp_eta_framing）
- **档 2 维度对比（负结果）**：复现几何不选三维、偏向低维

## Exp6a：Hopf 荷 = odd Chern-Simons（2026-09-02）

`experiments/exp6a_hopf_charge.py`——把「拓扑荷 = 谱不变量」从 1D/2D（绕数/陈数 = 零模，指标定理）推进到 3D（Hopf 荷 = odd Chern-Simons 绕数）。

方法：给定 $n(\mathbf x):\mathbb R^3\to S^2$，算 emergent 磁场 $B_i=n\cdot(\partial_j n\times\partial_k n)$（谱微分），FFT 库仑规范解 $\nabla\times A=B$（$A_\text{hat}=-i(k\times B_\text{hat})/|k|^2$），$Q_H=\frac{1}{16\pi^2}\int A\cdot B\,d^3x$。

| 场 | 预期 | 数值 Q_H |
|---|---|---|
| 平凡场 $n=(0,0,1)$ | 0 | 0（精确） |
| 标准 Hopf 映射（Q=1） | 1 | 0.9962（N=64,L=5）→ **0.99996**（N=128,L=10） |
| Q=2 场 | 4 | 3.99–4.00 |

- Q=2 → 4 不是 bug，是 **Hopf 不变量的复合律** $H(g\circ f)=(\deg g)^2 H(f)$：$n=\mathrm{inv.stereo}((z_1/z_2)^Q)$ 在目标空间复合了度 Q 映射。这**独立验证了算法能分辨不同拓扑荷**（1 vs 4）。
- Hopf 荷对盒子敏感：skyrmion 密度在无穷远 $1/r^4$ 衰减，$L=5$ 截断尾部（0.996），$L=10$ 才收敛。
- 诚实边界：Exp6a 是**预设三维格点**，第三层（D 自发长出 $n(\mathbf x)$）仍开放；谱流版（族指标/Bott）是 Exp6b。

## 谱流前置（exp_spectral_flow，2026-09-02 修正）

旧版用「均匀质量 m 从 -m0→+m0」算谱流得 flow=0——物理上对（1D 周期 Dirac 能级成对只触碰 0 不穿越）。修正为**畴壁 + $\lambda\sigma_y$ + 能级追踪**：畴壁束缚精确零模（E≈1e-16），$\lambda\sigma_y$ 给零模能量偏移，零模随 $\lambda$ 穿越 0 被能级追踪（贪心最近邻）捕获。「净谱流=绕数」作为整数不变量需闭环参数族 + 手征破缺，留 Exp6b。

## 空间涌现·维度对比（exp_space_dim_compare，2026-09-02）

`experiments/exp_space_dim_compare.py`——把目标距离 L 换成 2D/3D/4D 超立方曼哈顿距离，比较 D 复现各维度的难易（同一框架 $S=\mathrm{Tr}(D^2)+\lambda\cdot\mathrm{geo}$，只优化模长、相位固定 0，因作用量只依赖 $|z|$）。

| 维度 | N | 暖启动 hit | 冷启动 excess（归一化） |
|---|---|---|---|
| 2D | 25 | 1.0（稳定） | 0.64 |
| 3D | 27 | 1.0（稳定） | 0.59 |
| 4D | 16 | 1.0（稳定） | 0.55（baseline 修正后） |

**负结果**：① 暖启动 2D/3D/4D 全 hit=1.0——「三维是稳定解」不特殊；② 冷启动 2D>3D>4D——维度越低越好找（配位数越少 $\mathrm{Tr}(D^2)$ 越便宜）。**结论：复现几何框架不选三维、反而偏向低维；「三维特殊」必须来自 Hopf 拓扑压力项（第三层），不在复现几何这一层。** 印证「维数压力项 = Hopf 拓扑」是必要的。

坑：4D n=2 超立方全角点、无内部态，真实配位是 4 不是 2d=8；baseline 应按 `mean(len(neighbors))` 算。

## 3/2 = 谱流 + $\eta$ 分数部分（exp_eta_framing，2026-09-02）

`experiments/exp_eta_framing.py`——1D 数值演示「总拓扑荷 = 整数 Hopf（谱流）+ 分数自旋（$\eta$ framing）= 半整数 3/2」。

模型：$S^1$ 上 1D Dirac 算子 $D_A=-i\,d/d\theta+A$（周期边界），本征值 $\lambda_n=n+A$（$n$ 整数，精确）。

- **谱流（整数）**：一族 $D(t)=-i\,d/d\theta+t$，$t:0 \to 1$ 的本征值净上穿 0 次数 = **1**（数值逐点追踪）；
- **$\eta$-不变量（分数）**：热核正则化 $\eta_\varepsilon=\sum_n \mathrm{sign}(n+A)\,e^{-\varepsilon|n+A|}$，$\varepsilon \to 0$ 收敛到闭式 $1-2A$——$A=1/4 \to +1/2$、$A=3/4 \to -1/2$（自旋 framing 的 $\pm \frac{1}{2}$，数值偏差降到 $4\times10^{-7}$）；
- **合成**：$1+\tfrac12=\tfrac32$。

诚实边界：这是 1D 类比（演示机制）；真正的 3D Hopf 荷（整数=1）已由 Exp6a 钉死（0.99996），这里补的是 $\eta$ 的分数部分（1/2）。两者在 framed/APS 框架里合成 3/2。真正的 $S^3$ 上自旋 1/2 Hopfion 的 $\eta=3/2$ 需在具体 Dirac 算子上把 $\eta$ 完整算出来（未做）。

## 谱维数区分维度（exp_spectral_dim_compare，2026-09-02）

`experiments/exp_spectral_dim_compare.py`——谱层第一步「读」：谱维数 $d_s$ 从热核 $K(t)=\mathrm{Tr}(e^{-tL})\sim t^{-d_s/2}$ 读出，区分 2D/3D/4D，不预设坐标。

用周期边界 d 维超立方晶格的拉普拉斯 L 的解析本征值 $\lambda=4\sum_j\sin^2(\pi m_j/n)$，热核迹 $K(t)=[\sum_m e^{-4t\sin^2(\pi m/n)}]^d$，在中间 t 区间（$1\ll t\ll n^2$）拟合 $\log K$ vs $\log t$，斜率 $=-d_s/2$。

| 维度 | n | t 区间 | $d_s$ | 目标 |
|---|---|---|---|---|
| 2D | 100 | [10,100] | 2.009 | 2 ✓ |
| 3D | 60 | [10,100] | 3.014 | 3 ✓ |
| 4D | 40 | [10,60] | 4.018 | 4 ✓ |

**两个坑（重要）**：① $K(t)$ 里的 $D^2$ 必须是**图拉普拉斯**，不是邻接矩阵——邻接矩阵的谱在能带中心（$\sum\cos=0$）DOS 是常数，$\mathrm{Tr}(e^{-tA^2})\sim t^{-1/2}$，所有维度都读出 $d_s\approx1$（错）；② $t$ 要取**中间区**（$1\ll t\ll n^2$），不是「小 $t$」——小 $t$ 时 $K\approx N$ 常数（斜率 $\approx 0$）。

**诚实定位**：这是「读」步（谱维数读出 D 的维数），不是「选」步。谱维数本身不选三维，只是给了一个**不预设坐标的维度读数**，作为后续「维数压力项」的原料；「选」步（偏好 $d_s=3$）预设了三维，是循环，落在第三层（为什么 SU(2)）。

## 弦网凝聚 sanity check（exp_string_net_condensation，2026-09-04）

`experiments/exp_string_net_condensation.py`——把「自指闭弦 → SU(2)_k 弦网相」这个 if-then 的两个「if」**分开**做 sanity check，明确不组装成证明。理论路线（vault 笔记）是：自指闭弦 → 弦网相 → 谱三元组 → GR，其中「闭弦 → 弦网相」需要两个 if：① 张力 $T \to 0$（弦凝聚）；② 融合结构 = SU(2)_k（=「为什么 SU(2)」公设）。

**Part A：SU(2)_k 拓扑纠缠熵目标（解析）**——量子维数 $d_j=[2j+1]_q$、总量子维数 $D=\sum_j d_j^2$、TEE $\gamma=\log D$：

| k | 量子维数 $d_j$ | $D$ | TEE $=\log D$ |
|---|---|---|---|
| 1 | 1, 1 | 2 | 0.6931 |
| 2 | 1, √2, 1 | 4 | 1.3863 |
| 3 | 1, $\varphi$, $\varphi$, 1 | 7.2361 | 1.9791 |
| 4 | 1, √3, 2, √3, 1 | 12 | 2.4849 |
| 5 | 1, 1.802, 2.247, 2.247, 1.802, 1 | 18.5918 | 2.9227 |

（k=4 时 $D=12$、$\log12=2.4849$ 恰好撞 1.0 的「ln12」，是**巧合**，勿追。）

**Part B：toy 环气体凝聚（张力 T → 0）**——3×3 环面、1024 个闭合环构型，$H=T\cdot(\text{长度})-\lambda\sum_p(\text{面翻转})$。扫 T 大→小：环密度 0.002→0.500（单调），参与率 1.02→460（局域→离域）。**这是 if-1「张力→0 触发弦凝聚」的机制签名**，但它是**等权**（单一环型，$D=2$，toric-code 型），不是 SU(2)_k 加权凝聚。

**Part C：加权弦网凝聚（B_p 算符，单面，一步插入）**——弦网动能项 $B_p=\frac1D\sum_j d_j B_p^j$ 的基态按量子维数加权（正确），对比等权 $B_p=\frac1N\sum_j B_p^j$。单面环型分布：

| k | 加权 $p_j=d_j^2/D$ | 等权 $p_j=1/N$ |
|---|---|---|
| 2 | 0.250, 0.500, 0.250 | 0.333, 0.333, 0.333 |
| 3 | 0.138, 0.362, 0.362, 0.138 | 0.250, 0.250, 0.250, 0.250 |
| 4 | 0.083, 0.250, 0.333, 0.250, 0.083 | 0.200, 0.200, 0.200, 0.200, 0.200 |

**加权凝聚把真空压到 $p_0=1/D$、抬高高量子维数自旋**——$d_j^2/D$ 分布正是弦网凝聚的指纹，也是 $\log D$（TEE）的来源。这是往 Levin-Wen 靠的第一步（实现了 B_p 的量子维数加权），但仍是单面/一步插入，**没有 F-符号、没有顶点融合规则**。

**定位与边界（诚实）**：Part A 算已知指纹、Part B 演示已知凝聚机制、Part C 验证已知加权分布——三者都在**已知理论（Levin-Wen 弦网）内部**复现公式，**没有让自指 D 自己去选 SU(2)**，因此**没有推进「自指闭弦 → SU(2)_k 弦网相」**。该路线**封存**（不再加 F-符号、不再扩 Levin-Wen），本实验仅作为「已知工具复现」的存档。自指闭环（自反性→断裂→SU(2)→三维/Hopf）这条推导链保持完整、依然成立，是理论最硬的部分。

## Born 偏离判据演示（exp_born_deviation，2026-09-04）

`experiments/exp_born_deviation.py`——把 [[Born规则偏离判据（时间戳短稿）]] 的「弱命题 + 偏离判据」落成可跑的量。

**弱命题（小定理）**：自指（无外部观察者）= 内部 SO(3)（断裂→SU(2) 伴随作用）无偏好方向 ⇒ 坍缩概率 SO(3) 协变 ⇒ $p=f(\hat m\cdot\hat n)$（只依赖态与测量的夹角）。

**演示结果**：

| 部分 | 内容 | 结果 |
|---|---|---|
| Part A | $p_{\text{Born}}=(1+\hat m\cdot\hat n)/2$ 在 50 个随机 SO(3) 转动下 | 偏离 6e-16（= 0，SO(3) 协变） |
| Part B | 引入偏好方向 $\hat e$：$p_{\text{corr}}=p_{\text{Born}}+0.3(\hat m\cdot\hat e)$ | 偏离 0.50（非零，且实测 = 解析 $\varepsilon|(R\hat m)\cdot\hat e-\hat m\cdot\hat e|$） |
| Part C | 偏离随 $\varepsilon$ 扫描 | 线性缩放，$\varepsilon=0$ 时归零 |

**判据**：转动整个系统（态 + 测量轴一起转），概率不变 = 自指（无外部观察者）；概率变 = 有绝对方向 $\hat e$ 混入（外部参考系/观察者），偏离度 = 外部性的度量。

**诚实边界**：演示的是弱命题（SO(3) 协变 ⇒ $f(\hat m\cdot\hat n)$）+ 偏离判据，**不是** Born 规则的推导——$p_{\text{Born}}=(1+x)/2$ 的 $\cos^2$ 形式是**作为输入**的。所以不撞 Gleason/Zurek（他们证明「为什么是 $\cos^2$」，这里证明「自指 ⇒ 概率不能偏好绝对方向」）。

## 纠缠熵·面积律 vs 体积律（exp_entanglement_area_law，2026-09-04）

`experiments/exp_entanglement_area_law.py`——纠缠线的第一个 toy：演示「纠缠熵的标度 = 局域性签名」（v3 六·六，详见 vault 笔记《纠缠线：从自指 D 到纠缠熵（面积律 vs 体积律）》）。

**公式**：约化密度矩阵 $\rho_A=\mathrm{Tr}_B|\psi\rangle\langle\psi|$，纠缠熵 $S(\rho_A)=-\mathrm{Tr}(\rho_A\ln\rho_A)$；自由费米子用 Peschel 公式（关联矩阵本征值 $\lambda$）：$S=-\sum_\lambda[\lambda\ln\lambda+(1-\lambda)\ln(1-\lambda)]$。

**结果**（1D 链 N=400）：

| L | S_local（Fermi 海，局域基态） | S_random（随机态） |
|---|---|---|
| 10 | 0.85 | 6.81 |
| 40 | 1.09 | 25.77 |
| 200 | 1.28 | 77.58 |

- 局域基态：$S\sim0.15\ln L$（≈ $\frac16\ln L$，边缘块 $c=1$）——**次体积，局域性签名**；
- 随机态：$S\sim0.38\,L$——**体积律，非局域**。

**演示**：$S(\rho_A)$ 的标度（对数 vs 线性）= 局域 vs 非局域的可测签名。

**诚实边界**：运动学纠缠（图结构划分 + ansatz 态），非背景独立；演示的是面积律/体积律的机制，不是「自指 D 自发产生面积律」（需局域 Hamiltonian + 基态 = 缺的动力学）。

## 纠缠熵·面积律 + 维度（exp_entanglement_area_law_dim，2026-09-04）

`experiments/exp_entanglement_area_law_dim.py`——纠缠线第二步：2D/3D 上 Fermi 海（基态 = Gibbs 的 $\beta\to\infty$ 极限）的纠缠熵，看面积律 + 维度。

**结果**（Fermi 海，角落子区 $l\times l$ [$\times l$]）：

| 维度 | S 标度 | $S/l^{d-1}$ | 结论 |
|---|---|---|---|
| 2D（40×40） | $S\sim l^{1.25}$ | 0.83→1.29（基本恒定） | 周长律（$S\sim l$），非面积 $l^2$ |
| 3D（12³） | $S\sim l^{2.25}$ | 0.76→1.00（基本恒定） | 面积律（$S\sim l^2$），非体积 $l^3$ |

（指数 1.25/2.25 = 1、2 加对数尾巴 = 费米面面积律的对数违反 $S\sim L^{d-1}\ln L$。）

**结论**：纠缠熵按**边界**（$L^{d-1}$）标度、非体积（$L^d$），且**指数 = 维度**——量子腿（纠缠）和几何腿（维度）第一次在同一计算里碰头。

**诚实边界**：图是预设的（2D 方格 / 3D 立方），态是基态（非有限温 Gibbs，有限温给体积律假信号）。证明的是「若 D 长在 d 维规则图上，基态纠缠熵指数 = d−1（面积律）」，**不是**「D 自发变成三维」（还需图结构从 D 涌现 = 第三层/缺口三）。

## 面积律 = 局域性判据（exp_area_law_locality，2026-09-04）

`experiments/exp_area_law_locality.py`——把「面积律」做成「空间像不像空间」的谱判据：对比有费米面（局域）vs 无费米面（非局域）的图。

**结果**（Fermi 海基态，子区 = 前 L 个点）：

| 图 | S(L=10) | S(L=100) | 增长率 | 律 |
|---|---|---|---|---|
| 1D 链（有费米面） | 0.85 | 1.17 | 1.38× | 面积律（对数） |
| 随机图 G(N,0.3)（无费米面） | 6.7 | 38.6 | 5.7× | 体积律（线性） |

**判据**：**面积律 $\Leftrightarrow$ 基态有费米面（平移不变/局域）；体积律 $\Leftrightarrow$ 无费米面（随机/非局域）**。

**诚实修正（重要）**：长程跳跃 $1/r^\alpha$ **仍是平移不变的**（Toeplitz，本征态还是平面波），所以不管 $\alpha$ 多小都是面积律——**「局域性」（纠缠意义上）≠「跳跃范围」，等于「有没有费米面」**。真正非局域（体积律）是打破平移不变的随机图。

**诚实边界**：预设图 + Fermi 海基态；证明的是「面积律 = 费米面签名」，不是「自指 D 自发变局域」（需动力学）。

## 裸重连动力学（exp_bare_reconnection，2026-09-04）

`experiments/exp_bare_reconnection.py`——4-正则图上「满填充圈（FPL）+ 局部重连」的裸玩具动力学，问「闭合弦在 D 上只允许重连，会不会自己组织」。

**结果（负结果，结局 B）**：N=30、10 万步。闭合弦数量 2~10 乱跳、最大弦长 10~30 乱跳、平均弦长 5.6~30 乱跳（$n\times\text{avg}=2N$ 恒成立，自洽），**没有「少数大弦逐渐支配」的趋势**。

**机制**：全部接受 ⇒ $3^N$ 配对构型上的对称随机游走 ⇒ 收敛到均匀分布 ⇒ 无偏好大弦 ⇒ 不组织。

**裸动力学三连**（同一核心判断的三个证据）：① 经典作用量不能生成拓扑；② 经典热不能打破 Z₂ 简并；③ **裸重连不能自发组织**。→ **「裸」的经典动力学（无能量偏好、无量子相干、无拓扑守恒）都不产生非平凡结构；要产生结构，需拓扑守恒（Hopf/链接 = 第三层）。**

**诚实边界**：底层图固定 ⇒ $d_s$ 常数（且随机正则图无格点式幂律热核，$d_s$ 不是好观测量）；N=30 偏小、涨落大，但「不组织」结论不依赖 N。

## 相位当联络（exp_phase_connection，2026-09-04）

`experiments/exp_phase_connection.py`——把 $D$ 的相位 $e^{i\theta_{ij}}=D_{ij}/|D_{ij}|$ 当离散 U(1) 联络，Wilson 环 / 面通量当拓扑对象（接 Exp6a 的 emergent 规范场 $A$）。

**结果**：3×3 网格 4 个面，极小化麦克斯韦项 $S_{\text{top}}=-\sum_p\cos\Phi_p$ → 面通量全部落到 0 mod $2\pi$（平坦）、$S_{\text{top}}$ 打到下界 $-N_{\mathrm{faces}}$（= 负的面数，此处 $-4$）。

**结论**：① 相位当联络是对的（「拓扑荷从相位里长，不从配对里硬找」= 第三层正确答案）；② 但面通量是**阿贝尔 U(1) 曲率 = 麦克斯韦项，极小化到平坦/平庸**，不逼出拓扑。→ **拓扑荷 = 相位的 U(1) 曲率（平坦）+ 非阿贝尔 Hopf/链接（三维），后者才是「选三维」的来源**（Exp6a 已碰 odd Chern-Simons/链接）。

## Wilson 全局不变量（exp_wilson_invariant，2026-09-04）

`experiments/exp_wilson_invariant.py`——4-正则图上，全体闭合弦 Wilson 荷的乘积 $Q=\prod_C W(C)=\prod_{\text{edges}}e^{i\theta}$ 是重连下的严格全局不变量（每条边在闭合弦分解里恰好出现一次，与配对无关）。

**结果**：$Q_{wilson}=Q_{direct}$ 精确相等（arg 0.172388），5 万步重连后偏差 $4\times10^{-15}$（浮点 0）——**严格不变**。

**链条**：$D_{ij}=D_{ji}^*\Rightarrow U_{ij}=e^{i\theta_{ij}}\Rightarrow W(C)\Rightarrow Q$。这是纯图 D 上第一个「D 自己的拓扑不变量」（局部规范不变、重连守恒、不预设空间），是拓扑荷的**种子**（整体相位守恒），Hopf/链接（相对相位）是下一步。

## SU(2) Wilson 环：相对链接（exp_su2_wilson，2026-09-04）

`experiments/exp_su2_wilson.py`——把边相位从 U(1) 升级成 SU(2)（$2\times2$ 矩阵），Wilson 环 $W(C)=\mathrm{tr}\prod U_{ij}$（非阿贝尔，顺序进入迹）。

**结果**：① 单条 SU(2) 环迹非平凡、重连后环数 8→2→6→4、trace 一路变（**非不变**，对比 U(1) 的 $Q$ 严格不变）；② Fierz/skein 恒等式 $\mathrm{tr}(A)\mathrm{tr}(B)=\mathrm{tr}(AB)+\mathrm{tr}(AB^\dagger)$ 成立（偏差 2e-16）。

**意义**：U(1)（阿贝尔）给整体种子 $Q$，SU(2)（非阿贝尔）给相对链接（skein → Jones/Hopf）。**「为什么 SU(2」= 自同构（断裂→二元→SU(2)）= 相对链接所需的非阿贝尔群——同一个 SU(2) 的两张脸**。

**诚实**：SU(2) 的不变量是 Jones 多项式（非单一数字），「trace 会变」是预期。

## 群交换子 L（exp_group_commutator，2026-09-04）

`experiments/exp_group_commutator.py`——两条闭合弦的群交换子 $L=1-\frac12\mathrm{Tr}(W_1 W_2 W_1^{-1} W_2^{-1})$。

**结果**：SU(2)（非阿贝尔）5000 对随机矩阵，L 均值 0.748、非零比例 100%；U(1)（阿贝尔）L 恒 0（max 4e-16）。

**确认**：「为什么 SU(2」再次钉死——U(1) 阿贝尔 → 交换子恒 0（整体种子）；SU(2) 非阿贝尔 → 交换子非零（相对结构）。L 是「相对非对易性」的正确度量。

**边界**：L 是「非对易」度量、**不是链接数**——一般 SU(2) 相位下 2D 和 3D 都非零；真正的「2D=0/3D≠0」（Hopf 链接）需要 Jones 多项式/量子迹（Markov 迹），比交换子高一层。

## 辫子 R 矩阵（exp_yang_baxter，2026-09-04）

`experiments/exp_yang_baxter.py`——把「链接数」换成「辫子」：交叉的本质是辫子，不是上下/链接。SU(2) R 矩阵给「交叉一次」的矩阵。

**结果**：QYBE $R_{12}R_{13}R_{23}=R_{23}R_{13}R_{12}$ 和 braid 关系 $B_{12}B_{23}B_{12}=B_{23}B_{12}B_{23}$（$B=PR$）对所有 k 成立（偏差 2.5e-16）；braid 矩阵本征值 $\{q,q,q,q^2\}$ = 3 对称 + 1 反对称（$2\otimes2=3\oplus1$）。

**意义**：Yang–Baxter 关系 $\sigma_i\sigma_{i+1}\sigma_i=\sigma_{i+1}\sigma_i\sigma_{i+1}$ = 「不分先后/对立统一」的数学形式；SU(2) R 矩阵满足它，所以「交叉 → 辫群」的第一块砖钉死。下一步 = 辫子闭包迹 = Jones 多项式。

**诚实**：完整 R 矩阵是 4×4；2×2 是它的对角化（本征值），约定不同但 $3\oplus 1$ 结构不变。

## Jones 多项式（exp_jones，2026-09-04）

`experiments/exp_jones.py`——Kauffman 括号 + writhe 修正 → Jones 多项式，验证「一个可算量区分 2D/3D」。

**结果**：对所有 k，⟨trivial⟩ ≠ ⟨Hopf⟩、V(trivial) ≠ V(Hopf)。k=1（$A=q^{1/4}$，$q=e^{i\pi/3}$）：Kauffman $-1.732$ vs $-1.0$；Jones $-1.732$ vs $+i$。

**意义**：Jones 多项式区分平凡链接（二维可画）和 Hopf 链接（只能三维）——第三层「只差伸手」。

**诚实**：① Jones = Kauffman + writhe 修正 $(-A)^{-3w}$（不是「辫子闭包的普通迹」）；② 2D/3D 是嵌入区分（非抽象图性质）；③ 接到抽象图 D 上需先定嵌入/辫子词（最后一道窄门）。

## 配对 ↔ 辫子交叉（exp_pairing_to_braid，2026-09-04）

`experiments/exp_pairing_to_braid.py`——补上「配对 → 交叉」这一概念缺口的精确版。核心：4 度节点（端口约定 1,2=下、3,4=上）3 种配对 → 3 种构型：$(13)(24)$=直穿（$1$）、$(12)(34)$=平局（杯+帽，$e$）、$(14)(23)$=交叉（crossing）；skein 关系 $\sigma=A\cdot1+A^{-1}e$ 把「交叉」展开成前两种配对的量子叠加。

**结果**：① skein 自洽（$\check R=A\cdot1+A^{-1}e$、$\check R^{-1}=A^{-1}\cdot1+Ae$ 互逆 iff $e^2=de$，$d=-A^2-A^{-2}$）；② $\check R$ 本征值 $3\oplus1$ = 通道（对称自旋-1 vs 反对称自旋-0），**不是** $\sigma/\sigma^{-1}$；③ 正负交叉 = 同一对配对 $\{1,e\}$ 上系数 $A\leftrightarrow A^{-1}$ 互换；④ 纯 FPL 构型 Kauffman 括号 = $d^{c-1}$（只依赖闭合弦数）。

**意义**：交叉（$(14)(23)$）是 FPL 的第三种配对，但它是**无方向交叉**（不指定 over/under）；$\sigma/\sigma^{-1}$ 配对相同，区别只在 over/under，而这在 skein 叠加系数 $A$ vs $A^{-1}$ 里。→ 经典纯配对无符号，符号需量子相干（= B 路线「经典不能生成拓扑」的精确化）。

## 辫子词 → Jones 多项式（exp_braid_word_to_jones，2026-09-04）

`experiments/exp_braid_word_to_jones.py`——正向链的完整闭环：辫子词 → skein 展开 → 配对（TL 弦图）→ Markov 迹（并查集数圈）→ Kauffman 括号 → writhe 修正 → Jones 多项式。

**结果**：Hopf $\sigma_1^2$（$-A^4-A^{-4}$）、Trefoil $\sigma_1^3$、Figure-8 全对教科书值；**R-II**（$\sigma_1\sigma_1^{-1}=\mathrm{id}$）与 **R-III**（$\sigma_1\sigma_2\sigma_1=\sigma_2\sigma_1\sigma_2$）不变性成立——链环不变量的充要条件，算法被独立钉死。

**坑（三个）**：① `coeff=1j` 写成虚数单位，所有结果乘 $i$；② transposition 环数模型在「连续 cap 同一弦对」时错（transposition$^2$=identity 把 $e_i^2=d e_i$ 的两个独立圈坍缩成一个），改弦图并查集才修对；③ trefoil 手性约定（$\sigma_1^3$ = 左旋 trefoil，镜像无害）。

**意义（正向 vs 逆向）**：正向（辫子词 → Jones）完全闭环，拓扑荷是严格链环不变量。**逆向（纯配对 → 唯一辫子词）不唯一**——交叉配对 $(14)(23)$ 无方向、丢失 over/under。所以第五步「拓扑荷生成」的唯一缺口 = 经典纯配对如何获得叠加系数 $A$（量子相干），不是「配对怎么变交叉」。

## A 从哪来（exp_quantization_condition，2026-09-05）

`experiments/exp_quantization_condition.py`——交叉振幅 $A$（$\sigma=A\cdot1+A^{-1}e$）是自由参数还是被自洽逼出？

**结果**：① 经典自洽不固定 A——任意 $A \neq 0$ 都满足 $R\cdot R^{-1}=I$（数值 1e-16），$A$ 是自由参数；② 量子化 = 单位根 $q^{k+2}=1$ 才固定 $A=q^{1/4}$，此时 $-d=A^2+A^{-2}=2\cos\!\left(\frac{\pi}{2(k+2)}\right)$ = 量子维度（k=1→1.732<2）；③ 经典极限 k→∞ 恢复 2。

**意义**：$-d=A^2+A^{-2}$ 才是量子维度（**不是 $A$ 本身**）；$A$ 是自由参数，量子化（单位根）是**物理输入**，不能从经典自洽长出。第三层唯一剩下的开关 = 量子相干。

## 量子化真正来自 Jones-Wenzl 截断（exp_jones_wenzl，2026-09-05）

`experiments/exp_jones_wenzl.py`——量子化（单位根）的真正来源：Hecke 关系还是 Jones-Wenzl 幂等元？

**结果**：① Hecke 关系 $\sigma^2=(q-1)\sigma+q$ 对任意 $q$ 自洽（不逼出单位根）；② Jones-Wenzl 幂等元 $f_n$ 存在 $\Leftrightarrow$ Chebyshev $\Delta_{n-1} \neq 0$；③ SU(2)_k 的 $\delta=-2\cos\!\left(\frac{\pi}{k+2}\right)$ 使 $\Delta_{k+1}=0$（$f_{k+2}$ 消失）= 截断 = 单位根。

**意义（修正「Hecke 逼出单位根」的直觉）**：量子化机制 = **Jones-Wenzl 幂等元存在性（截断）**，不是 Hecke 关系。「绕自己一圈不多不少」（$f_n^2=f_n$）才是「整个一不能看出破绽」的精确数学形式。

## 关系自指 → 拓扑自指（exp_relation_to_topology，2026-09-05）

`experiments/exp_relation_to_topology.py`——$D_{ij}=D_{ji}^*$（关系自指/厄米）能否推出 $f^2=f$（拓扑自指/幂等）？

**结果**：① 厄米 $D$ 的谱投影 $P_k=q_kq_k^\dagger$ **自动幂等**（$|P_k^2-P_k|\sim10^{-16}$，本征值 $\{0,1\}$）= 谱定理，白送；② 两步分解——厄米→投影（免费/可推导），投影→单位根截断（付费/额外）；③ 截断 $\Delta_{k+1}=0$ 是 Chebyshev 零点（额外条件）。

**意义**：关系自指 → 拓扑自指是**两座桥**——免费桥（厄米→谱投影幂等，平凡 $\{0,1\}$）+ 付费桥（投影→单位根截断）。**第三层的墙 = 付费桥**，名字是「本征值尺度为什么恰好 $\delta=2\cos\!\left(\frac{\pi}{k+2}\right)$」（= 耦合强度为什么有界、界为什么是 2）。

## 有限 N → 单位根（exp_finite_quantization，2026-09-05）

`experiments/exp_finite_quantization.py`——从"有限"（N 条弦）推出量子化 $\delta=2\cos\!\left(\frac{\pi}{k+2}\right)$ 的桥。

**结果**：① **Chebyshev 零点恒等式**——$\delta=2\cos\!\left(\frac{\pi}{N+1}\right)$ 时 $\Delta_{N-1}=1$（$f_N$ 存在）、$\Delta_N=0$（$f_{N+1}$ 消失），N=3..8 全部精确；② 这个 $\delta$ 精确 = SU(2)_k 的 $2\cos\!\left(\frac{\pi}{k+2}\right)$，level $k=N-1$；③ $A_k$ 结构——N 个自旋 $0,\frac12,\dots,\frac k2$ 的量子维度满足融合规则 $[\frac12]^2=[0]+[1]$。

**意义（突破）**：**「有限 $N$ 条最大弦 $\Rightarrow$ $\delta$=单位根」是证明的数学恒等式**（Chebyshev 零点），所以「有限 → 量子化」从猜想变成了定理。唯一剩下的开放输入 = 「D 用满 N 条弦」（= 区分所有 N 个节点）= v3 压缩公设「观察 = 一次有向区分」。量子化的开关 = 这个公设本身。

## Layout

```
self_ref_spacetime/
  src/           # D 参数化、最短路、作用量、指标、闭合环磁通
  experiments/   # Exp1 环近邻；Exp2 Tr(D^4)/相位；Exp3 Dirac 号差；Exp4 拓扑零模；Exp5 陈数；Exp1 v2/v3/v4 磁通项；exp_trace_anomaly 迹反常；exp_quantization_selfref 量子化自指；exp_space_3d(+anneal) 空间涌现；exp_spectral_flow 谱流；exp6a_hopf_charge Hopf 荷（odd Chern-Simons）；exp_space_dim_compare 维度对比（负结果）；exp_eta_framing 3/2=谱流+$\eta$；exp_spectral_dim_compare 谱维数区分维度；exp_cold_start 冷启动（负结果）；exp_spectral_dim 谱维数早期版（被 compare 取代）；exp_string_net_condensation 弦网凝聚 sanity check（Part A log D / Part B 张力凝聚 / Part C 加权 B_p）；exp_born_deviation Born 偏离判据演示（弱命题 SO(3) 协变 + 偏离判据）；exp_entanglement_area_law 纠缠熵面积律 vs 体积律；exp_entanglement_area_law_dim 纠缠熵面积律 + 维度（2D/3D）；exp_area_law_locality 面积律 = 局域性判据（费米面 vs 随机图）；exp_bare_reconnection 裸重连动力学（负结果：裸动力学不组织）；exp_phase_connection 相位当联络（U(1) 麦克斯韦→平坦，三维靠非阿贝尔 Hopf/链接）；exp_wilson_invariant Wilson 全局不变量（Q=∏W(C) 严格守恒 = 拓扑荷种子）；exp_su2_wilson SU(2) Wilson 环（U(1) 种子 → SU(2) 相对链接）；exp_group_commutator 群交换子 L（非对易度量，非链接数）；exp_yang_baxter 辫子 R 矩阵（QYBE/braid，交叉→辫群）；exp_jones Jones 多项式（区分平凡 2D vs Hopf 3D）；exp_pairing_to_braid 配对↔辫子交叉（skein：交叉=配对{1,e}的叠加，正负=系数互换）；exp_braid_word_to_jones 辫子词→skein→配对→Markov迹→Kauffman/Jones（Hopf/Trefoil/Figure-8 + R-II/R-III）；exp_quantization_condition A 从哪来（经典自洽不固定 A、量子化=单位根、量子维度、经典极限）；exp_jones_wenzl 量子化真正来自 Jones-Wenzl 截断（Hecke 不自洽逼出、幂等元存在性逼出单位根）；exp_relation_to_topology 关系自指→拓扑自指（厄米→谱投影幂等=免费桥、投影→单位根=付费桥）；exp_finite_quantization 有限 N→单位根（N 最大弦→$\delta=2\cos\!\left(\frac{\pi}{N+1}\right)$，Chebyshev 恒等式，level $k=N-1$）
  tests/         # 已知环 D → d 与解析一致
```
