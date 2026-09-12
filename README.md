# Self-Referential Spacetime — Numerical Closure Tests

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](#)
[![numpy](https://img.shields.io/badge/numpy-%3E%3D1.24-green.svg)](requirements.txt)
[![scipy](https://img.shields.io/badge/scipy-%3E%3D1.10-green.svg)](requirements.txt)

「自指时空」可证伪数值实验合集：厄米耦合矩阵 D（D_ij = D_ji*）如何读出距离、号差、拓扑、空间维数与量子骨架。

> **状态**：研究原型 / **封存期**（约至 2026-11）。数值是可复现计算证据，**不是**已完成物理证明；不声称证明 3+1。  
> **联系**：王超 · 1186306891@qq.com  
> **理论**：vault《量子潮水理论行动指南 v3》《断裂与自指：量子性的自指来源》  
> **公式约定（GitHub 稳妥）**：正文与表格用 Unicode/ASCII；关键公式只用独立 `$$` 块。**禁止行内 `$...$`**（GitHub/部分预览会整段乱码）。

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
- 付费桥 2 精确化：3D π 磁通磁平移生成 Cl(3)，不可约分解证明 SU(2) 在**动量空间**（非实空间局域）

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
| 纯迹作用量选维 | `exp_pure_spectral` | **负结果** | 度约束解塌缩；模长项偏好 dimer 非环；Tr(D⁴) 恒等式 err=0 |
| 纯迹·锁4reg+曲率 | `exp_pure_spectral_4reg` | **正结果** | 冷启动自发 toroidal 3×3 + π 磁通 |
| **结构选择门·自指锁模长** | `exp_structure_gate_proof` | **正结果（全局最优已证）** | 模长相等（无外部观察者推论）→ Tr(D⁴) 模长项 Nd(2d-1)r⁴ 退化为常数 448 → 只剩 4 环 holonomy 选 π 磁通；把 8 条不可缩回环也翻 π（24 环全 π）→ Tr(D⁴)=**256**；**严格全局最优**：Tr(D⁴)≥ Nd²=256（Cauchy–Schwarz），4×4 torus 全 frustrate 取等（本征值全 ±2、D²=4I）→ 无需穷举。取等图唯一性（分类 W(16,4)）开放 |
| 纯迹·全局最优 | `exp_pure_spectral_anneal` | **正结果** | 暖启动 5/5 回 toroidal；basinhopping 命中；toroidal S=−71.1 全局最优 |
| 纯迹·谱维 vs N | `exp_pure_spectral_dim16` | **修正** | 谱维用拉普拉斯（非 D²）；随 N 向 2 逼近 |
| π 磁通零模 | `exp_pi_flux_dirac` | **正结果** | 4 个孤立 Dirac 点：零模 4/0 交替（N≡0 mod 4→4，N≡2 mod 4→0）；无磁通 Fermi 面 2N−2 |
| π 磁通质量 | `exp_pi_flux_mass` | **正结果** | 交错质量 m → 能隙 = m 精确 |
| π 磁通拓扑 | `exp_pi_flux_chern3/4`（+`edge`） | **负结果** | 陈数恒 0（dₓ=2cos² kₓ≥0）；边界态 0；t2 分数陈数 0.36–0.70（半金属） |
| 纠缠·面积律 | `exp_toroidal_area_law` | **正结果** | toroidal+π磁通走面积律 S~l^1.15；随机图体积律 |
| 链接 2D/3D | `exp_linking_2d_3d` | **正结果** | Gauss 链接数 2D=0、3D(Hopf)=±1（缠绕只能 3D） |
| 自指 π 磁通 | `exp_pure_spectral_selfref` | **正结果** | 纯迹墨西哥帽（无曲率项）自然产生 π 磁通 0.5——自指（割裂）产生曲率，非设计 |
| Tr(D⁶)/Tr(D⁸) | `exp_pure_spectral_selfref2` | **负结果** | 更长闭合弦也不选维（谱维 0.2–0.3） |
| SU(2) 版纯迹 | `exp_pure_spectral_su2` | **负结果** | 每条边=2×2 厄米矩阵（断裂最小单元），偏向完全图，不选维 |
| 谱 k2 简并（Exp7） | `exp7_k2_degeneracy` | **负（随机相位）/正（π磁通偶尺寸）** | 随机相位→无 k2（GUE 全单态）；π 磁通 N=16/36/64 自对偶（T²=−1，J 验证~0）、N=9/25/49 否；纯环面/随机相位从不自对偶 |
| 谱 k2 交错质量 | `exp7_k2_staggered_mass` | **正结果** | Semenoff 质量保 T²=−1（全程自对偶，Kramers 内禀）；N=36 无零模仍自对偶；E=±√(m²+E₀²) 精确验证；但 ±2 的 4 重简并抬不掉（需破 C₄） |
| 谱 k2 破 C₄ | `exp7_k2_break_c4` | **正结果（纯 k2 需两者）** | 各向异性 t_y=1+ε：±2(4) 劈成 ±2t_x(2)+±2t_y(2)，公式 E=±√((2t_y cos kx)²+(2t_x cos ky)²) 精确验证；但 0(4) Dirac 点不劈 → **纯 k=2 = 破 C₄ + 交错质量 两者** |
| 谱 k2 自发 vs 手放 | `exp7_k2_spontaneous` | **正结果（手放=诊断）** | 墨西哥帽在对称点曲率>0（各向异性+260、质量+92）→ 破 C₄/加质量都是**手放**；但 π 磁通奖励∝t_x²t_y² 锁死各向同性 → **T²=−1 自发、清理靠手放（诊断性）** |
| 谱 k2 完整 SU(2) 代数 | `exp7_k2_su2_algebra` | **正结果（自动定理）** | Kramers ⟹ spin-1/2 SU(2) 代数（J_i=V⨁σ_i/2V†，四验证~1e-15，J²=3/4，J_z=±1/2）；形式/涌现（pseudospin、稠密）非物理/输入；命门收窄为唯一性 |
| 唯一性（方向 1） | `exp7_k2_uniqueness` | **正结果（靠推定理）** | 号差（γ⁰）与 SU(2)（J₀）是同一对象 [[0,1],[-1,0]] 的两面，是「有向区分」（J²=−I）唯一最小实现（a=±1）；收窄到「有向⟹反对称」 |
| S²→Hopf→三维 | `exp7_k2_hopf` | **正结果（标准拓扑接上）** | SU(2)≅S³（四元数）→ 断裂→S²（Hopf map）→ Hopf 荷（π₃(S²)=ℤ）→ 只在 3D；三维内建于 SU(2)≅S³ |
| Clifford 代数/自旋联络 | `exp7_k2_clifford` | **正结果（定位缺口1）** | 有向区分+SU(2)≅S³ → 4×4 Dirac，Clifford 号差 −+++ 误差 0；2×2 装不下 Cl(3,1) → 绳子=「2×2→4×4」（自旋几何桥） |
| 自旋联络 2×2→4×4 自发 | `exp7_k2_chiral` | **正结果（缺口1 已打通）** | π 磁通 D 自发携带手征 Γ（bipartite，谱 ±E 误差 ~1e-15）+ Kramers T²=−1 → AZ 组合 = DIII 类 = 4×4 Dirac；零模四重 = 2 手征 × 2 Kramers；Dirac 型 J（ΓJ=−JΓ）构造验证 ~1e-15。下一道门 = 焊成规范协变 ω_μ |
| 代数 A·一阶条件 | `exp_algebra_A_firstorder` / `_bimodule` / `_edgespace` | **正结果（三定理 + 全扫描）** | ① 模长相等 = 无偏好方向（S_N 传递 ⟹ 常数）；② 一阶条件 ⟹ 连通图上对角交换 A 平凡（[[D,a],b]ᵢⱼ=Dᵢⱼ(aⱼ-aᵢ)(bⱼ-bᵢ)）；③ M_N ⟹ D=0（万物统一→单代数过强）。数值全谱：A=ℂ^16 违例 34、A=M_16 违例 322、A=ℂ 平凡 0；双模 J=K 退化成单边、J=T(Kramers) 与点代数不兼容；路径代数嵌顶点→D=0；边空间 ∂+∂* 四候选全灭 ⟹ **标准一阶条件对离散关系无解，唯一出路 = 标准模型有限部分（A=ℂ⊕ℍ⊕ M₃ + 电荷共轭 J）** |
| 自旋联络·方向1 | `exp_spin_connection_lattice` / `_localJ` / `exp_pi_flux_3D` | **正结果（规范律验证 + 区分 + 框架攻坚 + S²序参量 + 3D）** | 框架联络 Uᵢⱼ=gᵢgⱼ⁻¹：格点规范律 err 6.7e-16、连续极限 ω→ gω g⁻¹+gdg⁻¹ err 4e-10。**关键区分**：D 边相位是 U(1) π 磁通（16 plaquette holonomy=−1，Chern=8）、**不是** SU(2) 自旋联络。**框架攻坚**：最小 Kramers 对象 = plaquette（面）；局部 J_p（从局部 D_p 算）存在但 gauge，物理=均匀 π 磁通（Z₂）。**S² 序参量**：winding=0（Z₂=SU(2)中心→单点），非平凡 S² 是 3D 对象（Hopf 荷）。**3D π 磁通**（4×4×4）：Kramers + 三对反对易磁平移（3D 上下）。**总结论：代数 A/方向 1/框架/S² 四条线撞同一堵墙 = 全局 SU(2) 局部化 = 第三层（2D→3D，Hopf）** |
| 付费桥2·不可约分解 | `exp_pi_flux_3D_irreps` | **正结果（严格证明 SU(2) 在动量空间）** | 磁平移 Tₓ,T_y,T_z 生成 Cl(3)：三对反对易、Tᵢ² 对易（本征值 ±1）、χ=TₓT_yT_z 中心。不可约分解 = 32 个 2 维泡利块，全由动量标签 (sₓ,s_y,s_z,χ) 组织；位置判据（PR=32、var(x)=1.25=均匀）证明**全是平面波（动量空间），非实空间局域**。→ 付费桥 2 = 动量空间 SU(2) → 实空间 S² 场 = Wannier 局域化（拓扑障碍 = Hopf 荷来源，未证） |
| 付费桥2·生成元局域化 | `exp_pi_flux_3D_localize` | **正结果（SU(2) 生成元 = 磁平移 = 非局域自旋）** | 全局 SU(2) 生成元（32 块泡利）变换到位置基 = 磁平移 Tₓ,T_y,T_z **本身**（恒等式；PR=1、对角线 fraction=0 = 无固定点置换）。→ π 磁通 SU(2) = **动量 SU(2)**（平移生成元），非局域自旋；付费桥 2 = 动量 SU(2) → 自旋 SU(2) = 自旋联络 ω_μ |
| 付费桥2·第一堵墙 | `exp_pi_flux_3D_first_wall` | **正结果（SU(2) 生成元 = link 算子，非 site）** | PₓTᵢPₓ=0 精确（单点投影 = 0），Tᵢ 对角权重 = 0、非零 = N 条边 → 磁平移是 **link 算子**（边），非 site 算子（点）。付费桥 2 = link SU(2) → site SU(2) = Wannier 局域化（拓扑障碍 = Hopf 荷） |
| 付费桥2·曲率/ω_μ | `exp_pi_flux_3D_curvature` | **正结果（ω_μ 已内建在磁平移里）** | 曲率 F_μν=[T_μ,T_ν]=2T_μ T_ν（反对易，无迹、SU(2)-like），F²=-4s_μ s_ν I。→ T_μ = 平移 + U(1) 磁通 + SU(2) 自旋的**统一**；ω_μ **已内建**，付费桥 2 = **分离**（SU(2) 从 U(1) 拆出并局域化），非「从无到有构造」 |
| 付费桥2·分离不可分 | `exp_pi_flux_3D_separation` | **正结果（SU(2) 是 U(1) 的涌现，不可分）** | 带磁通 {T_x,T_y}=0（反对易→SU(2)），无磁通 {T_x,T_y}=16（对易→无 SU(2)）。→ SU(2) 是 U(1) 相位的涌现，去掉相位就没有 SU(2)；**纯 SU(2) 的 ω_μ 在 π 磁通里不存在**，非平凡 ω_μ 需换非平凡拓扑体系 |
| 逼三维·Jones 判别器 | `exp_force3d_jones_toy` | **正结果（绕远但可用）** | Jones 区分 2D 平凡（σ⁰→d）vs 3D Hopf（σ²→−A⁴−A⁻⁴）；TL₂ Markov 迹递推，err ~1e-15 |
| 逼三维·作用量项 | `exp_force3d_jones_action` | **正结果（绕远但可用）** | S_link=−ν\|V(σ^k)−V(σ⁰)\|：σ⁰（2D）奖励 0、σ²（3D）奖励 −1.85 |
| 逼三维·读 k | `exp_force3d_read_k` | **正结果（绕远但可用）** | 环结构→交叉→Gauss 链接（2D=0/3D=−1）→k→S_link |
| **逼三维·收口（上下）** | `exp_force3d_updown` | **✅ 最终结论** | **上下 = π 磁通反对易 T_xT_y=−T_yT_x**（|anti|≈0、|comm|=2、holonomy=−1），在 D 自己结构里、不预设；反对易=非阿贝尔=SU(2)=内部三维（SU(2)≅S³）→ **逼三维闭环** |
| 尺度·边缘间隔 | `exp_scale_chebyshev` | **正结果** | Chebyshev 零点边缘间隔 = 3π²/N² ~ N⁻²（尺度生成）；中间 2π/N |
| 圈气体配分 | `exp_partition_unitroot` | **负结果** | 圈气体 Z(A)=Σd(A)^c 在单位根处无峰值，偏向经典极限（阶段 2a 收口证据） |
| 相位干涉 | `exp_phase_interference` | **负结果** | plaquette 相位积分 e^{iS} 可分离 → Φ 均匀（阶段 2a 收口证据） |
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
- **负结果钉子**：虚相位不给号差；经典曲率/热难产非平凡拓扑；裸重连无组织；纯 Tr(D²) 不选维；复现几何偏向低维（文档）；纯迹 Tr(D⁴) 压不出环（度约束下模长项偏好 dimer）；随机相位不自发 k=2（Exp7，GUE 全单态）
- **已钉死（归档）**：π 磁通偶数尺寸（N=16/36/64）自对偶（T²=−1 Kramers 结构，J 验证误差~0）——「k=2 二元」藏在 π 磁通的涌现时间反演里，非随机相位（Exp7）
- **已钉死（2026-09-11）**：3D π 磁通 SU(2) 在**动量空间**（磁平移 Tₓ,T_y,T_z 生成 Cl(3)，不可约分解 = 32 个泡利块全由动量标签组织、位置弥散 var(x)=均匀值）；且全局 SU(2) 生成元 = 磁平移 Tᵢ **本身**（恒等式，无固定点置换 = 非局域自旋）——付费桥 2 = 动量 SU(2) → 自旋 SU(2) = 自旋联络 ω_μ
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

### 纯迹作用量选维（`exp_pure_spectral_last_run.json`）

作用量 S=-α Tr(D²)+β Tr(D⁴)+γΣᵢ[(D²)ᵢᵢ-c]²（冷启动、无预设距离）。

- 恒等式验证（环 C6 / K6 / 随机厄米）err 全 = 0；
- 原始参数（N=6, α=2, β=0.5, γ=10, c=1）：4 seed 全停 S≈-9.14，d_h≈1.048，结构乱跳，**非近邻环**；
- 解析比较（dimer / ring / complete）：N=6,c=1 → −9.0 / −7.5 / +0.6；N=6,c=2 → −12.0 / −6.0 / +26.4；**dimer 恒最优**。

结论：纯迹作用量解决了纯 Tr(D²) 塌缩（d_h 锁到 c），但选不出环——「选维」正杠杆是曲率（4 环 holonomy），不是 Tr(D⁴) 模长惩罚。详见 vault《纯迹作用量的精确分解与冷启动选维》。

### 纯迹·锁 4-regular + 曲率（`exp_pure_spectral_4reg_last_run.json`）

作用量 S=-α Tr(D²)+γΣ(d_h-c)²+δΣ(Σⱼ r_hj⁴-c²/m)²+νΣ rₚrod(1+cosΦ)（4 次度约束锁 4-regular + 曲率项奖励 π 磁通）。

- 杀 dimer 对照（`exp_pure_spectral2d`）：μΣ r⁴ 把结构推向完全图 K9（负结果）；
- 4 次度约束干净锁定 4-regular（排除 dimer/环/完全图）；
- 曲率项：8 seed 冷启动，seed=1 自发长出 **toroidal 3×3 + 9 个 π 磁通 plaquette**（度数全 4、89% π 磁通）；4 seed 4-regular、3 seed 局部极小。

结论：**「从纯代数 D 自发选二维 + π 磁通」可实现**，但收敛难（局部极小）+ 曲率项不唯一锁定 toroidal + 谱维数小图伪影（0.5–1.6，需 N=16 重测）。

### π 磁通支线（`exp_pi_flux_*.py`）

π 磁通 toroidal 的拓扑性质探查（严谨数值见 vault《纯迹作用量：精确恒等式与自发二维涌现（严谨数值记录）》第九节）。

- **零模数**（`exp_pi_flux_dirac`）：π 磁通 N=4→4、6→0、8→4、10→0、12→4（4/0 交替）；无磁通 2N−2 随 N 线性 → **4 个孤立 Dirac 点**（vs 无磁通 Fermi 面）。
- **能隙**（`exp_pi_flux_mass`，N=8）：交错质量 m → **能隙 = m 精确**（0.2→0.2、0.5→0.5、1.0→1.0、2.0→2.0）。
- **陈数**（`exp_pi_flux_chern3/4`，QWZ 校准 m=−1→−1、0→+1、1→+1、3→0）：交错质量及 σ_z/σ_y 型项陈数 = 0（dₓ=2cos²kₓin[0,2] 恒 ≥0，winding 恒 0）；边界态 0（`exp_pi_flux_edge`）；加实数次近邻 t₂ → 陈数分数 0.36–0.70（半金属）。

结论：**π 磁通正方格点陈数恒 0（拓扑平凡），要陈数拓扑须换 honeycomb（Haldane）或加自旋（Kane-Mele）**——「拓扑保护」这一环排除了正方格点。

### 自指产生曲率 + 逼三维（`exp_pure_spectral_selfref*`、`exp_toroidal_area_law`、`exp_linking_2d_3d`）

> 2026-09-08 后半段。严谨数值见 vault《纯迹作用量：精确恒等式与自发二维涌现（严谨数值记录）》第十～十五节。

- **面积律**（`exp_toroidal_area_law`）：toroidal + π 磁通费米海纠缠熵走**面积律**（Ssim l¹.15）、随机图体积律 → 长出的 toroidal 是「真空间」。
- **链接 2D/3D**（`exp_linking_2d_3d`）：Gauss 链接数 2D=0、3D(Hopf)=±1 → **缠绕（螺旋=双链）只能在 3D**。
- **自指产生 π 磁通**（`exp_pure_spectral_selfref`）：纯迹墨西哥帽 S=-αTr(D²)+βTr(D⁴)（**无度约束、无曲率项**），冷启动 π 磁通 0.48–0.67 → **自指（割裂）自然产生曲率/π 磁通，不需要设计曲率项**。
- **Tr(D⁶)/Tr(D⁸)**（`exp_pure_spectral_selfref2`）：加更长闭合弦，谱维仍 0.2–0.3 → 不选维（Tr(Dⁿ) 只数单条闭合弦）。
- **SU(2) 版纯迹**（`exp_pure_spectral_su2`）：每条边 = 2×2 厄米矩阵（断裂最小单元 = SU(2)），冷启动偏向完全图 → 不选维。

结论：**纯迹到不了三维（要 Hopf 拓扑 / 第三层）。** 三段机制完整——自指（割裂）产生曲率/π 磁通（U(1)）、断裂给 SU(2)（最小单元 2×2）、链接只在 3D——但「三维的自发生成」= 开放问题：纯迹不选维，要 Hopf 拓扑；「曲率/π 磁通的自发动力学」= 自指（割裂）产生曲率，动力 = 观察者缺陷视角（观察者代数 Type III→II）。见 vault [[纯迹作用量：精确恒等式与自发二维涌现（严谨数值记录）]]、[[字典观测点]] 十六号种子。

### 谱 k2 简并（Exp7，`exp7_k2_degeneracy`）

> v4 第十节「逼三维 / 第三层」的 ① Exp7。判据 = 谱的客观简并结构（无预设聚类度量），严格化为「自对偶/四元数」（∃反幺 T、T²=−1、[T,D]=0 ⟺ 全偶重数，Kramers 定理）+ 构造性 J 验证。严谨数值见 vault《D谱k2简并：本征态是否自发聚成二元SU(2)（严谨数值记录）》。

- **随机相位**（toroidal 模长 + 每条边随机 θ_ij，300 实例，N=9..64）：精确 k2 占比恒 0（全单态，⟨r⟩≈0.60 = GUE）；近似 2% 简并与 GUE 零假设持平 → **无 k=2 信号**。
- **π 磁通**（纯迹作用量长出的结构）：**偶数尺寸自对偶**——N=16/36/64 全偶重数、J 验证（反对称/幺正/对易/T²=−1）误差 ~0；N=9/25/49 含奇重数 → 非自对偶。规律：**自对偶 ⟺ n 偶**（Dirac 点在 (±π/2,±π/2)，落动量格 ⟺ n≡0 mod 4 → 零模；但自对偶只要 n 偶）。
- **交错质量分离**（`exp7_k2_staggered_mass`）：Semenoff 质量 m(-1)^i+j 保 T²=−1 → 扫 m 全程自对偶（**T²=−1 内禀**）；N=36 无零模仍自对偶（**T²=−1 不依赖 Dirac 点**）。能带 E₀=±2√(cos² kₓ+cos² k_y)，质量后 E=±√(m²+E₀²)（精确验证：0→±m、±2→±√(m²+4)、±2√2→±√(m²+8)）。但 **±2 的 4 重格点简并纹丝不动**（Semenoff 质量抬不动，要破 C₄/各向异性跳变）。
- **破 C₄ / 各向异性跳变**（`exp7_k2_break_c4`）：t_y→1+ε，±2(4) 劈成 ±2t_x(2)+±2t_y(2)，能带 E=±√((2t_ycos kₓ)²+(2tₓcos k_y)²) 精确验证（ε=0.3 实测 2.0/2.6/3.28 全中）。**但 0(4) Dirac 点不劈**（被 Dirac 结构保护、不靠 C₄）→ **「纯 k=2」= 破 C₄ + 交错质量 两者**（破 C₄ 抬 ±2 格点、质量 gap 0），得 8 个 Kramers 对。
- **自发 vs 手放**（`exp7_k2_spontaneous`）：墨西哥帽 S=-αTr(D²)+βTr(D⁴) 在对称点曲率 >0（各向异性 +260、质量 +92）→ 破 C₄/加质量都是**手放**。但 π 磁通奖励 ∝ tₓ²t_y²（AM-GM 各向同性最大）→ **同一个机制既自发产生 T²=−1，又锁死 C₄**。结论：**T²=−1 自发、清理（破 C₄+质量）靠手放（诊断性，非生成性）**——「纯 k=2」是露出已自发的 Kramers 结构，不是造出 SU(2)。
- **完整 SU(2) 代数**（`exp7_k2_su2_algebra`）：Kramers（全偶重数）⟹ spin-1/2 SU(2) 代数，**自动定理**——Jᵢ=V(⊕σᵢ/2)V^†，四验证（厄米/[Jᵢ,H]=0/[Jₓ,J_y]=iJ_z/Casimir J²=tfrac34）~1e-15，J_z=±tfrac12。但它是「形式/涌现」（pseudospin、J 稠密非局域），非「物理/输入」（spinless 无物理自旋）。**命门收窄为唯一性（为什么是 SU(2））。**
- **π 磁通 + 相位扰动**：eps=0.05 即碎掉 k2 → 二重简并是离散对称保护，非泛型。
- **纯环面（零相位）**：从不自对偶（格点平移简并混奇重单态，非 SU(2)）。
- **对照**：四元数正控制 k2=1.000、⟨r⟩=0.674（GSE）；GUE k1=1.000。

结论：**「k=2 自发涌现」不来自随机相位（GUE），而藏在 π 磁通的涌现时间反演 T²=−1 里**——且恰好是纯迹作用量已选出的结构；交错质量证明 T²=−1 **内禀**（不依赖 Dirac 点）。**纯 k=2 的配方已找到 = 破 C₄ + 交错质量**。诚实边界：自对偶只证「存在 T²=−1」（Kramers 特征），≠「SU(2) 已涌现」；且破 C₄/交错质量都是「可能手放」的项 → 真问题 = 二者是**自发还是手放**（对应自发晶格畸变 / 自发质量生成，都可测）。

### 量子化与尺度（`exp_scale_chebyshev`、`exp_partition_unitroot`、`exp_phase_interference`）

> 主线二（量子化 + 尺度）。严谨数值见 vault《量子化与尺度：有限N到单位根与边缘间隔（严谨数值记录）》。

- **尺度生成**（`exp_scale_chebyshev`）：Chebyshev 零点 2cos(π m/(N+1)) 的相邻间隔——边缘间隔 = 3π²/N²（系数精确 3π²，N→∞ → 0 即墙三）、中间间隔 = 2π/N。**尺度 = 有限 N 的边缘间隔**。
- **量子化收口**（`exp_partition_unitroot` + `exp_phase_interference`，两个负结果）：① 圈气体配分函数 Z(A)=Σ d(A)^c 在单位根处**无峰值**（偏向经典极限 d=2）；② plaquette 相位积分 e^iS 可分离 → Φ 均匀。结论：「量子涨落选单位根」是**伪问题**——量子化 = 有限性 + 截断（不是选择），公设完成闭环，不需要路径积分。

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

**时间 = 号差 + 旋转 + Wick + H（2026-09-09）**：有向区分 J=[[0,1],[-1,0]]=γ⁰=J₀（号差 + SU(2) 旋转，唯一性），生成 SO(2) 旋转 exp(tJ)；exp(2πJ)=I。Wick 转动 = γ⁰→iγ⁰（反对称↔虚耦合，号差 −↔+）；号差 −+++ 编码在 Dirac 结构 D²=-D_t²+D_x²，纯迹作用量 S[D] 自动给号差。**H ≠ γ⁰（纠正）**：H=γ⁰(γ·p+m)，三块落地——γ^μ 已推 / p=-i∇（Dirac 点线性色散 v=2 + 连续极限 p=-(i/2)(T-T⁻¹) 本征值 sin k→k）/ m=交错质量能隙=m，= 标准有质量 Dirac（石墨烯同款）。**π-flux 自发**：Tr(D⁴)=8r⁴(3+cosΦ) 在 Φ=π 极小。⚠️ 剩四道：① S[D] 自发选 π-flux toroidal 结构（M2）——**部分已证**（固定 4-regular+等模时 Tr(D⁴) 全局最小 = π-flux toroidal，Tr(D⁴)≥ Nd²=256 严格、无需穷举；仍开的是「S[D] 是否自发选到 4-regular+等模」这步 = 维度手选）；② Connes 一阶条件对 A=ℂ^N 不满足（**已定位**：对角交换 A ⟹ 平凡、M_N ⟹ D=0、标准一阶条件无解 → 需替换一阶条件本身，见 vault《代数A：模长相等的形式化与一阶条件的真坎》）；③ 「奖励 T²=−1」失败（dimer 也自对偶 [8,8]）；④ m→物理质量=尺度读出（数字巧合）。**诚实结论：磁通自发、维度手选。**

**Tr(D⁴)**：精确分解 Tr(D⁴)=2Σ d_h² − Σ r⁴ + Σ_{四不同指标} D_ij D_jk D_kl D_li。等模环上 =6n，相位盲（Exp2，因无 4 环）；有 4 环的图相位敏感，π 磁通让正方形 Tr(D⁴)=8r⁴(3+cosΦ) 从 32r⁴ 降到 16r⁴（奖励二维 plaquette）。**全局下界（2026-09-10）**：等模 d-正则图 Tr(D²)=Nd、Tr(D⁴)=Σ λ⁴ ≥ (Nd)²/N = Nd²（Cauchy–Schwarz）；4×4 torus 全 frustrate（24 环全 π）本征值全 ±2、D²=4I，取等 → 全局最优。

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
python -m experiments.exp7_k2_chiral
python -m experiments.exp_force3d_jones_toy
python -m experiments.exp_force3d_jones_action
python -m experiments.exp_force3d_read_k
python -m experiments.exp_force3d_updown
python -m experiments.exp_structure_gate_proof
python -m experiments.exp_structure_gate_minflux
python -m experiments.exp_algebra_A_firstorder
python -m experiments.exp_algebra_A_bimodule
python -m experiments.exp_algebra_A_edgespace
python -m experiments.exp_spin_connection_lattice
python -m experiments.exp_spin_connection_localJ
python -m experiments.exp_pi_flux_3D
python -m experiments.exp_pi_flux_3D_irreps
python -m experiments.exp_pi_flux_3D_localize
python -m experiments.exp_pi_flux_3D_first_wall
python -m experiments.exp_pi_flux_3D_curvature
python -m experiments.exp_pi_flux_3D_separation
```

依赖：`numpy>=1.24`，`scipy>=1.10`。结果写入对应 `experiments/<name>_last_run.json`。

> **Windows 注意**：若 `python` 命令指向 Windows Store 别名（无输出），改用 `py -m`。含中文/Unicode 输出的脚本（如 `exp7_k2_chiral`、`exp_force3d_jones_*`）在 GBK 控制台下可能报 `UnicodeEncodeError`，先在 PowerShell 设置环境变量 `PYTHONIOENCODING=utf-8` 再跑。

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
| 纯迹 Tr(D⁴) 压出近邻环 | ❌ 负结果（度约束下偏好 dimer；见纯迹作用量选维） |
| 纯迹冷启动自发选二维+π磁通 | ✅ 钉死：暖启动 5/5 稳定、basinhopping 全局最优 S=−71.1 |
| π 磁通正方格点陈数拓扑 | ❌ 负结果：陈数恒 0（dₓ=2cos²kₓ≥0），须换 honeycomb/自旋 |
| 自指产生曲率/π磁通（割裂） | ✅ 数值：纯迹墨西哥帽（无曲率项）π 磁通 0.5，非设计曲率项 |
| 面积律（toroidal 是空间） | ✅ 数值：面积律 S~l^1.15（vs 随机图体积律） |
| 纯迹自发选三维 | ❌ 负结果：Tr(D⁶)/Tr(D⁸)/SU(2) 版都不选维，要 Hopf 拓扑/第三层 |
| 自指 D → SU(2)_k 弦网 | ⛔ 封存 |
| 证明 3+1 | out of scope |
