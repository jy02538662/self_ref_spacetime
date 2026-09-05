# Exp7 Design: Spontaneous Clustering / SU(2) from $D$

> **Status**: design draft (code not written) · 2026-09-01  
> **Place**: strong-proposition probe only — core claims of this repo do **not** depend on Exp7  
> **Related**: vault「字典观测点」八号种子

---

## Goal

Do **not** preset spins or cluster partitions. Let $N$ indistinguishable pure states, under action minimization, spontaneously form “$k$ states per cluster”. The core check:

> Does **$k=2$** (SU(2) spin) spontaneously dominate?

This would be **numerical evidence** for “internal symmetry is an output of $D$, not an input” — not an analytic proof.

---

## Action skeleton

Keep Exp1’s self-referential skeleton:

$$
S[D]=\mathrm{Tr}(D^{2})+\lambda_{1}\sum_{i<j}\bigl(d_{ij}-L_{ij}\bigr)^{2}
$$

- $\mathrm{Tr}(D^{2})$: coupling cost (favors fewer / weaker links)
- Geometric term: self-ref distances $d_{ij}$ should match targets $L_{ij}$

**Hard point**: define $L_{ij}$ **without** presupposing a cluster partition.

### Scheme 1 — implementable, weakly presupposes $k=2$

- Take $N=2M$ states; target geometry is $M$ points (1D ring or 2D grid) with distances $L^{\mathrm{geom}}_{ab}$
- Randomly partition the $N$ states into $M$ blocks of size 2; within a block $L_{ij}=0$; between blocks $L_{ij}=L^{\mathrm{geom}}_{c_i c_j}$
- After minimization, check whether $D$ keeps 2-state clusters (intra-cluster $|z|$ ≫ inter-cluster)

**Pros**: tests whether “2-state cluster = SU(2)” is a stable solution of $D$.  
**Cons**: presupposes block size 2 — not true spontaneous $k=2$.

### Scheme 2 — truly unpresupposed (needs exploration)

- $L_{ij}$ is a rank-$M$ low-rank lift of an $M$-point geometry, but the lift (which states share a cluster) is organized by $D$
- Soft assignment / OT: cluster labels $c_i$ from current $d_{ij}$ each step (e.g. spectral clustering), then $L_{ij}=L^{\mathrm{geom}}_{c_i c_j}$
- Iterate; check whether cluster-size distribution converges to $k=2$

**Pros**: genuine spontaneous-$k$ test.  
**Cons**: implementation details of soft assignment need exploration.

---

## Clustering diagnostics (post-min)

On $d_{ij}$ or $|D|$:

1. **Intra / inter coupling ratio**:
   $$\frac{\mathrm{mean}\,|z|_{\mathrm{intra}}}{\mathrm{mean}\,|z|_{\mathrm{inter}}}>10^{3}$$
2. **Cluster-size histogram**: $k$-means or spectral clustering on $d_{ij}$; does size $k=2$ dominate?
3. **Silhouette score** $>0.9$
4. **Pauli decomposition** (after Scheme 1 holds): for each $2\times 2$ reduced block, expand in $\{I,\sigma_x,\sigma_y,\sigma_z\}$ and check whether the three $\sigma$ components are continuously tunable (spin direction)

---

## Verdict table

| Outcome | Verdict |
|---------|---------|
| $k=2$ most frequent and lowest energy | ✅ SU(2) spontaneous; strong prop has numeric support |
| $k=1$ (no internal structure) dominates | Classical $k=1$ cheaper; need ZPE / fluctuation correction |
| $k\ge 3$ dominates | Strong prop fails (or needs topo constraint) |
| Sizes do not converge / look random | Action insufficient; add constraints |

---

## Open points

1. Scheme 2 soft-assignment / spectral iteration
2. Quantum fluctuation / ZPE correction of the form $E\propto k(k-1)-c\sqrt{k}$ (note: $\sqrt{k}$ scaling **unverified**)
3. Topological constraint making $k=2$ the minimal irreducible unit (the old “$k\ge 3$ topo charge fractionalizes” argument was judged **wrong** — need another reason)

---

## Honest placement

- Exp7 can only deliver **numeric evidence**, not an analytic proof of “why $k=2$”
- “Why spin-$1/2$” is a deep foundational question; the minimal model need not answer it
- Core contribution of this repo (same $D$ → distance + signature + topology) **already stands without** the strong proposition

---

## Code skeleton (reuse `src/`)

```
experiments/exp7_clustering.py
  - Scheme 1: N=2M + random block L + minimize_action (src/optimize.py)
  - Clustering: sklearn KMeans / SpectralClustering + silhouette_score
  - Pauli: expand each 2×2 block in I, σx, σy, σz
```

Deps: existing `requirements.txt` (`numpy`, `scipy`) + `scikit-learn` for clustering.
