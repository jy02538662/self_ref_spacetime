# Exp7 Design: Spontaneous Clustering / SU(2) from D

> **Status**: design draft (code not written) · 2026-09-01  
> **Place**: strong-proposition probe only — core claims do **not** depend on Exp7  
> **Math**: GitHub-safe — Unicode in prose/tables; rare `$$` blocks only  
> **Related**: vault「字典观测点」八号种子

---

## Goal

Do not preset spins or cluster partitions. Let N indistinguishable pure states, under action minimization, spontaneously form “k states per cluster”. Core check: does **k=2** (SU(2) spin) spontaneously dominate?

This would be **numerical evidence** only — not an analytic proof that internal symmetry is an output of D.

---

## Action skeleton

Keep Exp1’s skeleton:

$$
S[D]=\mathrm{Tr}(D^{2})+\lambda_{1}\sum_{i \lt j}\bigl(d_{ij}-L_{ij}\bigr)^{2}
$$

- Tr(D²): coupling cost  
- Geometric term: self-ref distances d_ij match targets L_ij  

Hard point: define L_ij **without** presupposing clusters.

### Scheme 1 — implementable, weakly presupposes k=2

- N=2M states; target geometry = M points (ring or grid) with L^geom_ab  
- Random partition into M blocks of size 2; within block L_ij=0; between blocks use L^geom  
- After min: check whether D keeps 2-state clusters (intra |z| ≫ inter)

Pros: tests stability of “2-state cluster”. Cons: presupposes block size 2.

### Scheme 2 — unpresupposed (needs exploration)

- L_ij is a rank-M lift of an M-point geometry; assignment organized by D  
- Soft assignment / spectral clustering each step from current d_ij → update L  
- Check whether cluster-size histogram converges to k=2  

Pros: true spontaneous-k test. Cons: implementation details open.

---

## Clustering diagnostics (post-min)

1. Intra / inter |z| ratio > 10³  
2. Cluster-size histogram (k-means / spectral); does size 2 dominate?  
3. Silhouette score > 0.9  
4. After Scheme 1: Pauli expand each 2×2 block in {I, σx, σy, σz}; check tunable spin direction  

---

## Verdict table

| Outcome | Verdict |
|---------|---------|
| k=2 most frequent and lowest energy | SU(2) spontaneous — numeric support |
| k=1 dominates | Classical k=1 cheaper; need ZPE / fluctuation |
| k≥3 dominates | Strong prop fails (or needs topo constraint) |
| Sizes random / no converge | Action insufficient |

---

## Open points

1. Scheme 2 soft-assignment iteration  
2. ZPE-style correction ~ k(k−1)−c√k (√k scaling **unverified**)  
3. Topological reason for k=2 (old “k≥3 fractionalizes” argument was judged **wrong**)  

---

## Honest placement

- Numeric evidence only; “why spin-1/2” is deeper than the minimal model must answer  
- Core repo contribution (same D → distance + signature + topology) stands **without** Exp7  

---

## Code skeleton

```
experiments/exp7_clustering.py
  - Scheme 1: N=2M + random block L + minimize_action (src/optimize.py)
  - Clustering: sklearn KMeans / SpectralClustering + silhouette
  - Pauli: expand 2×2 blocks
```

Deps: `requirements.txt` + `scikit-learn`. No `*_last_run.json` until implemented.
