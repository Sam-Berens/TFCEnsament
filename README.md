# TFCEnsament

*TFCEnsament* is a small, self-contained MATLAB toolbox that performs **Threshold-Free Cluster Enhancement (TFCE)**–based statistical inference for volumetric (3-D) data, such as fMRI or other neuro-imaging modalities.  
It implements one-sample, correlation, and general linear model (GLM) contrasts with non-parametric null-distribution generation via bootstrapping or permutation.

---

## Features

| Function (file) | What it does |
|-----------------|--------------|
| **`tfce_getStat.m`** | Converts any voxel-wise test-statistic map (_t_, _F_, etc.) into its TFCE score. |
| **`tfce_Xcon.m`** | Computes GLM contrasts and returns both the raw statistic map and the TFCE-enhanced version. |
| **`tfce_nullBoot_oneSample.m`** | One-sample test against *μ = h₀* (default 0) using sign-flip bootstrapping. |
| **`tfce_nullBoot_corr.m`** | Tests a univariate correlation (*x* vs. each voxel) with sign-flip bootstrapping. |
| **`tfce_nullBoot_permuteRows.m`** | Arbitrary GLM contrast with a row-permutation bootstrap of _**X**_. |

All functions default to **10 000** bootstrap iterations and can run in parallel (`UseParFor = true`) when the Parallel Computing Toolbox is available.

---

## Installation

```bash
git clone https://github.com/<your-username>/TFCEnsament.git
cd TFCEnsament
# Inside MATLAB:
addpath(genpath(pwd))
savepath          # (optional) make it permanent
```

---

## Requirements
 - MATLAB R2018b or newer.
 - Statistics & Machine Learning Toolbox.
 - (Optional) Parallel Computing Toolbox — enables faster bootstrapping.
 - (Optional) SPM12 — only needed if you want to pass NIfTI filenames directly instead of 4-D matrices.

---

## Quick start
```MATLAB
%% One-sample test (H0: mean = 0)
[pVal, tfceStat] = tfce_nullBoot_oneSample(Y);   % Y is X×Y×Z×N (double)

%% Voxel-wise correlation
[pVal, tfceStat] = tfce_nullBoot_corr(Y, predictor);  % predictor is [N×1]

%% Arbitrary GLM contrast
[pVal, tfceStat] = tfce_nullBoot_permuteRows(Y, X, H);  % X is [N×p], H is [1×p]
```
Each call returns:
 - tfceStat — TFCE value at every voxel.
 - pVal — voxel-wise, family-wise-error–controlled p-values.
Masking or cluster-forming thresholds are never required — that’s the point of TFCE!

---

## How it works
 1. GLM / test statistic: <code>tfce_Xcon.m</code> fits OLS, computes a _t_- or _F_-map, then calls <code>tfce_getStat.m</code>.
 2. TFCE conversion: For each height threshold h (step = 0.1), voxels above _h_ are clustered (26-neighbour connectivity). Cluster size _E_ and height _H_ are integrated with the exponents _E^0.5 × H²_ and summed over _h_ (Smith & Nichols, 2009).
 3. Null distribution: Depending on the design, rows of _Y_ or _X_ are sign-flipped / permuted 10 000 times, re-computing TFCE each time.
 4. Family-wise error: Observed TFCE scores are compared against the null to yield voxel-wise _p_-values.

---

## Performance tips
 - Parallelise: Leave UseParFor at its default true if you have the Parallel Computing Toolbox.
 - Fewer iterations: For exploratory work, reduce nBoot inside the functions, but note that accuracy suffers.
 - Memory: Bootstrapping stores TFCE vectors only (not whole volumes) to keep RAM usage low.

## Citing this toolbox
If you use TFCEnsament in your work, please cite:

> Berens, S. (2025). TFCEnsament: MATLAB tools for threshold-free cluster enhancement. GitHub repository. doi:10.5281/zenodo.xxxxx

and the original TFCE paper:

> Smith, S.M., & Nichols, T.E. (2009). Threshold-Free Cluster Enhancement: Addressing problems of smoothing, threshold dependence and localisation in cluster inference. NeuroImage, 44(1), 83-98.

**Sources used:** :contentReference[oaicite:0]{index=0}
