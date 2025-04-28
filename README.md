# TFCEnsament

*TFCEnsament* is a small, self-contained MATLAB toolbox that performs **Threshold-Free Cluster Enhancement (TFCE)**–based statistical inference for volumetric (3-D) data, such as fMRI or other neuro-imaging modalities.  
It implements one-sample, correlation, and general linear-model (GLM) contrasts with non-parametric null-distribution generation via bootstrapping or permutation.

<p align="center">
  <img src="https://user-images.githubusercontent.com/placeholder/tfce-logo.png" width="320" alt="TFCE illustration">
</p>

---

## Features

| Function (file) | What it does |
|-----------------|--------------|
| **`tfce_getStat.m`** | Converts any voxel-wise test-statistic map (t, F, etc.) into its TFCE score. |
| **`tfce_Xcon.m`** | Computes GLM contrasts and returns both the raw statistic map and the TFCE-enhanced version. |
| **`tfce_nullBoot_oneSample.m`** | One-sample test against *μ = h₀* (default 0) using sign-flip bootstrapping. |
| **`tfce_nullBoot_corr.m`** | Tests a univariate correlation (*x* vs. each voxel) with sign-flip bootstrapping. |
| **`tfce_nullBoot_permuteRows.m`** | Arbitrary GLM contrast with a row-permutation bootstrap of **X**. |

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
