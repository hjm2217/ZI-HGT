---
title: ZI-HGT
layout: default
---

# The Zero-Inflated Hierarchical Generalized Transformation Model (ZI-HGT)

![The ZI-HGT is a Bayesian auxiliary technique to noisily transform spatial transcriptomics data to approximate normality prior to analysis and to enable uncertainty quantification even when the analysis technique does not. We apply the ZI-HGT to CARD ([Ma et al, 2022](https://www.nature.com/articles/s41587-022-01273-7)) to perform spatially-informed cell-type deconvolution with data that better matches CARD's normal assumption and conduct uncertainty quantification.
In the ZI-HGT, the raw spatial transcriptomic data $X$ is input, generating $C$ posterior replicates ($H^{[1]}$, ..., $H^{[C]}$) that are interpreted as transformations of the original data that better follow model assumptions.  CARD is then applied to each transformed dataset and a corresponding single-cell RNA sequencing (scRNA-seq) dataset $B$, yielding cell-type proportion estimates ($\widehat{V}^{[1]}$, ..., $\widehat{V}^{[C]}$). These estimates are amalgamated into a single estimate $\widehat{V}$.](Figure1.pdf)

See our tutorial for running the ZI-HGT with CARD [here](Tutorial.html)!
