# Pearl Millet Biofortification Diallel
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23101654.svg)](https://doi.org/10.5281/zenodo.23101654)

Companion code for:

> Punia M, Sharma LD, Gothwal DK, Kajla SL, Rolaniya LK, Sharma V, Jat RL (2026)
> **Combining ability and gene action for grain yield and biofortification traits
> in pearl millet [*Pennisetum glaucum* (L.) R. Br.]: implications for breeding
> high-yielding biofortified hybrids in arid regions.**
> *Theoretical and Applied Genetics* 139:265.
> https://doi.org/10.1007/s00122-026-05378-4

*This repository is dedicated to the memory of Dr. Monika Punia, the first
author, who conceived and led this research as her doctoral work. The paper is
published in her memory, in fulfilment of her wish to see it in print.*

## What is here

| Path | Contents |
|------|----------|
| `R/griffing_engine.R` | Griffing (1956) Method 2, Model 1 combining ability analysis pooled over environments. Base R, no packages |
| `R/01_reproduce_tables.R` | Reproduces Table 3, Table 4 and Supplementary Table S2 of the paper from the Table S1 means |
| `R/02_figures.R` | Example figures (GCA effects, Fe-Zn association) |
| `data/input_format.md` | Expected input format, trait codes and design details |
| `figures_supplementary/` | Supplementary Figures S1-S9 as published, with captions |
| `app/` | Shiny app: upload replicated half-diallel data (or entry means) and get the full pooled analysis |

The dataset itself is not redistributed here. The entry x environment means
are published as Supplementary Table S1 on the article page; place them in
`data/diallel_means.csv` (format in `data/input_format.md`) to run the
reproduction scripts.

## Quick start

```r
# from the repository root, after preparing data/diallel_means.csv
source("R/01_reproduce_tables.R")   # Tables 3, 4 and S2 -> results/tables/
source("R/02_figures.R")            # figures -> results/figures/

# interactive app (upload your own CSV)
shiny::runApp("app")
# or without cloning:
shiny::runGitHub("pearlmillet-biofortification-diallel", "lalitrolaniya")
```

A browser version of the app (no R needed) is deployed with shinylive at:
`https://lalitrolaniya.github.io/pearlmillet-biofortification-diallel/`

## Statistical methods

Griffing (1956) Method 2, Model 1 for a half diallel (parents + F1s, no
reciprocals), pooled over environments (Singh and Chaudhary 1979). GCA is
tested against GCA x E, SCA against SCA x E, and the interactions against the
pooled error (df = 216). Variance components follow eqs. 6-8 of the paper:

- sigma2_GCA = (MS_GCA - MS_GCAxE) x p / [r x e x (p + 2)]
- sigma2_SCA = (MS_SCA - MS_SCAxE) / (r x e)
- Baker's ratio = 2 sigma2_GCA / (2 sigma2_GCA + sigma2_SCA)

## Reproducibility

Running `R/01_reproduce_tables.R` on the Table S1 means reproduces the
published Tables 3 and 4: GCA effects, standard errors, significance and
GCA x E stability match exactly, and all mean squares and Baker's ratios
agree within 0.5%, the residual being due to the two-decimal rounding of the
published means. The pooled error mean squares (which require
replication-level data) are taken from Table 3 of the paper.

## Using the app with your own data

Preferred input is replicated plot-level data: columns `Genotype, P1, P2,
Env, Rep` followed by one column per trait, one row per plot. The app then
computes the pooled error itself and gives the complete analysis with
significance tests. Entry x environment means (same columns without `Rep`)
are also accepted; significance of effects then needs a separate CSV of
pooled error mean squares (`Trait,MS`). Any number of parents,
environments, replications and traits is supported.

## Citing

Please cite the paper (above). Citation metadata for this repository is in
`CITATION.cff`.
To cite this repository itself: https://doi.org/10.5281/zenodo.23101654

## License

Code is released under the MIT License. The supplementary figures are from
the published article's electronic supplementary material; please cite the
paper when using them.
