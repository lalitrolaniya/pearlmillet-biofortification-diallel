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

## The app

**https://lalitrolaniya.github.io/pearlmillet-biofortification-diallel/**

Runs entirely in the browser (no R installation needed). Upload replicated
half-diallel data and get the complete combining ability analysis pooled over
environments, following Griffing (1956) Method 2, Model 1:

- **ANOVA and genetic parameters**: GCA, SCA, GCA x E, SCA x E and pooled
  error mean squares with F-tests; variance components; Baker's ratio;
  average degree of dominance; percent contributions of GCA and SCA
- **GCA effects and stability**: effects with significance, SE(gi),
  SE(gi-gj), and GCA x E stability variance per parent
- **SCA effects**: all crosses with significance
- **GCA plots**: bar plots of GCA effects per trait, downloadable as PNG

All result tables can be downloaded as CSV. Any number of parents,
environments, replications and traits is supported.

## Data format

One CSV file, one row per plot:

| Genotype | P1 | P2 | Env | Rep | GY | Fe | ... |
|----------|----|----|-----|-----|------|------|-----|
| P1 | P1 | P1 | E1 | 1 | 19.89 | 45.58 | ... |
| P1 x P2 | P1 | P2 | E1 | 1 | 16.45 | 51.22 | ... |
| P1 x P2 | P1 | P2 | E1 | 2 | 15.80 | 49.75 | ... |
| P1 x P2 | P1 | P2 | E2 | 1 | 14.91 | 52.10 | ... |

- `Genotype`: entry name (for a parent, its own name; for a cross, any label)
- `P1`, `P2`: the two parental lines (`P1` = `P2` for a parent entry)
- `Env`: environment label (E1, E2, ...)
- `Rep`: replication number within the environment
- then one column per trait, any number of traits, any names

A ready example is in `data/example_replicated_data.csv` (4 parents,
2 environments, 3 replications, 2 traits) and can also be downloaded from
inside the app. Full details are in `data/input_format.md`.

The app also accepts entry x environment means (same columns without `Rep`);
significance of effects then needs a separate CSV of pooled error mean
squares (`Trait,MS`).

## Statistical methods

GCA is tested against GCA x E, SCA against SCA x E, and the interactions
against the pooled error from the replication-level RCBD ANOVA pooled over
environments. Variance components follow eqs. 6-8 of the paper:

- sigma2_GCA = (MS_GCA - MS_GCAxE) x p / [r x e x (p + 2)]
- sigma2_SCA = (MS_SCA - MS_SCAxE) / (r x e)
- Baker's ratio = 2 sigma2_GCA / (2 sigma2_GCA + sigma2_SCA)

The engine (`R/griffing_engine.R`) is base R with no package dependencies
and can be used directly in scripts.

## Run locally

```r
shiny::runApp("app")
# or without cloning:
shiny::runGitHub("pearlmillet-biofortification-diallel", "lalitrolaniya")
```

## Citing

Please cite the paper (above). Citation metadata for this repository is in
`CITATION.cff`. To cite this repository itself:
https://doi.org/10.5281/zenodo.23101654

## License

Code is released under the MIT License. Please cite the paper when using
this work.
