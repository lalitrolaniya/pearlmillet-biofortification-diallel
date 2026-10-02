# Input format

## Replicated plot-level data (recommended)

One CSV, one row per plot:

| Column   | Description |
|----------|-------------|
| Genotype | Entry name: a parent's own name, or a label like "P1 x P2" for a cross |
| P1, P2   | The two parental lines (P1 = P2 for a parent entry) |
| Env      | Environment label (E1, E2, ...) |
| Rep      | Replication number within the environment |
| traits   | One column per trait, any number, any names |

The design must be a half diallel: all p parents plus the p(p-1)/2 F1
crosses, no reciprocals, evaluated in an RCBD in every environment. The
app computes the pooled error itself from this file; nothing else is
needed. See `example_replicated_data.csv` in this folder.

## Entry x environment means (fallback)

Same columns without `Rep`, one row per entry per environment (values are
means over replications). Significance of effects then requires the pooled
error mean squares in a separate CSV with columns `Trait,MS`, plus the
error degrees of freedom and replication number entered in the app.

## The dataset of the paper

Entry x environment means for the 55 pearl millet genotypes are published
as Supplementary Table S1 on the article page
(https://doi.org/10.1007/s00122-026-05378-4) and are not redistributed
here.
