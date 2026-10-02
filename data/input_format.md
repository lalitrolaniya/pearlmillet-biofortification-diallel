# Input format - data/diallel_means.csv (not distributed)

The scripts expect entry x environment means (means over
replications): one row per entry per environment. The dataset of the
paper (55 genotypes: 10 parents + 45 F1 hybrids, half diallel without
reciprocals, RCBD with 3 replications, two sowing-date environments at
Jaipur) is published as Supplementary Table S1 on the article page and
is not redistributed here.

| Column   | Description |
|----------|-------------|
| Genotype | Parent name, or "P1 x P2" for an F1 cross |
| P1, P2   | Parental lines (P1 = P2 for a parent entry) |
| Env      | Environment: E1 (timely sowing), E2 (late sowing) |
| DF       | Days to 50% flowering |
| DM       | Days to maturity |
| PH       | Plant height (cm) |
| Tillers  | Productive tillers per plant |
| PL       | Panicle length (cm) |
| PG       | Panicle girth (cm) |
| TW       | Test weight (g) |
| DFY      | Dry fodder yield per plant (g) |
| GY       | Grain yield per plant (g) |
| HI       | Harvest index (%) |
| Fe       | Grain iron content (mg/kg) |
| Zn       | Grain zinc content (mg/kg) |
| Protein  | Grain protein content (%) |

Parents: RIB-9178, RIB-9184, RIB-9185, RIB-9205, RIB-15131,
RIB-16300, RIB-16324, J-2340, 20K86, RIB-192.
