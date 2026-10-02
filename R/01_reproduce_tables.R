# =====================================================================
# 01_reproduce_tables.R
# Reproduces Table 3, Table 4 and Supplementary Table S2 of
# Punia et al. (2026) Theor Appl Genet 139:265 from the entry x
# environment means in data/diallel_means.csv (Supplementary Table S1).
#
# Run from the repository root:  Rscript R/01_reproduce_tables.R
# =====================================================================

source("R/griffing_engine.R")

r_reps   <- 3
error_df <- 216   # e(r-1)(g-1) = 2 x 2 x 54

# Pooled error mean squares (plot basis) from Table 3 of the paper.
# They come from the replication-level RCBD ANOVA pooled over
# environments and cannot be recomputed from entry means. If you have
# plot-level data, compute them per environment with
# aov(trait ~ Rep + Genotype) and pool the residual SS and df.
error_ms <- c(DF = 7.45, DM = 21.55, PH = 186.99, Tillers = 0.04,
              PL = 2.66, PG = 0.25, TW = 0.41, DFY = 101.10,
              GY = 2.66, HI = 3.78, Fe = 9.30, Zn = 7.63,
              Protein = 0.56)

traits <- c("DF","DM","PH","Tillers","PL","PG","TW",
            "DFY","GY","HI","Fe","Zn","Protein")

# The dataset is NOT distributed with this repository. Build
# data/diallel_means.csv from Supplementary Table S1 of the paper
# (one row per entry per environment; see data/input_format.md).
dat <- read.csv("data/diallel_means.csv", check.names = FALSE)
out <- build_tables(dat, traits, r_reps, error_ms, error_df)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(out$table3, "results/tables/Table3_combining_ability_pooled.csv")
write.csv(out$table4, "results/tables/Table4_GCA_effects_stability.csv")
write.csv(out$sca,    "results/tables/TableS2_SCA_effects_pooled.csv",
          row.names = FALSE)

cat("Tables written to results/tables/\n")
cat("\nTable 3 (first six traits):\n"); print(out$table3[, 1:6])
cat("\nTable 4 (first six traits):\n"); print(out$table4[, 1:6])
