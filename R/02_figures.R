# =====================================================================
# 02_figures.R
# Example figures from the diallel means data.
# Run from the repository root:  Rscript R/02_figures.R
# =====================================================================

source("R/griffing_engine.R")

dat    <- read.csv("data/diallel_means.csv", check.names = FALSE)
r_reps <- 3
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

## Figure 1: GCA effects of the 10 parents for GY, Fe, Zn, Protein ----
traits4 <- c(GY = "Grain yield (g/plant)", Fe = "Iron (mg/kg)",
             Zn = "Zinc (mg/kg)", Protein = "Protein (%)")
fits <- lapply(names(traits4), function(tr) pooled_griffing(dat, tr, r_reps))
names(fits) <- names(traits4)
parents <- fits[[1]]$parents

png("results/figures/Fig_GCA_effects.png", 2000, 1500, res = 220)
par(mfrow = c(2, 2), mar = c(7, 4, 2.5, 1), mgp = c(2.4, 0.7, 0))
for (tr in names(traits4)) {
  g <- fits[[tr]]$g
  barplot(g, names.arg = parents, las = 2, cex.names = 0.8,
          col = ifelse(g > 0, "#1F4023", "#B98A1E"), border = NA,
          main = traits4[[tr]], ylab = "GCA effect")
  abline(h = 0)
}
dev.off()

## Figure 2: Fe vs Zn, pooled means of all 55 entries ----------------
pool <- aggregate(dat[c("Fe", "Zn", "GY")],
                  by = dat[c("Genotype", "P1", "P2")], FUN = mean)
is_parent <- pool$P1 == pool$P2
best <- pool$Genotype == "RIB-9184 x RIB-15131"

png("results/figures/Fig_Fe_Zn_scatter.png", 1800, 1500, res = 220)
par(mar = c(4, 4, 2.5, 1))
plot(pool$Fe, pool$Zn, pch = ifelse(is_parent, 17, 16),
     col = ifelse(best, "#8C3B2E",
           ifelse(is_parent, "#B98A1E", "#1F402380")),
     cex = ifelse(best, 1.6, 1.1),
     xlab = "Grain Fe (mg/kg)", ylab = "Grain Zn (mg/kg)",
     main = "Fe-Zn association across 55 entries (pooled means)")
abline(lm(Zn ~ Fe, data = pool), lty = 2, col = "grey40")
text(pool$Fe[best], pool$Zn[best], "RIB-9184 x RIB-15131",
     pos = 2, col = "#8C3B2E", cex = 0.85)
legend("topleft", pch = c(16, 17, 16), bty = "n", cex = 0.85,
       col = c("#1F4023", "#B98A1E", "#8C3B2E"),
       legend = c("F1 hybrids", "Parents", "Best hybrid"))
r <- cor(pool$Fe, pool$Zn)
mtext(sprintf("r = %.2f", r), side = 3, adj = 1, cex = 0.9)
dev.off()

cat("Figures written to results/figures/\n")
