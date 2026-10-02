# =====================================================================
# griffing_engine.R
# Griffing (1956) Method 2, Model 1 combining ability analysis,
# pooled over environments. Base R only.
#
# Companion code for:
#   Punia M, Sharma LD, Gothwal DK, Kajla SL, Rolaniya LK, Sharma V,
#   Jat RL (2026) Combining ability and gene action for grain yield
#   and biofortification traits in pearl millet. Theor Appl Genet
#   139:265. https://doi.org/10.1007/s00122-026-05378-4
#
# Key formulas (paper, Methods):
#   sigma2_GCA = (MS_GCA - MS_GCAxE) * p / (r * e * (p + 2))   [eq. 6]
#   sigma2_SCA = (MS_SCA - MS_SCAxE) / (r * e)                 [eq. 7]
#   Baker's ratio = 2 s2GCA / (2 s2GCA + s2SCA)                [eq. 8]
#   GCA tested against GCA x E; SCA against SCA x E;
#   interactions tested against the pooled error (df = e(r-1)(g-1)).
#   Effect SEs use MSe/(r e) as error variance:
#     Var(g_i)        = (p-1)/(p(p+2))            * MSe/(re)
#     Var(g_i - g_j)  = 2/(p+2)                   * MSe/(re)
#     Var(s_ij), i!=j = (p^2+p+2)/((p+1)(p+2))    * MSe/(re)
#     Var(s_ii)       = p(p-1)/((p+1)(p+2))       * MSe/(re)
# =====================================================================

# Griffing Method 2 on a p x p symmetric matrix of entry means
# (diagonal = parents). Returns effects and SS on the MEAN basis.
griffing_m2 <- function(M) {
  p    <- nrow(M)
  X..  <- sum(M[upper.tri(M, diag = TRUE)])
  Ti   <- rowSums(M) + diag(M)                 # X_i. + X_ii
  g    <- (Ti - 2 * X.. / p) / (p + 2)
  S    <- matrix(NA_real_, p, p, dimnames = dimnames(M))
  for (i in 1:p) for (j in i:p)
    S[i, j] <- M[i, j] - (Ti[i] + Ti[j]) / (p + 2) +
               2 * X.. / ((p + 1) * (p + 2))
  SSg <- (sum(Ti^2) - 4 * X..^2 / p) / (p + 2)
  SSs <- sum(M[upper.tri(M, diag = TRUE)]^2) -
         sum(Ti^2) / (p + 2) + 2 * X..^2 / ((p + 1) * (p + 2))
  list(g = g, s = S, SSg = SSg, SSs = SSs)
}

stars <- function(pval)
  ifelse(is.na(pval), "",
  ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**",
  ifelse(pval < 0.05, "*", ""))))

# Pooled Griffing analysis for one trait.
# dat: data.frame with columns Genotype, P1, P2, Env and the trait;
#      one row per entry per environment (entry means over reps).
# error_ms may be NA: ANOVA F-tests of GCA and SCA (vs their
# interactions) are still exact; stars needing the error MS are blank.
pooled_griffing <- function(dat, trait, r_reps, error_ms = NA,
                            error_df = NA) {
  parents <- unique(c(dat$P1, dat$P2))
  p  <- length(parents)
  es <- sort(unique(dat$Env)); e <- length(es)

  mats <- lapply(es, function(env) {
    M <- matrix(NA_real_, p, p, dimnames = list(parents, parents))
    d <- dat[dat$Env == env, ]
    for (k in seq_len(nrow(d))) {
      i <- match(d$P1[k], parents); j <- match(d$P2[k], parents)
      M[i, j] <- M[j, i] <- d[[trait]][k]
    }
    M
  })
  fit_env    <- lapply(mats, griffing_m2)
  M_pooled   <- Reduce(`+`, mats) / e
  fit_pooled <- griffing_m2(M_pooled)

  # Sums of squares on the plot basis
  SS_gca <- r_reps * e * fit_pooled$SSg
  SS_sca <- r_reps * e * fit_pooled$SSs
  SS_gxe <- r_reps * sum(sapply(fit_env, `[[`, "SSg")) - SS_gca
  SS_sxe <- r_reps * sum(sapply(fit_env, `[[`, "SSs")) - SS_sca
  df_g <- p - 1; df_s <- p * (p - 1) / 2
  df_ge <- df_g * (e - 1); df_se <- df_s * (e - 1)
  MS <- c(GCA = SS_gca/df_g, SCA = SS_sca/df_s,
          GCAxE = SS_gxe/df_ge, SCAxE = SS_sxe/df_se, Error = error_ms)

  pv <- c(GCA   = pf(MS[["GCA"]]/MS[["GCAxE"]], df_g, df_ge, lower.tail = FALSE),
          SCA   = pf(MS[["SCA"]]/MS[["SCAxE"]], df_s, df_se, lower.tail = FALSE),
          GCAxE = if (is.na(error_ms)) NA else
                  pf(MS[["GCAxE"]]/error_ms, df_ge, error_df, lower.tail = FALSE),
          SCAxE = if (is.na(error_ms)) NA else
                  pf(MS[["SCAxE"]]/error_ms, df_se, error_df, lower.tail = FALSE))

  # Variance components and genetic parameters (eqs. 6-8)
  s2gca <- max(0, (MS[["GCA"]] - MS[["GCAxE"]]) * p / (r_reps * e * (p + 2)))
  s2sca <- max(0, (MS[["SCA"]] - MS[["SCAxE"]]) / (r_reps * e))
  s2A <- 2 * s2gca; s2D <- s2sca
  baker <- 2 * s2gca / (2 * s2gca + s2sca)
  add   <- sqrt(2 * s2D / s2A)
  pctG  <- 100 * SS_gca / (SS_gca + SS_sca)

  # GCA effects, SEs, t-tests, GCA x E stability
  g <- fit_pooled$g
  v_unit <- if (is.na(error_ms)) NA else error_ms / (r_reps * e)
  se_g   <- sqrt((p - 1) / (p * (p + 2)) * v_unit)
  se_gg  <- sqrt(2 / (p + 2)             * v_unit)
  p_g    <- if (is.na(error_ms)) rep(NA, p) else
            2 * pt(abs(g / se_g), error_df, lower.tail = FALSE)
  g_env  <- sapply(fit_env, `[[`, "g")
  stab   <- apply(g_env, 1, var)        # denominator (e - 1)

  # SCA effects with t-tests (variances derived from Griffing M2)
  sca   <- fit_pooled$s
  se_sij <- sqrt((p^2 + p + 2) / ((p + 1) * (p + 2)) * v_unit)
  se_sii <- sqrt(p * (p - 1)   / ((p + 1) * (p + 2)) * v_unit)
  p_sca <- sca
  for (i in 1:p) for (j in i:p) {
    se <- if (i == j) se_sii else se_sij
    p_sca[i, j] <- if (is.na(error_ms)) NA else
                   2 * pt(abs(sca[i, j] / se), error_df, lower.tail = FALSE)
  }

  list(parents = parents, MS = MS, pval = pv,
       s2gca = s2gca, s2sca = s2sca, s2A = s2A, s2D = s2D,
       baker = baker, ratio = s2gca / s2sca, add = add,
       pctG = pctG, pctS = 100 - pctG,
       g = g, p_g = p_g, se_g = se_g, se_gg = se_gg, stab = stab,
       sca = sca, p_sca = p_sca,
       df = c(GCA = df_g, SCA = df_s, GCAxE = df_ge, SCAxE = df_se,
              Error = error_df))
}

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# Assemble the three output tables for a set of traits.
build_tables <- function(dat, traits, r_reps, error_ms, error_df) {
  res <- lapply(traits, function(tr)
    pooled_griffing(dat, tr, r_reps,
                    if (tr %in% names(error_ms)) error_ms[[tr]] else NA,
                    error_df))
  names(res) <- traits
  parents <- res[[1]]$parents

  t3 <- data.frame(row.names = c(
    "GCA", "SCA", "GCA x E", "SCA x E", "Error",
    "sigma2 GCA", "sigma2 SCA", "sigma2 A (additive)",
    "sigma2 D (dominance)", "Baker's ratio", "s2GCA/s2SCA",
    "Avg. degree of dominance", "% Contribution GCA",
    "% Contribution SCA"))
  for (tr in traits) {
    x <- res[[tr]]
    t3[[tr]] <- c(paste0(fmt(x$MS[1:4]), stars(x$pval)),
                  ifelse(is.na(x$MS[5]), "-", fmt(x$MS[5])),
                  fmt(c(x$s2gca, x$s2sca, x$s2A, x$s2D)),
                  fmt(c(x$baker, x$ratio, x$add)),
                  fmt(c(x$pctG, x$pctS)))
  }

  t4 <- data.frame(row.names = c(paste0("GCA_", parents),
                                 "SE(gi)", "SE(gi-gj)",
                                 paste0("Stability_", parents)))
  for (tr in traits) {
    x <- res[[tr]]
    t4[[tr]] <- c(paste0(fmt(x$g), stars(x$p_g)),
                  ifelse(is.na(x$se_g), "-", fmt(x$se_g, 3)),
                  ifelse(is.na(x$se_gg), "-", fmt(x$se_gg, 3)),
                  fmt(x$stab, 3))
  }

  crosses <- which(upper.tri(res[[1]]$sca, diag = FALSE), arr.ind = TRUE)
  ts2 <- data.frame(Cross = paste(parents[crosses[, 1]], "x",
                                  parents[crosses[, 2]]))
  for (tr in traits) {
    x <- res[[tr]]
    ts2[[tr]] <- paste0(fmt(x$sca[crosses]), stars(x$p_sca[crosses]))
  }

  list(res = res, table3 = t3, table4 = t4, sca = ts2)
}

# ---------------------------------------------------------------------
# Replicated (plot-level) data preparation.
# plot_dat: Genotype, P1, P2, Env, Rep, then one column per trait
# (one row per plot). For each environment an RCBD ANOVA
# (trait ~ Rep + Genotype) is fitted; residual SS and df are pooled
# over environments to give the pooled error MS (df = e(r-1)(g-1)).
# Returns entry x environment means, error MS per trait, error df,
# and the number of replications.
# ---------------------------------------------------------------------
prep_from_plot_data <- function(plot_dat, traits = NULL) {
  id <- c("Genotype", "P1", "P2", "Env", "Rep")
  if (is.null(traits))
    traits <- names(plot_dat)[!(names(plot_dat) %in% id) &
                              sapply(plot_dat, is.numeric)]
  ems <- numeric(0); edf <- NA
  for (tr in traits) {
    sse <- 0; dfe <- 0
    for (env in unique(plot_dat$Env)) {
      d <- plot_dat[plot_dat$Env == env & !is.na(plot_dat[[tr]]), ]
      fit <- aov(d[[tr]] ~ factor(d$Rep) + factor(d$Genotype))
      an  <- anova(fit)
      sse <- sse + an["Residuals", "Sum Sq"]
      dfe <- dfe + an["Residuals", "Df"]
    }
    ems[tr] <- sse / dfe
    edf <- dfe
  }
  means <- aggregate(plot_dat[traits],
                     by = plot_dat[c("Genotype", "P1", "P2", "Env")],
                     FUN = mean, na.rm = TRUE)
  list(means = means, error_ms = ems, error_df = edf,
       r = length(unique(plot_dat$Rep)), traits = traits)
}
