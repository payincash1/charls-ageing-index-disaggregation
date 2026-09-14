# 08_single_marker.R — single-marker + single-dimension benchmarks (the ladder)
rm(list = ls()); library(dplyr); library(survival)
base <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT
load(file.path(base, "output/metrics_final.RData"))
X8 <- read.csv(file.path(base, "processed/charls_std_for_vae.csv"))
ids8 <- read.csv(file.path(base, "processed/charls_ids.csv"), colClasses = c(ID = "character"))
X8$ID <- ids8$ID
d <- analytic_df %>% left_join(X8, by = "ID")
c_base <- as.numeric(summary(cox_base)$concordance["C"])

row_one <- function(v) {
  f <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +", v)), data = d, cluster = communityID)
  s <- summary(f)$coefficients
  data.frame(model = v, HR = round(exp(s[v,"coef"]), 3), C = round(as.numeric(summary(f)$concordance["C"]), 4),
             dC = round(as.numeric(summary(f)$concordance["C"]) - c_base, 4)) }
cands <- c("cystatin_c_rev","crp_rev","hba1c_rev","total_chol_rev",
           "physio_proxy","tics_cog","cesd_rev","self_health_rev")
res <- bind_rows(lapply(cands, row_one)) %>% arrange(desc(dC))
res <- bind_rows(
  data.frame(model = "dims4", C = round(as.numeric(summary(cox_dims)$concordance["C"]), 4),
             dC = round(as.numeric(summary(cox_dims)$concordance["C"]) - c_base, 4), HR = NA),
  data.frame(model = "MPAI", C = round(as.numeric(summary(cox_mpai)$concordance["C"]), 4),
             dC = round(as.numeric(summary(cox_mpai)$concordance["C"]) - c_base, 4), HR = NA),
  res)
print(res)
write.csv(res, file.path(base, "output/single_marker_benchmark.csv"), row.names = FALSE)

# dimension-alone ΔC
for (v in c("z_renal","z_funcinfl","z_metabolic","z_central")) {
  f <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +", v)), data = d, cluster = communityID)
  cat(v, "alone: dC =", round(as.numeric(summary(f)$concordance["C"]) - c_base, 4), "\n") }

# Cox-Snell R2 for scalar comparators
d$eq8 <- rowMeans(d[, c("tics_cog","self_health_rev","cesd_rev","physio_proxy",
                        "hba1c_rev","crp_rev","cystatin_c_rev","total_chol_rev")])
pc <- prcomp(d[, c("tics_cog","self_health_rev","cesd_rev","physio_proxy",
                   "hba1c_rev","crp_rev","cystatin_c_rev","total_chol_rev")])
d$pc1 <- pc$x[, 1]
f_eq8 <- coxph(Surv(surv_time, died) ~ age + sex + eq8, data = d, x = TRUE)
f_pc1 <- coxph(Surv(surv_time, died) ~ age + sex + pc1, data = d, x = TRUE)
n <- nrow(d)
cs_r2 <- function(fit) 1 - exp(-2 * (fit$loglik[2] - fit$loglik[1]) / n)
cat("PC1 R2 =", round(cs_r2(f_pc1), 4), " eq8 R2 =", round(cs_r2(f_eq8), 4), "\n")

# loader-set synergy benchmark
f_syn <- coxph(Surv(surv_time, died) ~ age + sex + cystatin_c_rev +
                 self_health_rev + tics_cog + crp_rev, data = d, cluster = communityID)
cat("loader-set dC:", round(as.numeric(summary(f_syn)$concordance["C"]) - c_base, 4),
    " (renal-dim alone 0.0163)\n")
