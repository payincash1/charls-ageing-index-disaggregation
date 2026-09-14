# 14_second_round_analyses.R —— 分析A（解释变异+替代标量）+ C（胱抑素C缺失电池）+ D（维度-年龄相关）

rm(list = ls()); library(dplyr); library(haven); library(survival)
base <- "/Users/a1234/Desktop/桌面/chalrs"
load(file.path(base, "output/metrics_final.RData"))

X8 <- read.csv(file.path(base, "processed/charls_std_for_vae.csv"))
ids8 <- read.csv(file.path(base, "processed/charls_ids.csv"), colClasses = c(ID = "character"))
X8$ID <- ids8$ID
d <- analytic_df %>% left_join(X8, by = "ID")
z8 <- c("tics_cog","self_health_rev","cesd_rev","physio_proxy","hba1c_rev","crp_rev","cystatin_c_rev","total_chol_rev")
c_base <- as.numeric(summary(cox_base)$concordance["C"])
getC <- function(f) as.numeric(summary(f)$concordance["C"])
dC <- function(f) round(getC(f) - c_base, 4)

# ========== A1. Cox-Snell explained variation ==========
n <- cox_base$n
cs_r2 <- function(fit) 1 - exp(-2 * (fit$loglik[2] - fit$loglik[1]) / n)
r_base <- cs_r2(cox_base); r_mpai <- cs_r2(cox_mpai); r_dims <- cs_r2(cox_dims)
r_raw <- cs_r2(cox_raw8); r_w8 <- cs_r2(cox_w8)
cat("=== A1 Cox-Snell R2 ===\n")
cat("base", round(r_base, 4), " MPAI", round(r_mpai, 4), " dims4", round(r_dims, 4),
    " raw8", round(r_raw, 4), " weighted8", round(r_w8, 4), "\n")
cat("incremental retention vs raw8: MPAI", round((r_mpai - r_base)/(r_raw - r_base), 3),
    " dims4", round((r_dims - r_base)/(r_raw - r_base), 3), "\n")

# ========== A2. 替代标量 composite ==========
d$eq8  <- rowMeans(d[z8])
d$eq4  <- rowMeans(d[DIMS])
pc <- prcomp(d[z8])
d$pc1 <- pc$x[, 1]
f_eq8 <- coxph(Surv(surv_time, died) ~ age + sex + eq8, data = d, cluster = communityID)
f_eq4 <- coxph(Surv(surv_time, died) ~ age + sex + eq4, data = d, cluster = communityID)
f_pc1 <- coxph(Surv(surv_time, died) ~ age + sex + pc1, data = d, cluster = communityID)
pc_hr <- summary(f_pc1)$coefficients["pc1", "coef"]
if (pc_hr < 0) { d$pc1 <- -d$pc1; f_pc1 <- coxph(Surv(surv_time, died) ~ age + sex + pc1, data = d, cluster = communityID) }
cat("=== A2 alternative scalars ===\n")
cat("equal-weight-8: dC", dC(f_eq8), " | equal-weight-4-dims: dC", dC(f_eq4),
    " | PCA-PC1: dC", dC(f_pc1), " (MPAI 0.0064, weighted8 0.0229, dims4 0.0185)\n")

# ========== C. 胱抑素C缺失电池 ==========
mk <- function(raw_id) {
  s <- gsub('"', '', gsub("\\s", "", as.character(raw_id)), fixed = TRUE)
  len <- nchar(s)
  if (all(len[!is.na(len)] == 12)) return(s)
  paste0(sprintf("%09d", suppressWarnings(as.numeric(substr(s, 1, pmax(len - 2, 0))))), "0",
         substr(s, len - 1, len))
}
bfile <- list.files(file.path(base, "chalrs2011"), pattern = "blood.*\\.dta$",
                    ignore.case = TRUE, full.names = TRUE)[1]
bl <- read_dta(bfile); bl$ID <- mk(bl$ID)
bl <- bl %>% transmute(ID, cyst_obs = as.numeric(cystatinc)) %>% filter(!is.na(cyst_obs))

dc <- analytic_df %>% left_join(bl, by = "ID") %>% mutate(cys_assayed = !is.na(cyst_obs))
cat("=== C1 complete-case (cystatin assayed) ===\n")
dcc <- dc %>% filter(cys_assayed)
cat("N =", nrow(dcc), " events =", sum(dcc$died), "\n")
f_cc <- coxph(Surv(surv_time, died) ~ age + sex + z_renal, data = dcc, cluster = communityID)
s_cc <- summary(f_cc)$coefficients
cat("complete-case z_renal HR:", round(exp(s_cc["z_renal","coef"]), 3),
    "(", round(exp(s_cc["z_renal","coef"] - 1.96*s_cc["z_renal","robust se"]), 3), "-",
    round(exp(s_cc["z_renal","coef"] + 1.96*s_cc["z_renal","robust se"]), 3), ")\n")

cat("=== C2 missingness model (predictors of cystatin missingness) ===\n")
gm <- glm(cys_assayed ~ age + sex + education + chronic_count + smoking + alcohol + died,
          data = dc, family = binomial)
sg <- summary(gm)$coefficients
print(data.frame(OR = round(exp(sg[rownames(sg) != "(Intercept)", "Estimate"]), 2),
                 p = signif(sg[rownames(sg) != "(Intercept)", "Pr(>|z|)"], 2)))

# ========== D. 维度-年龄相关 ==========
cat("=== D dimension-age correlations ===\n")
print(round(sapply(DIMS, function(v) cor(analytic_df[[v]], analytic_df$age, use = "complete.obs")), 3))
