# 16_test_set_eval.R —— test-set-only ΔC（Cox 在 test 集内重拟合）
rm(list = ls()); library(dplyr); library(survival)
base <- "/Users/a1234/Desktop/桌面/chalrs"
load(file.path(base, "output/metrics_final.RData"))
test_ids <- read.csv(file.path(base, "processed/test_ids.csv"), colClasses = c(ID = "character"))
dt <- analytic_df %>% semi_join(test_ids, by = "ID")
cat("test N =", nrow(dt), " events =", sum(dt$died), "\n")

f_b <- coxph(Surv(surv_time, died) ~ age + sex, data = dt, cluster = communityID)
f_m <- coxph(Surv(surv_time, died) ~ age + sex + MPAI, data = dt, cluster = communityID)
f_d <- coxph(Surv(surv_time, died) ~ age + sex + z_metabolic + z_central + z_funcinfl + z_renal,
             data = dt, cluster = communityID)
cb <- as.numeric(summary(f_b)$concordance["C"])
cat(sprintf("test-set: base C=%.4f | MPAI dC=%.4f | dims4 dC=%.4f\n",
            cb,
            as.numeric(summary(f_m)$concordance["C"]) - cb,
            as.numeric(summary(f_d)$concordance["C"]) - cb))
cat("(full cohort reference: MPAI 0.0064, dims4 0.0185)\n")
