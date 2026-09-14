# 15_dim_sensitivity 的 R 端计算（python 跑完后执行）
for (k in c(1,2,3,5,6)) {
  zk <- read.csv(sprintf("%s/processed/vae_dim%d.csv", base, k), colClasses = c(ID = "character"))
  dk <- analytic_df %>% left_join(zk, by = "ID")
  zc <- grep("^z\\d+$", names(dk), value = TRUE)
  for (v in zc) {   # 符号对齐
    f1 <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +", v)), data = dk)
    if (coef(f1)[3] < 0) dk[[v]] <- -dk[[v]]
  }
  dk$MPAIk <- sqrt(rowSums(dk[zc]^2))
  f_dim <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +", paste(zc, collapse="+"))), data = dk, cluster = communityID)
  f_sca <- coxph(Surv(surv_time, died) ~ age + sex + MPAIk, data = dk, cluster = communityID)
  cat(sprintf("k=%d  dims dC=%.4f  scalar dC=%.4f\n", k,
      as.numeric(summary(f_dim)$concordance["C"]) - c_base,
      as.numeric(summary(f_sca)$concordance["C"]) - c_base))
}