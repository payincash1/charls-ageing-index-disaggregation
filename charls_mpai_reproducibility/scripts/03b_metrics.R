# 03b_metrics.R — main survival analysis (all patches applied)
rm(list = ls()); options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))
pkgs <- c("dplyr","survival","mice","survIDINRI","rms","riskRegression","dcurves","boot","ggplot2")
for (p in pkgs) { if (!require(p, character.only = TRUE, quietly = TRUE)) {
  install.packages(p); library(p, character.only = TRUE) } }

base <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT
PROC_DIR <- file.path(base, "processed"); OUT_DIR <- file.path(base, "output")

covars <- read.csv(file.path(PROC_DIR, "charls_covars.csv"), colClasses = c(ID = "character"))
vae    <- read.csv(file.path(PROC_DIR, "vae_scores.csv"), colClasses = c(ID = "character"))
stopifnot(all(c("z1","z2","z3","z4","MPAI") %in% names(vae)))

# dimension renaming per m01 loadings + sign alignment
DIM_MAP <- c(z_metabolic = "z1", z_central = "z2", z_funcinfl = "z3", z_renal = "z4")
DIMS <- names(DIM_MAP)
df <- covars %>% left_join(vae, by = "ID") %>% rename(!!!DIM_MAP)
stopifnot(all(DIMS %in% names(df)), !any(is.na(df$MPAI)))

flips <- sapply(DIMS, function(v) {
  f <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +", v)), data = df)
  ifelse(coef(f)[3] < 0, -1, 1) })
df[DIMS] <- sweep(df[DIMS], 2, flips, `*`)
cat("sign-flipped:", paste(names(flips[flips == -1])), "\n")

analytic_vars <- c("age","sex","surv_time","died","MPAI",DIMS,
                   "communityID","education","smoking","alcohol","chronic_count")
analytic_df <- df %>% filter(complete.cases(across(all_of(analytic_vars))))
analytic_df <- analytic_df %>% mutate(edu_cat3 = case_when(
  education <= 3 ~ "low", education <= 6 ~ "mid", TRUE ~ "high"))
cat(sprintf("Analytic N = %d, Events = %d\n", nrow(analytic_df), sum(analytic_df$died)))
extract_c <- function(m) as.numeric(summary(m)$concordance["C"])
f_dims <- paste("Surv(surv_time, died) ~ age + sex +", paste(DIMS, collapse = " + "))

cox_base <- coxph(Surv(surv_time, died) ~ age + sex, data = analytic_df, cluster = communityID)
cox_mpai <- coxph(Surv(surv_time, died) ~ age + sex + MPAI, data = analytic_df, cluster = communityID)
cox_dims <- coxph(as.formula(f_dims), data = analytic_df, cluster = communityID)
cox_full <- coxph(Surv(surv_time, died) ~ age + sex + edu_cat3 + smoking + alcohol +
                    chronic_count + MPAI, data = analytic_df, cluster = communityID)

# control models (raw 8 + Cox-weighted)
X8 <- read.csv(file.path(PROC_DIR, "charls_std_for_vae.csv"))
ids8 <- read.csv(file.path(PROC_DIR, "charls_ids.csv"), colClasses = c(ID = "character"))
X8$ID <- ids8$ID
z8 <- c("tics_cog","self_health_rev","cesd_rev","physio_proxy","hba1c_rev","crp_rev",
        "cystatin_c_rev","total_chol_rev")
adf <- analytic_df %>% left_join(X8, by = "ID")
cox_raw8 <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +",
                                   paste(z8, collapse = " + "))), data = adf, cluster = communityID)
f8 <- coxph(as.formula(paste("Surv(surv_time, died) ~", paste(z8, collapse = " + "))), data = adf)
adf$lp_w8 <- predict(f8)
cox_w8 <- coxph(Surv(surv_time, died) ~ age + sex + lp_w8, data = adf, cluster = communityID)

cat("\n=== Model Comparison ===\n")
mods <- list(base = cox_base, MPAI = cox_mpai, dims4 = cox_dims, full = cox_full,
             raw8 = cox_raw8, weighted8 = cox_w8)
cmp <- data.frame(model = names(mods), C = sapply(mods, extract_c),
                  dC = sapply(mods, extract_c) - extract_c(cox_base),
                  n = sapply(mods, function(m) m$n))
print(cmp); write.csv(cmp, file.path(OUT_DIR, "model_comparison.csv"), row.names = FALSE)
lrt <- 2 * (cox_mpai$loglik[2] - cox_base$loglik[2])
cat(sprintf("LRT MPAI vs base: chi2 = %.2f, p = %.4g\n", lrt, pchisq(lrt, 1, lower.tail = FALSE)))
cat("\n=== MPAI (robust SE) ===\n"); print(summary(cox_mpai)$coefficients)
cat("\n=== 4 dimensions (robust SE) ===\n"); print(summary(cox_dims)$coefficients)

# cluster bootstrap dC, B=1000
cat("\n=== Cluster bootstrap dC (B=1000) ===\n")
set.seed(2024); comms <- unique(analytic_df$communityID); dC <- numeric(1000)
for (b in 1:1000) {
  idx <- unlist(lapply(sample(comms, replace = TRUE),
                       function(cc) which(analytic_df$communityID == cc)))
  d <- analytic_df[idx, ]
  dC[b] <- tryCatch({
    cb <- coxph(Surv(surv_time, died) ~ age + sex, data = d)
    cm <- coxph(Surv(surv_time, died) ~ age + sex + MPAI, data = d)
    extract_c(cm) - extract_c(cb) }, error = function(e) NA_real_)
  if (b %% 200 == 0) cat("  boot", b, "\n") }
dC_ci <- quantile(dC, c(0.025, 0.5, 0.975), na.rm = TRUE)
cat(sprintf("dC = %.4f (95%% CI %.4f-%.4f)\n", dC_ci[2], dC_ci[1], dC_ci[3]))
write.csv(data.frame(dC = dC), file.path(OUT_DIR, "bootstrap_dC.csv"), row.names = FALSE)

# pooled across 20 aligned imputations
cat("\n=== Pooled across 20 aligned imputations ===\n")
pool_one <- function(i) {
  vm <- read.csv(file.path(PROC_DIR, sprintf("vae_scores_aligned_m%02d.csv", i)),
                 colClasses = c(ID = "character")) %>% rename(!!!DIM_MAP)
  vm[DIMS] <- sweep(vm[DIMS], 2, flips, `*`)
  vm$MPAI <- sqrt(rowSums(vm[DIMS]^2))
  di <- covars %>% left_join(vm, by = "ID") %>%
    filter(complete.cases(across(all_of(analytic_vars))))
  fm <- coxph(Surv(surv_time, died) ~ age + sex + MPAI, data = di)
  fd <- coxph(as.formula(f_dims), data = di)
  sm <- summary(fm)$coefficients; sd_ <- summary(fd)$coefficients
  out <- data.frame(m = i, logHR_mpai = sm["MPAI","coef"], se_mpai = sm["MPAI","se(coef)"])
  for (v in DIMS) { out[[paste0("logHR_", v)]] <- sd_[v,"coef"]
                    out[[paste0("se_", v)]]    <- sd_[v,"se(coef)"] }
  out }
pres <- bind_rows(lapply(1:20, pool_one))
write.csv(pres, file.path(OUT_DIR, "pooled_20imputations.csv"), row.names = FALSE)
pool_line <- function(est, se, name) {
  p <- pool.scalar(Q = est, U = se^2)
  ci <- exp(p$qbar + c(-1,1) * qt(0.975, df = p$df) * sqrt(p$t))
  cat(sprintf("%s pooled HR = %.3f (95%% CI %.3f-%.3f)\n", name, exp(p$qbar), ci[1], ci[2]))
  data.frame(term = name, HR = exp(p$qbar), lo = ci[1], hi = ci[2]) }
pooled_res <- pool_line(pres$logHR_mpai, pres$se_mpai, "MPAI")
for (v in DIMS) pooled_res <- bind_rows(pooled_res,
  pool_line(pres[[paste0("logHR_", v)]], pres[[paste0("se_", v)]], v))
write.csv(pooled_res, file.path(OUT_DIR, "pooled_HRs.csv"), row.names = FALSE)

# NRI/IDI (secondary; errors surfaced)
cat("\n=== NRI / IDI at 5 years ===\n")
lp_base <- predict(cox_base); lp_mpai <- predict(cox_mpai); lp_dims <- predict(cox_dims)
tt <- analytic_df$surv_time; ss <- analytic_df$died
nri_idi <- function(lp_b, lp_n) {
  n <- length(tt)
  get_risk <- function(lp) { bh <- basehaz(coxph(Surv(tt, ss) ~ lp)); 
    1 - exp(-approx(bh$time, bh$hazard, 5, method="constant", rule=2, f=0)$y * exp(lp)) }
  rb <- get_risk(lp_b); rn <- get_risk(lp_n)
  cens <- survfit(Surv(tt, ss == 0) ~ 1)
  Ga <- approx(cens$time, cens$surv, xout = pmin(tt, max(cens$time)),
               method = "constant", rule = 2, f = 0)$y
  G5 <- approx(cens$time, cens$surv, 5, method = "constant", rule = 2, f = 0)$y
  Ga[Ga < 0.05] <- 0.05
  w_case <- ifelse(ss == 1 & tt <= 5, 1/Ga, 0); w_ctl <- ifelse(tt > 5, 1/G5, 0)
  est <- function(perturb) { w <- if (perturb) rexp(n) else rep(1, n)
    wc <- w*w_case; wk <- w*w_ctl
    nri <- (sum(wc*(rn>rb)) - sum(wc*(rn<rb)))/sum(wc) + (sum(wk*(rn<rb)) - sum(wk*(rn>rb)))/sum(wk)
    idi <- sum(wc*(rn-rb))/sum(wc) - sum(wk*(rn-rb))/sum(wk); c(nri, idi) }
  set.seed(2024); point <- est(FALSE); reps <- replicate(200, est(TRUE))
  ci <- apply(reps, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
  data.frame(measure = c("cNRI","IDI"), estimate = round(point, 4),
             lo = round(ci[1,], 4), hi = round(ci[2,], 4)) }
r_mpai <- nri_idi(lp_base, lp_mpai); r_dims <- nri_idi(lp_base, lp_dims)
print(r_mpai); print(r_dims)
write.csv(rbind(cbind(comparison = "MPAI_vs_base", r_mpai),
                cbind(comparison = "dims_vs_base", r_dims)),
          file.path(OUT_DIR, "nri_idi_manual.csv"), row.names = FALSE)

# calibration
cat("\n=== Calibration ===\n")
lp <- predict(cox_mpai)
slope_fit <- coxph(Surv(surv_time, died) ~ lp, data = analytic_df)
cat(sprintf("Calibration slope = %.3f (95%% CI %.3f-%.3f)\n", coef(slope_fit)["lp"],
            confint(slope_fit)["lp",1], confint(slope_fit)["lp",2]))
suppressMessages({ dd <- datadist(analytic_df); options(datadist = "dd") })
cph_m <- cph(Surv(surv_time, died) ~ age + sex + MPAI, data = analytic_df,
             x = TRUE, y = TRUE, surv = TRUE)
cal <- calibrate(cph_m, u = 5, B = 200)
png(file.path(OUT_DIR, "calibration_plot_5y.png"), width = 2000, height = 1600, res = 200)
plot(cal); dev.off()

# DCA
rs <- tryCatch(riskRegression::predictRisk(cox_mpai, analytic_df, times = 5), error = function(e) NULL)
rb <- tryCatch(riskRegression::predictRisk(cox_base, analytic_df, times = 5), error = function(e) NULL)
if (!is.null(rs) && !is.null(rb)) {
  analytic_df$risk_mpai <- as.numeric(rs); analytic_df$risk_base <- as.numeric(rb)
  dc <- dca(Surv(surv_time, died) ~ risk_base + risk_mpai, data = analytic_df,
            time = 5, thresholds = seq(0.05, 0.40, by = 0.01))
  ggsave(file.path(OUT_DIR, "dca_5y.png"), plot(dc, smooth = TRUE),
         width = 8, height = 6, dpi = 300) }

# survey-weighted
cat("\n=== Survey-weighted Cox ===\n")
dfw <- analytic_df %>% filter(!is.na(weight_blood))
cox_w <- coxph(Surv(surv_time, died) ~ age + sex + MPAI, data = dfw,
               weights = weight_blood, cluster = communityID)
print(summary(cox_w)$coefficients["MPAI", ])

# interactions
cat("\n=== Interactions ===\n")
int_terms <- c(paste0("age:", DIMS), paste0("sex:", DIMS))
cox_int <- coxph(as.formula(paste(f_dims, "+", paste(int_terms, collapse = " + "))),
                 data = analytic_df, cluster = communityID)
si <- summary(cox_int)$coefficients
print(si[grepl(":", rownames(si)), c("coef", "robust se", "Pr(>|z|)")])

# P0-4 sensitivity battery
cat("\n=== Sensitivity battery ===\n")
fit3 <- function(d, label) {
  f <- coxph(as.formula(f_dims), data = d, cluster = communityID)
  s <- summary(f)$coefficients
  data.frame(scenario = label, term = rownames(s), HR = exp(s[,"coef"]),
             lo = exp(s[,"coef"] - 1.96*s[,"robust se"]),
             hi = exp(s[,"coef"] + 1.96*s[,"robust se"]), p = s[,"Pr(>|z|)"], row.names = NULL) }
bat <- bind_rows(fit3(analytic_df, "overall"),
                 fit3(filter(analytic_df, age >= 60), "age60plus"),
                 fit3(filter(analytic_df, age >= 65), "age65plus"),
                 fit3(filter(analytic_df, died == 0 | surv_time > 2), "exclude_first2y_deaths"))
ad2 <- analytic_df %>% mutate(t2 = pmin(surv_time, 3.5),
                              e2 = as.integer(died == 1 & surv_time <= 3.5))
cox_p1 <- coxph(as.formula(paste("Surv(t2, e2) ~ age + sex +", paste(DIMS, collapse = " + "))),
                data = ad2, cluster = communityID)
late <- analytic_df %>% filter(surv_time > 3.5) %>% mutate(t3 = surv_time - 3.5)
cox_p2 <- coxph(as.formula(paste("Surv(t3, died) ~ age + sex +", paste(DIMS, collapse = " + "))),
                data = late, cluster = communityID)
mk <- function(s, se_col, label) data.frame(scenario = label, term = rownames(s),
  HR = exp(s[,"coef"]), lo = exp(s[,"coef"] - 1.96*s[, se_col]),
  hi = exp(s[,"coef"] + 1.96*s[, se_col]), p = s[,"Pr(>|z|)"], row.names = NULL)
bat <- bind_rows(bat, mk(summary(cox_p1)$coefficients, "se(coef)", "period_2011_2014"),
                      mk(summary(cox_p2)$coefficients, "se(coef)", "period_2014_2018"))
dd_mid <- rep(as.Date(NA), nrow(analytic_df))
ok_mid <- analytic_df$died == 1 & !is.na(analytic_df$death_year) &
          analytic_df$death_year >= 1950 & analytic_df$death_year <= 2020
dd_mid[ok_mid] <- as.Date(sprintf("%04d-07-01", as.integer(analytic_df$death_year[ok_mid])))
analytic_df$dd_mid <- dd_mid
dm <- analytic_df %>%
  mutate(last2 = if_else(died == 1 & !is.na(dd_mid), dd_mid, as.Date("2018-12-31")),
         st2 = as.numeric(difftime(last2, as.Date("2011-07-01"), units = "days"))/365.25) %>%
  filter(st2 > 0)
cox_mid <- coxph(as.formula(paste("Surv(st2, died) ~ age + sex +", paste(DIMS, collapse = " + "))),
                 data = dm, cluster = communityID)
bat <- bind_rows(bat, mk(summary(cox_mid)$coefficients, "robust se", "death_midpoint"))
write.csv(bat, file.path(OUT_DIR, "sensitivity_battery.csv"), row.names = FALSE)
print(bat[bat$term %in% DIMS, ])

tab_adl <- analytic_df %>% mutate(func_q = ntile(func_score, 4)) %>% group_by(func_q) %>%
  summarise(n = n(), events = sum(died), event_rate = round(mean(died), 3),
            func_mean = round(mean(func_score), 1), .groups = "drop")
write.csv(tab_adl, file.path(OUT_DIR, "adl_distribution.csv"), row.names = FALSE)

save(cox_base, cox_mpai, cox_dims, cox_full, cox_raw8, cox_w8, analytic_df,
     pres, pooled_res, bat, dC, flips, DIMS,
     file = file.path(OUT_DIR, "metrics_final.RData"))
cat("\n03b done\n")
