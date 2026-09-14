# 04_calibration.R — 5-year AUC/Brier, calibration-in-the-large, decile calibration plot
rm(list = ls()); options(stringsAsFactors = FALSE)
library(survival); library(dplyr); library(ggplot2); library(riskRegression)
base <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT
OUT_DIR <- file.path(base, "output")
load(file.path(OUT_DIR, "metrics_final.RData"))
d <- analytic_df

f_base <- coxph(Surv(surv_time, died) ~ age + sex, data = d, x = TRUE)
f_mpai <- coxph(Surv(surv_time, died) ~ age + sex + MPAI, data = d, x = TRUE)
f_dims <- coxph(Surv(surv_time, died) ~ age + sex + z_metabolic + z_central +
                  z_funcinfl + z_renal, data = d, x = TRUE)
d$risk_base <- as.numeric(predictRisk(f_base, d, times = 5))
d$risk_mpai <- as.numeric(predictRisk(f_mpai, d, times = 5))
d$risk_dims <- as.numeric(predictRisk(f_dims, d, times = 5))
stopifnot(!anyNA(d$risk_dims))

km5 <- function(dd) 1 - summary(survfit(Surv(surv_time, died) ~ 1, data = dd), times = 5)$surv[1]
citl <- data.frame(model = c("base","MPAI","dims4"),
  mean_pred_5y = round(c(mean(d$risk_base), mean(d$risk_mpai), mean(d$risk_dims)), 4),
  obs_KM_5y = round(km5(d), 4))
print(citl)
write.csv(citl, file.path(OUT_DIR, "calibration_in_large.csv"), row.names = FALSE)

dec_cal <- function(risk) {
  dd <- d %>% mutate(g = ntile(!!sym(risk), 10))
  dd %>% group_by(g) %>% group_modify(~{
    sf <- survfit(Surv(surv_time, died) ~ 1, data = .x)
    data.frame(pred = mean(.x[[risk]]), obs = 1 - summary(sf, times = 5)$surv[1],
               n = nrow(.x)) }) %>% ungroup() %>% mutate(model = risk) }
cal_all <- bind_rows(dec_cal("risk_base"), dec_cal("risk_mpai"), dec_cal("risk_dims"))
write.csv(cal_all, file.path(OUT_DIR, "calibration_deciles.csv"), row.names = FALSE)
lim <- max(cal_all$pred, cal_all$obs, na.rm = TRUE) * 1.1
p <- ggplot(cal_all, aes(x = pred, y = obs, color = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_line() + geom_point(aes(size = n), alpha = 0.85) +
  scale_color_manual(values = c(risk_base="grey55", risk_mpai="#D7301F", risk_dims="#2171B5"),
                     labels = c("Age + sex","MPAI","4 dimensions")) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(x = "Predicted 5-year risk", y = "Observed 5-year risk (Kaplan-Meier)",
       color = NULL, size = "N") + theme_bw(base_size = 13)
ggsave(file.path(OUT_DIR, "calibration_deciles_5y.png"), p, width = 7, height = 6, dpi = 300)

sc <- tryCatch(Score(list(base = f_base, MPAI = f_mpai, dims4 = f_dims),
    formula = Surv(surv_time, died) ~ 1, data = d, times = 5,
    metrics = c("auc","brier"), B = 200),
  error = function(e) { message("Score error: ", conditionMessage(e)); NULL })
if (!is.null(sc)) {
  print(sc$AUC); print(sc$Brier)
  capture.output(print(sc$AUC), print(sc$Brier),
                 file = file.path(OUT_DIR, "auc_brier_5y.txt")) }
cat("04 done\n")
