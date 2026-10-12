# ============================================================
# 20_crp_acute_phase_sensitivity.R
# CRP acute-phase sensitivity: does restricting CRP to <=10 mg/L change the
# eight-indicator model dC-index, and the dimension hierarchy?  (appendix S6)
# Spec: identical to the main-text cox_raw8 model (imputation-1 standardised
# eight indicators, same analytic sample).
# ============================================================
rm(list = ls()); options(stringsAsFactors = FALSE)
suppressMessages({library(dplyr); library(haven); library(survival)})

base_path <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT: path to local CHARLS root
PROC_DIR <- file.path(base_path, "processed")
OUT_DIR  <- file.path(base_path, "output")

load(file.path(OUT_DIR, "metrics_final.RData"))   # analytic_df, DIMS, flips
X8   <- read.csv(file.path(PROC_DIR, "charls_std_for_vae.csv"))
ids8 <- read.csv(file.path(PROC_DIR, "charls_ids.csv"), colClasses = c(ID = "character"))
X8$ID <- ids8$ID

z8 <- c("tics_cog","self_health_rev","cesd_rev","physio_proxy",
        "hba1c_rev","crp_rev","cystatin_c_rev","total_chol_rev")

d <- analytic_df %>% left_join(X8, by = "ID")
stopifnot(nrow(d) == nrow(analytic_df), !any(is.na(d$crp_rev)))

extract_c <- function(m) as.numeric(summary(m)$concordance["C"])

# ---------- 1. Read raw CRP and recover the raw scale ----------
mk <- function(raw_id) {
  s <- gsub('"', '', gsub("\\s", "", as.character(raw_id)), fixed = TRUE)
  len <- nchar(s)
  if (all(len[!is.na(len)] == 12)) return(s)
  paste0(sprintf("%09d", suppressWarnings(as.numeric(substr(s, 1, pmax(len - 2, 0))))),
         "0", substr(s, len - 1, len))
}
bl <- read_dta(file.path(base_path, "chalrs2011", "Blood_20140429.dta"))
bl$ID <- mk(bl$ID)
bl <- bl %>% transmute(ID, crp_raw_obs = suppressWarnings(as.numeric(newcrp))) %>%
  filter(!is.na(crp_raw_obs))
d <- d %>% left_join(bl, by = "ID")

cat("raw CRP observed:", sum(!is.na(d$crp_raw_obs)), "/", nrow(d), "\n")

# Recover the standardisation map: crp_rev = -(log(crp+1) - mu)/sd
#   => log(crp+1) = -crp_rev*sd + mu   (exact for observed cases)
ok <- !is.na(d$crp_raw_obs)
lm_map <- lm(log(d$crp_raw_obs[ok] + 1) ~ d$crp_rev[ok])
sd_c <- unname(-coef(lm_map)["d$crp_rev[ok]"]); mu_c <- unname(coef(lm_map)["(Intercept)"])
cat(sprintf("recovered map: sd=%.4f mu=%.4f  R2=%.6f\n", sd_c, mu_c, summary(lm_map)$r.squared))
d$crp_raw <- exp(-d$crp_rev * sd_c + mu_c) - 1

cat(sprintf("check (observed): max abs err %.3e mg/L, median err %.3e\n",
            max(abs(d$crp_raw[ok] - d$crp_raw_obs[ok]), na.rm = TRUE),
            median(abs(d$crp_raw[ok] - d$crp_raw_obs[ok]), na.rm = TRUE)))

cat(sprintf("\nraw CRP > 10 mg/L: observed %d, full sample %d (%.1f%%)\n",
            sum(d$crp_raw_obs > 10, na.rm = TRUE),
            sum(d$crp_raw > 10, na.rm = TRUE),
            100 * mean(d$crp_raw > 10, na.rm = TRUE)))

# ---------- 2. Baseline (main-text specification) ----------
f_raw8 <- as.formula(paste("Surv(surv_time, died) ~ age + sex +", paste(z8, collapse = " + ")))
c_base <- extract_c(coxph(Surv(surv_time, died) ~ age + sex, data = d, cluster = communityID))
f_orig <- coxph(f_raw8, data = d, cluster = communityID)
dC_orig <- extract_c(f_orig) - c_base
cat(sprintf("\n[main spec] base C=%.4f  raw8 C=%.4f  dC=%.4f\n",
            c_base, extract_c(f_orig), dC_orig))

# ---------- 3. Variant A: cap CRP at 10 mg/L (sample unchanged) ----------
d$crp_rev_capped <- ifelse(d$crp_raw > 10, -(log(11) - mu_c) / sd_c, d$crp_rev)
z8_cap <- z8; z8_cap[z8_cap == "crp_rev"] <- "crp_rev_capped"
f_cap <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +", paste(z8_cap, collapse = " + "))),
               data = d, cluster = communityID)
dC_cap <- extract_c(f_cap) - c_base
cat(sprintf("[A cap@10]  raw8 C=%.4f  dC=%.4f  (delta dC=%+.4f)\n",
            extract_c(f_cap), dC_cap, dC_cap - dC_orig))

# ---------- 4. Variant B: exclude participants with CRP > 10 mg/L ----------
d_ex <- d %>% filter(is.na(crp_raw) | crp_raw <= 10)
c_base_ex <- extract_c(coxph(Surv(surv_time, died) ~ age + sex, data = d_ex, cluster = communityID))
f_ex <- coxph(f_raw8, data = d_ex, cluster = communityID)
dC_ex <- extract_c(f_ex) - c_base_ex
cat(sprintf("[B exclude>10] n=%d events=%d  base C=%.4f  raw8 C=%.4f  dC=%.4f  (delta dC=%+.4f)\n",
            nrow(d_ex), sum(d_ex$died), c_base_ex, extract_c(f_ex), dC_ex, dC_ex - dC_orig))

# ---------- 5. Change in the CRP coefficient ----------
s_o <- summary(f_orig)$coefficients["crp_rev", ]
s_c <- summary(f_cap)$coefficients["crp_rev_capped", ]
cat(sprintf("\nCRP coef: original beta=%.4f (robust SE %.4f) | capped beta=%.4f (%.4f)\n",
            s_o["coef"], s_o["robust se"], s_c["coef"], s_c["robust se"]))

# ---------- 6. Is the dimension hierarchy preserved? ----------
dim_of <- function(dat, label) {
  fb <- extract_c(coxph(Surv(surv_time, died) ~ age + sex, data = dat, cluster = communityID))
  fd <- coxph(as.formula(paste("Surv(surv_time, died) ~ age + sex +",
                               paste(DIMS, collapse = " + "))), data = dat, cluster = communityID)
  sd_ <- summary(fd)$coefficients
  hr <- data.frame(dimension = DIMS, HR = round(exp(sd_[DIMS, "coef"]), 3),
                   lo = round(exp(sd_[DIMS,"coef"] - 1.96*sd_[DIMS,"robust se"]), 3),
                   hi = round(exp(sd_[DIMS,"coef"] + 1.96*sd_[DIMS,"robust se"]), 3),
                   scenario = label, row.names = NULL)
  hr$dC_dims <- round(extract_c(fd) - fb, 4)
  cat(sprintf("\n[%s] four-dim dC=%.4f  dimension HRs (descending):\n", label, hr$dC_dims[1]))
  print(hr[order(-hr$HR), c("dimension","HR","lo","hi")])
  hr
}
h_orig <- dim_of(d, "original")
h_ex   <- dim_of(d_ex, "crp_le_10")

cat("\ndimension ordering identical: original=", paste(h_orig$dimension[order(-h_orig$HR)], collapse=" > "),
    " | crp<=10=", paste(h_ex$dimension[order(-h_ex$HR)], collapse=" > "), "\n")

# ---------- 7. Output ----------
out <- data.frame(
  scenario = c("original", "crp_capped_at_10", "crp_le_10_excluded"),
  n        = c(nrow(d), nrow(d), nrow(d_ex)),
  events   = c(sum(d$died), sum(d$died), sum(d_ex$died)),
  base_C   = round(c(c_base, c_base, c_base_ex), 4),
  raw8_C   = round(c(extract_c(f_orig), extract_c(f_cap), extract_c(f_ex)), 4),
  raw8_dC  = round(c(dC_orig, dC_cap, dC_ex), 4))
write.csv(out, file.path(OUT_DIR, "crp_acute_phase_sensitivity.csv"), row.names = FALSE)
write.csv(bind_rows(h_orig, h_ex), file.path(OUT_DIR, "crp_sensitivity_dims.csv"), row.names = FALSE)
cat("\nwrote: output/crp_acute_phase_sensitivity.csv, output/crp_sensitivity_dims.csv\n")
