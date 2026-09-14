# 10_table1_fill.R — Table 1 stratified by vital status + manuscript first-paragraph numbers
rm(list = ls()); library(dplyr); library(haven)
base <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT
cov <- read.csv(file.path(base, "processed/charls_covars.csv"), colClasses = c(ID = "character"))

cat("\n===== manuscript first-paragraph numbers =====\n")
cat("N =", nrow(cov), " deaths =", sum(cov$died), "\n")
cat("mean age:", round(mean(cov$age), 1), " women%:", round(100*mean(cov$sex == 2), 1),
    " FU median:", round(median(cov$surv_time), 1), "\n")
pres <- read.csv(file.path(base, "output/pooled_20imputations.csv"))
cat("MPAI pooled HR range:", paste(round(range(exp(pres$logHR_mpai)), 3), collapse = "-"), "\n")

d0 <- cov %>% filter(died == 0); d1 <- cov %>% filter(died == 1)
rcont <- function(v, label, digits = 1) {
  m0 <- mean(d0[[v]], na.rm = TRUE); s0 <- sd(d0[[v]], na.rm = TRUE)
  m1 <- mean(d1[[v]], na.rm = TRUE); s1 <- sd(d1[[v]], na.rm = TRUE)
  p <- suppressWarnings(tryCatch(t.test(d0[[v]], d1[[v]])$p.value, error = function(e) NA))
  data.frame(var = label,
             alive = sprintf(paste0("%.", digits, "f (%.", digits, "f)"), m0, s0),
             dead  = sprintf(paste0("%.", digits, "f (%.", digits, "f)"), m1, s1),
             p = ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))) }
rcat <- function(v, label, lev) {
  tb <- table(cov[[v]], cov$died); pr <- prop.table(tb, 2)
  rows <- lapply(seq_along(lev), function(i) data.frame(
    var = paste0("  ", names(lev)[i]),
    alive = sprintf("%d (%.1f%%)", tb[lev[i], "0"], 100 * pr[lev[i], "0"]),
    dead  = sprintf("%d (%.1f%%)", tb[lev[i], "1"], 100 * pr[lev[i], "1"]), p = ""))
  p <- suppressWarnings(tryCatch(chisq.test(tb)$p.value, error = function(e) NA))
  bind_rows(data.frame(var = label, alive = "", dead = "",
      p = ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))), bind_rows(rows)) }

cov$edu3 <- cut(cov$education, c(0, 3, 6, 99), labels = c("low", "mid", "high"))
tabA <- bind_rows(
  rcont("age", "Age, years"),
  rcat("sex", "Sex", c(Male = "1", Female = "2")),
  rcat("edu3", "Education", c("<= primary/sishu" = "low", "elementary-high" = "mid", "vocational+" = "high")),
  rcat("smoking", "Smoking", c(Yes = "1", No = "0")),
  rcat("alcohol", "Alcohol use", c(Yes = "1", No = "0")),
  rcont("chronic_count", "Chronic conditions, count"),
  rcont("self_health", "Self-rated health, 1-5"),
  rcont("cesd_score", "CES-D-10"),
  rcont("tics_score", "Word recall (see 11_tics_check for corrected scoring)"),
  rcont("func_score", "ADL/IADL, 15-60"),
  rcont("surv_time", "Follow-up, years"))

# blood biomarkers (raw observed values, blood subsample)
mk <- function(raw_id) {
  s <- gsub('"', '', gsub("\\s", "", as.character(raw_id)), fixed = TRUE)
  len <- nchar(s)
  if (all(len[!is.na(len)] == 12)) return(s)
  paste0(sprintf("%09d", suppressWarnings(as.numeric(substr(s, 1, pmax(len - 2, 0))))), "0",
         substr(s, len - 1, len)) }
bfile <- list.files(file.path(base, "chalrs2011"), pattern = "blood.*\\.dta$",
                    ignore.case = TRUE, full.names = TRUE)[1]
bl <- read_dta(bfile); bl$ID <- mk(bl$ID)
bl <- bl %>% transmute(ID, HbA1c = as.numeric(newhba1c), CRP = as.numeric(newcrp),
    cystatin = as.numeric(cystatinc), cholesterol = as.numeric(newcho)) %>%
  inner_join(cov %>% select(ID, died), by = "ID")
rbio_mean <- function(v, label) {
  x0 <- bl[[v]][bl$died == 0]; x1 <- bl[[v]][bl$died == 1]
  p <- suppressWarnings(tryCatch(t.test(x0, x1)$p.value, error = function(e) NA))
  data.frame(var = paste0(label, ", mean (SD) [n]"),
    alive = sprintf("%.2f (%.2f) [%d]", mean(x0, na.rm = T), sd(x0, na.rm = T), sum(!is.na(x0))),
    dead  = sprintf("%.2f (%.2f) [%d]", mean(x1, na.rm = T), sd(x1, na.rm = T), sum(!is.na(x1))),
    p = ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))) }
rbio_med <- function(v, label) {
  x0 <- bl[[v]][bl$died == 0]; x1 <- bl[[v]][bl$died == 1]
  p <- suppressWarnings(tryCatch(wilcox.test(x0, x1)$p.value, error = function(e) NA))
  fmt <- function(x) sprintf("%.2f (%.2f-%.2f) [%d]", median(x, na.rm = T),
    quantile(x, .25, na.rm = T), quantile(x, .75, na.rm = T), sum(!is.na(x)))
  data.frame(var = paste0(label, ", median (IQR) [n]"), alive = fmt(x0), dead = fmt(x1),
    p = ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))) }
tabB <- bind_rows(rbio_mean("HbA1c", "HbA1c, %"), rbio_med("CRP", "CRP, mg/L"),
                  rbio_mean("cystatin", "Cystatin C, mg/L"), rbio_mean("cholesterol", "Total cholesterol, mg/dL"))
cat("\n===== Table 1. By vital status =====\n")
print(bind_rows(tabA, tabB), right = FALSE, row.names = FALSE)
