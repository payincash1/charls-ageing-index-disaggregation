# 01b_impute_m20.R — CHARLS data preparation, death linkage, m=20 imputation
rm(list = ls()); options(stringsAsFactors = FALSE)
pkgs <- c("dplyr","haven","mice","survival","readr")
for (p in pkgs) { if (!require(p, character.only = TRUE, quietly = TRUE)) {
  install.packages(p); library(p, character.only = TRUE) } }

base <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT: local CHARLS root
years <- c(2011, 2013, 2015, 2018)
PROC_DIR <- file.path(base, "processed"); OUT_DIR <- file.path(base, "output")
dir.create(PROC_DIR, showWarnings = FALSE); dir.create(OUT_DIR, showWarnings = FALSE)

# ---- ID conversion: raw 10/11-digit numeric ID -> 12-char standard ----
make_std_id <- function(raw_id) {
  s <- gsub('"', '', gsub("\\s", "", as.character(raw_id)), fixed = TRUE)
  len <- nchar(s)
  if (all(len[!is.na(len)] == 12)) {
    if (!all(substr(s, 10, 10) == "0", na.rm = TRUE)) stop("12-char ID, pos 10 not 0")
    return(s) }
  if (!all(len %in% c(9,10,11) | is.na(len))) stop("Unhandled ID length")
  pn <- substr(s, len-1, len); hh <- substr(s, 1, pmax(len-2, 0))
  out <- paste0(sprintf("%09d", suppressWarnings(as.numeric(hh))), "0", pn)
  out[len == 9 | is.na(len)] <- NA; out
}
convert_ids <- function(df, fname) {
  if (!"ID" %in% names(df)) stop(fname, " no ID")
  df$ID <- make_std_id(df$ID)
  if (!all(nchar(df$ID) == 12 | is.na(df$ID))) stop(fname, " non-12char ID")
  if (anyDuplicated(df$ID[!is.na(df$ID)])) stop(fname, " duplicate ID")
  df
}
all_dirs <- list.dirs(base, recursive = FALSE); dir_names <- basename(all_dirs)
match_year_dir <- function(yr) { idx <- grep(paste0(".*", yr, ".*"), dir_names, ignore.case = TRUE)
  if (length(idx) == 0) return(NA); all_dirs[idx[which.min(nchar(dir_names[idx]))]] }
year_dir <- sapply(years, match_year_dir); names(year_dir) <- paste0("w", seq_along(years))
find_dta <- function(dp, pat) {
  fs <- list.files(dp, pattern = paste0(pat, ".*\\.dta$"), ignore.case = TRUE, full.names = TRUE)
  if (length(fs) == 0) { ad <- list.files(dp, "\\.dta$", ignore.case = TRUE, full.names = TRUE)
    h <- grep(pat, basename(ad), ignore.case = TRUE, value = TRUE)
    if (length(h) == 0) stop("not found: ", pat); return(ad[basename(ad) == h[1]]) }
  fs[1]
}

cat("Reading Wave 1 (2011)...\n")
demog   <- convert_ids(read_dta(find_dta(year_dir["w1"], "demographic_background")), "demog")
health  <- convert_ids(read_dta(find_dta(year_dir["w1"], "health_status_and_functioning")), "health")
blood   <- convert_ids(read_dta(find_dta(year_dir["w1"], "blood")), "blood")
weights <- convert_ids(read_dta(find_dta(year_dir["w1"], "weight")), "weights")

# alcohol variable auto-detect (label contains "drink")
all_lbl <- sapply(names(health), function(v) { x <- attr(health[[v]], "label")
  ifelse(is.null(x), "", tolower(as.character(x))) })
alc_cand <- names(health)[grepl("drink", all_lbl)]
alc_var <- if (length(alc_cand) > 0)
  alc_cand[which.min(sapply(alc_cand, function(v) mean(is.na(health[[v]]))))] else NULL
cat("alcohol variable:", alc_var, "\n")

find_tics_vars <- function(px, df) {
  es <- intersect(paste0(px, "s", 1:10), names(df)); eu <- intersect(paste0(px, "_", 1:10), names(df))
  if (length(es) >= length(eu)) es else eu }
imm <- find_tics_vars("dc006", health); del <- find_tics_vars("dc027", health)
cesd_vars <- intersect(paste0("dc", sprintf("%03d", 9:18)), names(health))
adl_vars  <- intersect(paste0("db", sprintf("%03d", 1:15)), names(health))
chr_vars  <- intersect(paste0("da007_", 1:14, "_"), names(health))

health <- health %>% mutate(
  tics_score = rowSums(across(all_of(c(imm, del)), ~ as.numeric(.)), na.rm = TRUE),
  cesd_score = rowSums(across(all_of(cesd_vars), ~ as.numeric(.)), na.rm = TRUE),
  func_score = rowSums(across(all_of(adl_vars), ~ 5 - as.numeric(.)), na.rm = TRUE),
  chronic_count = rowSums(across(all_of(chr_vars), ~ as.numeric(.) == 1), na.rm = TRUE),
  self_health = coalesce(as.numeric(da001), as.numeric(da002)),
  smoking = if_else(as.numeric(da059) == 1, 1, 0))
health$alcohol <- if (!is.null(alc_var)) if_else(as.numeric(health[[alc_var]]) == 1, 1, 0) else NA_real_

demog <- demog %>% mutate(age = 2011 - as.numeric(ba002_1), sex = as.numeric(rgender),
                          education = as.numeric(bd001))
blood <- blood %>% mutate(hba1c = as.numeric(newhba1c), crp = as.numeric(newcrp),
                          cystatin_c = as.numeric(cystatinc), total_cholesterol = as.numeric(newcho))
weights <- weights %>% mutate(weight_blood = as.numeric(bio_weight1))

charls_raw <- demog %>% select(ID, age, sex, education) %>%
  left_join(health %>% select(ID, self_health, cesd_score, tics_score, func_score,
                              chronic_count, smoking, alcohol), by = "ID") %>%
  left_join(blood %>% select(ID, hba1c, crp, cystatin_c, total_cholesterol), by = "ID") %>%
  left_join(weights %>% select(ID, weight_blood, communityID), by = "ID") %>%
  filter(age >= 45 & age <= 110)
stopifnot(all(nchar(charls_raw$ID) == 12), !anyDuplicated(charls_raw$ID))

# ---- death linkage (three waves) ----
read_death <- function(dp, src) {
  if (is.na(dp) || !dir.exists(dp)) return(NULL)
  ef <- list.files(dp, "Exit_Interview.*\\.dta$", ignore.case = TRUE, full.names = TRUE)
  if (length(ef) > 0) { df <- read_dta(ef[1])
    yc <- intersect(c("exb001_1","exb0011","death_year"), names(df))
    mc <- intersect(c("exb001_2","exb0012","death_month"), names(df))
    if (length(yc) > 0) return(df %>% transmute(ID = as.character(ID), died = 1,
      death_year = as.numeric(.data[[yc[1]]]),
      death_month = if (length(mc) > 0) as.numeric(.data[[mc[1]]]) else 7) %>% filter(!is.na(death_year))) }
  sf <- list.files(dp, "Sample_Infor.*\\.dta$", ignore.case = TRUE, full.names = TRUE)
  if (length(sf) > 0) { df <- read_dta(sf[1])
    dv <- grep("^(died|death|dead)", names(df), ignore.case = TRUE, value = TRUE)
    if (length(dv) > 0) return(df %>% filter(as.numeric(.data[[dv[1]]]) == 1) %>%
      transmute(ID = as.character(ID), died = 1,
        death_year = as.numeric(gsub(".*(20[0-9]{2}).*", "\\1", src)), death_month = 7)) }
  NULL }
d13 <- read_death(year_dir["w2"], "exit2013"); d15 <- read_death(year_dir["w3"], "exit2015")
d18 <- read_death(year_dir["w4"], "exit2018")
all_deaths <- bind_rows(d13, d15, d18) %>% group_by(ID) %>%
  slice_min(death_year, n = 1, with_ties = FALSE) %>% ungroup() %>%
  filter(!is.na(death_year), !is.na(death_month), death_year >= 1950, death_year <= 2020,
         death_month >= 1, death_month <= 12)
all_deaths$ID <- as.character(all_deaths$ID)

safe_date <- function(y, m) { out <- rep(as.Date(NA), length(y))
  ok <- !is.na(y) & !is.na(m) & y >= 1950 & y <= 2020 & m >= 1 & m <= 12
  if (any(ok)) out[ok] <- as.Date(sprintf("%04d-%02d-15", as.integer(y[ok]), as.integer(m[ok])))
  out }
charls_surv <- charls_raw %>% left_join(all_deaths, by = "ID") %>%
  mutate(died = if_else(is.na(died), 0L, 1L),
         death_date = safe_date(death_year, death_month),
         baseline_date = as.Date("2011-07-01"), censor_date = as.Date("2018-12-31"),
         last_date = if_else(died == 1L & !is.na(death_date), death_date, censor_date),
         surv_time = as.numeric(difftime(last_date, baseline_date, units = "days"))/365.25) %>%
  filter(surv_time > 0) %>% filter(died == 0L | !is.na(death_date))
cat("N =", nrow(charls_surv), " deaths =", sum(charls_surv$died), "\n")

# ---- Nelson-Aalen + mice m=20 ----
cox0 <- coxph(Surv(surv_time, died) ~ age + sex, data = charls_surv)
sf1 <- survfit(Surv(surv_time, died) ~ 1, data = charls_surv)
charls_surv$nelson <- approx(sf1$time, sf1$cumhaz, xout = charls_surv$surv_time,
                             method = "constant", yleft = 0, rule = 2, f = 0)$y
imp_vars <- c("age","sex","education","self_health","cesd_score","tics_score","func_score",
              "hba1c","crp","cystatin_c","total_cholesterol","chronic_count","smoking",
              "alcohol","died","nelson")
id_vec <- charls_surv$ID
set.seed(2024)
imp <- mice(charls_surv %>% select(all_of(imp_vars)), m = 20, maxit = 20,
            method = "pmm", printFlag = FALSE)
pdf(file.path(OUT_DIR, "mice_trace_m20.pdf"), width = 12, height = 10)
plot(imp, layout = c(4,4)); dev.off()

build_vae_input <- function(d) d %>% mutate(
  self_health_rev = 6 - self_health, cesd_rev = -cesd_score, tics_cog = tics_score,
  physio_proxy = func_score, log_crp = log(crp + 1), hba1c_rev = -hba1c,
  crp_rev = -log_crp, cystatin_c_rev = -cystatin_c, total_chol_rev = -total_cholesterol) %>%
  select(tics_cog, self_health_rev, cesd_rev, physio_proxy,
         hba1c_rev, crp_rev, cystatin_c_rev, total_chol_rev)

for (i in 1:20) {
  imp_i <- complete(imp, i) %>% mutate(ID = id_vec) %>%
    left_join(charls_surv %>% select(ID, weight_blood, communityID,
                                     surv_time, death_year, death_month), by = "ID")
  z_i <- scale(as.matrix(build_vae_input(imp_i)))
  stopifnot(!any(is.na(z_i)))
  write.csv(data.frame(ID = imp_i$ID), file.path(PROC_DIR, sprintf("charls_ids_m%02d.csv", i)),
            row.names = FALSE, quote = TRUE)
  write.csv(z_i, file.path(PROC_DIR, sprintf("charls_std_for_vae_m%02d.csv", i)), row.names = FALSE)
  if (i == 1) {
    write.csv(data.frame(ID = imp_i$ID), file.path(PROC_DIR, "charls_ids.csv"),
              row.names = FALSE, quote = TRUE)
    write.csv(z_i, file.path(PROC_DIR, "charls_std_for_vae.csv"), row.names = FALSE)
    write.csv(imp_i %>% select(ID, age, sex, education, self_health, cesd_score, tics_score,
        func_score, chronic_count, smoking, alcohol, weight_blood, communityID,
        surv_time, died, death_year, death_month),
        file.path(PROC_DIR, "charls_covars.csv"), row.names = FALSE, quote = TRUE)
  }
}
cat("01b complete: 20 imputation sets + main files exported\n")
