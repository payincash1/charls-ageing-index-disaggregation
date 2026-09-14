# 11_tics_check.R — verify tics_score coding vs correct word-recall count
rm(list = ls()); library(dplyr); library(haven)
base <- "/Users/a1234/Desktop/桌面/chalrs"   # <-- EDIT
h <- read_dta(list.files(file.path(base, "chalrs2011"),
  pattern = "health_status.*\\.dta$", ignore.case = TRUE, full.names = TRUE)[1])
cov <- read.csv(file.path(base, "processed/charls_covars.csv"), colClasses = c(ID = "character"))
mk <- function(raw_id) {
  s <- gsub('"', '', gsub("\\s", "", as.character(raw_id)), fixed = TRUE)
  len <- nchar(s)
  if (all(len[!is.na(len)] == 12)) return(s)
  paste0(sprintf("%09d", suppressWarnings(as.numeric(substr(s, 1, pmax(len - 2, 0))))), "0",
         substr(s, len - 1, len)) }
h$ID <- mk(h$ID)
imm <- paste0("dc006s", 1:10); del <- paste0("dc027s", 1:10)
h2 <- h %>% mutate(
  tics_correct = rowSums(!is.na(across(all_of(c(imm, del))))),
  cesd_std = rowSums(across(all_of(paste0("dc", sprintf("%03d", 9:18))),
                          ~ as.numeric(.)), na.rm = TRUE) - 10) %>% select(ID, tics_correct, cesd_std)
cmp <- cov %>% select(ID, tics_score, cesd_score, died) %>% left_join(h2, by = "ID")
cat("matched:", sum(!is.na(cmp$tics_correct)), "/", nrow(cmp), "\n")
cat("cor(tics):", round(cor(cmp$tics_score, cmp$tics_correct, use = "complete.obs"), 3), "\n")
cat("cor(cesd):", round(cor(cmp$cesd_score, cmp$cesd_std, use = "complete.obs"), 3), "\n")
cmp %>% group_by(died) %>% summarise(tics_mean = round(mean(tics_correct, na.rm = T), 1),
  tics_sd = round(sd(tics_correct, na.rm = T), 1), cesd_mean = round(mean(cesd_std, na.rm = T), 1),
  cesd_sd = round(sd(cesd_std, na.rm = T), 1), .groups = "drop") %>% print()
