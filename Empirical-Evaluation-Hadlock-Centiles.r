
########## THIS R SCRIPT WAS USED FOR EMPIRICAL VALIDATION OF THE HADLOCK CHART ##########
####### BASED ON AN OBSERVED CLINCICAL DATASET #########################

# Hadlock equation parameters
meanEFW <- function(w){
  exp(0.578 + 0.332 * w - 0.00354 * w^2)
}

sd_logEFW <- 0.12


p_vec <- c(0.03, 0.10, 0.50, 0.90, 0.97)
w <- 24

efw_p <- sapply(p_vec, function(p) {
  mu_log <- log(meanEFW(w))
  z <- qnorm(p)
  exp(mu_log + z * sd_logEFW)
})

names(efw_p) <- paste0(p_vec*100, "%")
efw_p


# recreate analytic_simple if needed
weeks <- 10:40
percentiles <- c(0.03, 0.10, 0.50, 0.90, 0.97)

# Create empty dataframe with the right structure
correct_wide <- data.frame(WGA = weeks)

# Add each percentile column
for(p in percentiles) {
  # Calculate EFW for this percentile across all weeks
  correct_wide[[paste0("p_", p*100)]] <- sapply(weeks, function(w) {
    exp(log(meanEFW(w)) + qnorm(p) * sd_logEFW)
  })
}

# Rename columns to match wrong_chart format
colnames(correct_wide) <- c("WGA", "3rd", "10th", "50th", "90th", "97th")

View(correct_wide)

#The standard Hadlock growth chart
wrong_chart <- read.csv(text = 'WGA,3rd,10th,50th,90th,97th
10,26,29,35,41,44
11,34,37,45,53,56
12,43,48,58,68,73
13,55,61,73,85,91
14,70,77,93,109,116
15,88,97,117,137,146
16,110,121,146,171,183
17,136,150,181,212,226
18,167,185,223,261,279
19,205,227,273,319,341
20,248,275,331,387,414
21,299,331,399,467,499
22,359,398,478,559,598
23,426,471,568,665,710
24,503,556,670,784,838
25,589,652,785,918,981
26,685,758,913,1068,1141
27,791,876,1055,1234,1319
28,908,1004,1210,1416,1513
29,1034,1145,1379,1613,1724
30,1169,1294,1559,1824,1649
31,1313,1453,1751,2049,2189
32,1465,1621,1953,2285,2441
33,1622,1794,2162,2530,2703
34,1783,1973,2377,2781,2971
35,1946,2154,2595,3036,3244
36,2110,2335,2813,3291,3516
37,2271,2513,3028,3543,3785
38,2427,2686,3236,3786,4045
39,2576,2851,3435,4019,4294
40,2714,3004,3619,4234,4524')

head(wrong_chart)

write_xlsx(wrong_chart, "wrong_chart.xlsx")

# Remove any "X" prefixes from both datasets
clean_colnames <- function(df) {
  colnames(df) <- gsub("^X", "", colnames(df))
  return(df)
}

wrong_chart <- clean_colnames(wrong_chart)
correct_wide <- clean_colnames(correct_wide)

print("After cleaning:")
print(colnames(wrong_chart))
print(colnames(correct_wide))

library(dplyr)
library(tidyr)

# If the merge didn't work as expected, let's recreate it with explicit naming
library(dplyr)
library(tidyr)

# Clear approach with explicit column naming
wrong_long <- wrong_chart %>%
  pivot_longer(
    cols = c(`3rd`, `10th`, `50th`, `90th`, `97th`),
    names_to = "percentile", 
    values_to = "efw_wrong"
  )

correct_long <- correct_wide %>%
  pivot_longer(
    cols = c(`3rd`, `10th`, `50th`, `90th`, `97th`),
    names_to = "percentile",
    values_to = "efw_correct"
  )

# Merge with explicit column handling
merged_data_clean <- wrong_long %>%
  inner_join(correct_long, by = c("WGA", "percentile")) %>%
  mutate(
    bias = efw_correct - efw_wrong
  )

head(merged_data_clean)

# Use this clean version for the PFFR analysis
merged_data <- merged_data_clean

# PFFR preparation
bias_wide <- merged_data %>%
  select(WGA, percentile, bias) %>%
  pivot_wider(
    names_from = WGA,
    values_from = bias,
    names_prefix = "wga_"
  )

print("Bias data in wide format:")
print(dim(bias_wide))
head(bias_wide[, 1:6])



# Define the corrected model's CDF
Phi <- pnorm
F_corr <- function(v, w){
  mu <- log(meanEFW(w))
  Phi((log(v) - mu) / sd_logEFW)
}

# Compute corrected percentiles for each published value
corrected_p <- data.frame(WGA = wrong_chart$WGA)

for (pct in c("3rd","10th","50th","90th","97th")){
  corrected_p[[paste0(pct, "_pCorr")]] <-
    mapply(function(v, w) 100 * F_corr(v, w), 
           v = wrong_chart[[pct]], w = wrong_chart$WGA)
}

# Round for presentation
corrected_p <- corrected_p %>%
  mutate(across(-WGA, round, 0))

View(corrected_p)


write_xlsx(corrected_p, "Figure.xlsx")

# --- prerequisites: define meanEFW and sd_logEFW if not already ---
meanEFW <- function(w) exp(0.578 + 0.332*w - 0.00354*w^2)  # corrected model
sd_logEFW <- 0.12
z10 <- qnorm(0.10)

# weeks to evaluate (must match wrong_chart$WGA)
weeks <- wrong_chart$WGA

# Published percentile points available in the table
pub_pcts <- c(3,10,50,90,97)
pub_cols  <- c("3rd","10th","50th","90th","97th")  # adjust if needed

# Function to compute equivalent published percentile for a given week's corrected p (here 10%)
find_equiv_pct_for_true_p <- function(w, true_p = 0.10){
  # true EFW at true_p under corrected model
  v_true <- exp(log(meanEFW(w)) + qnorm(true_p) * sd_logEFW)
  
  # published EFW values for that week
  row_idx <- which(wrong_chart$WGA == w)
  if(length(row_idx) != 1) stop("Week ", w, " not found uniquely in wrong_chart")
  v_pub <- as.numeric(wrong_chart[row_idx, pub_cols])
  
  # if any NA or non-monotonic adjust
  if(any(is.na(v_pub))) {
    warning("NA in published values for week ", w, "; returning NA")
    return(NA_real_)
  }
  
  # interpolation on LOG scale (more appropriate for weights)
  log_v_pub <- log(v_pub)
  log_v_true <- log(v_true)
  
  # Make sure published percentiles are strictly increasing in value
  # If not strictly increasing, slightly jitter or fall back to linear approx on raw scale
  if(!all(diff(log_v_pub) > 0)){
    # fallback: use approx on raw values
    p_equiv <- approx(x = v_pub, y = pub_pcts, xout = v_true, rule = 2)$y
  } else {
    p_equiv <- approx(x = log_v_pub, y = pub_pcts, xout = log_v_true, rule = 2)$y
  }
  return(as.numeric(p_equiv))
}

# Apply across weeks for true 10th percentile
equiv_pub_pct_for_true10 <- sapply(weeks, find_equiv_pct_for_true_p, true_p = 0.10)

# Build a result table
result <- data.frame(
  WGA = weeks,
  true10_EFW_corr = round(exp(log(meanEFW(weeks)) + z10 * sd_logEFW)),  # grams
  pub_equiv_pct = round(equiv_pub_pct_for_true10, 0)
)

# Round/format the percentile
#result$pub_equiv_pct_round <- round(result$pub_equiv_pct, 0)  
result


install.packages("kableExtra")

install.packages("writexl")
library(writexl)
write_xlsx(result, "Hadlock_true10_equivalence_table.xlsx")


library(writexl)

# Simple direct matching without log transformation
direct_percentile_match <- function(w, true_p = 0.10) {
  true_10th <- exp(log(meanEFW(w)) + qnorm(true_p) * sd_logEFW)
  
  row_idx <- which(wrong_chart$WGA == w)
  v_pub <- as.numeric(wrong_chart[row_idx, pub_cols])
  
  # Direct linear interpolation on raw scale
  p_equiv <- approx(x = v_pub, y = pub_pcts, xout = true_10th, rule = 2)$y
  
  return(p_equiv)
}

direct_results <- sapply(weeks, direct_percentile_match, true_p = 0.10)

View(direct_results)

# Smooth the published curves first, then interpolate
library(mgcv)

smooth_interpolation_method <- function() {
  results <- numeric(length(weeks))
  
  # Smooth each published percentile curve
  smoothed_curves <- list()
  for(i in 1:length(pub_cols)) {
    gam_fit <- gam(wrong_chart[[pub_cols[i]]] ~ s(WGA, k = 10), data = wrong_chart)
    smoothed_curves[[i]] <- predict(gam_fit, newdata = data.frame(WGA = weeks))
  }
  
  for(j in 1:length(weeks)) {
    w <- weeks[j]
    true_10th <- exp(log(meanEFW(w)) + qnorm(0.10) * sd_logEFW)
    
    # Get smoothed published values at this WGA
    v_pub_smooth <- sapply(smoothed_curves, function(x) x[j])
    
    # Interpolate on log scale
    p_equiv <- approx(x = log(v_pub_smooth), y = pub_pcts, 
                      xout = log(true_10th), rule = 2)$y
    
    results[j] <- p_equiv
  }
  
  return(results)
}

smooth_results <- smooth_interpolation_method()

View(smooth_results)

setwd("~/Desktop/research work/Hadlock percentile")


#Importing the empirical dataset
dataset <- read.csv("Working clean dataset.csv")

dataset <- dataset[, c(
  "mrn_m",
  "thirtysix_1_examdate",
  "thirtysix_1_usgestationw",
  "thirtysix_1_usgestationd",
  "thirtysix_1_estweight",
  "thirtysix_1_estweight_centile",
    "thirtysix_1_hc", "thirtysix_1_ac", "thirtysix_1_fl",
    "mort_ga_36plus",
    "sb_ga_36plus",
  "peri_ga_36plus",
  "nd_ga_36plus",
  "ph_71",
  "apgar_7",
  "sepsis_confirmed_summary",
  "respiratory_distress_summary",
  "nnu_admissions"
)]

dataset <- dataset[
  !is.na(dataset$thirtysix_1_estweight) &
  !is.na(dataset$thirtysix_1_usgestationw),
]


nrow(dataset)

colnames(dataset)

dataset$WGA <- dataset$thirtysix_1_usgestationw

dataset$GA_decimal <- dataset$thirtysix_1_usgestationw + dataset$thirtysix_1_usgestationd / 7

dataset_valid <- dataset[dataset$WGA %in% 35:36, ]

merged$thirtysix_1_fl <- dataset$thirtysix_1_fl
merged$thirtysix_1_hc <- dataset$thirtysix_1_hc
merged$thirtysix_1_ac <- dataset$thirtysix_1_ac

interp_10th <- approxfun(wrong_chart$WGA, wrong_chart$`10th`, rule = 2)

merged <- merge(dataset_valid, wrong_chart, by = "WGA", all.x = TRUE)

merged$P10_interp <- interp_10th(merged$GA_decimal)
merged$SGA10_interp <- merged$thirtysix_1_estweight < merged$P10_interp
num <- sum(merged$SGA10_interp, na.rm = TRUE)
prop <- mean(merged$SGA10_interp, na.rm = TRUE)*100

cat("The number of SGA in the cohort is", num, "\n")
cat("The proportion of EFW below 10th centile in the cohort is", prop,"%", "\n")
week_summary_interp <- merged %>%
  filter(WGA %in% c(35, 36)) %>%
  group_by(WGA) %>%
  summarise(
    n = n(),
    n_SGA10 = sum(SGA10_interp, na.rm = TRUE),
    prop_SGA10 = mean(SGA10_interp, na.rm = TRUE),
    percentage = prop_SGA10 * 100
  )

week_summary_interp

x <- 917
n <- 21874
binom.test(x,n)$conf.int

interp_3rd <- approxfun(wrong_chart$WGA, wrong_chart$`3rd`, rule = 2)
merged$P3_interp <- interp_3rd(merged$GA_decimal)
merged$SGA3_interp <- merged$thirtysix_1_estweight < merged$P3_interp
num <- sum(merged$SGA3_interp, na.rm = TRUE)
prop <- mean(merged$SGA3_interp, na.rm = TRUE)*100

cat("The number of SGA in the cohort is", num, "\n")
cat("The proportion of EFW below 10th centile in the cohort is", prop,"%", "\n")
week_summary_interp <- merged %>%
  filter(WGA %in% c(35, 36)) %>%
  group_by(WGA) %>%
  summarise(
    n = n(),
    n_SGA3 = sum(SGA3_interp, na.rm = TRUE),
    prop_SGA3 = mean(SGA3_interp, na.rm = TRUE),
    percentage = prop_SGA3 * 100
  )

week_summary_interp

binom.test(167, 21874)$conf.int

# build interpolation function for corrected chart (10th or 3rd)
interp_10th_corr <- approxfun(correct_wide$WGA, correct_wide$`10th`, rule = 2)
interp_3rd_corr  <- approxfun(correct_wide$WGA, correct_wide$`3rd`,  rule = 2)

# decimal GA (weeks + days)
merged$GA_decimal <- merged$thirtysix_1_usgestationw + merged$thirtysix_1_usgestationd / 7

# interpolated centiles
merged$P10_corr_interp <- interp_10th_corr(merged$GA_decimal)
merged$P3_corr_interp  <- interp_3rd_corr(merged$GA_decimal)

# classification using interpolated centiles
merged$SGA10_corr_interp <- merged$thirtysix_1_estweight < merged$P10_corr_interp
merged$SGA3_corr_interp  <- merged$thirtysix_1_estweight < merged$P3_corr_interp

# totals and week-by-week (35 & 36 example)
cat("Number of EFW <10 (interp):", sum(merged$SGA10_corr_interp, na.rm=TRUE), "\n")
cat("Proportion of EFW <10 (interp):", mean(merged$SGA10_corr_interp, na.rm=TRUE)*100, "\n")
cat("Number of EFW <3 (interp):", sum(merged$SGA3_corr_interp, na.rm=TRUE), "\n")
cat("Proportion of EFW <3 (interp):", mean(merged$SGA3_corr_interp, na.rm=TRUE)*100, "\n")

week_summary_interp_corr <- merged %>%
  filter(WGA %in% c(35,36)) %>%
  group_by(WGA) %>%
  summarise(
    n = n(),
    n_SGA10_interp = sum(SGA10_corr_interp, na.rm = TRUE),
    prop_SGA10_interp = mean(SGA10_corr_interp, na.rm = TRUE),
    n_SGA3_interp = sum(SGA3_corr_interp, na.rm = TRUE),
    prop_SGA3_interp = mean(SGA3_corr_interp, na.rm = TRUE),
  )

week_summary_interp_corr

#10th centile
binom.test(1577, 21874)$conf.int

#3rd centile
binom.test(502, 21874)$conf.int

merged$peri_ga_36plus <- dataset$peri_ga_36plus

merged$peri_or_sb <- as.integer(
  merged$sb_ga_36plus == 1 | merged$peri_ga_36plus == 1
)

n_stillbirth <- sum(merged$peri_or_sb == 1, na.rm = TRUE)
n_stillbirth


sb_raw_highrisk <- sum(
  merged$peri_or_sb == 1 & merged$SGA10_interp == TRUE,
  na.rm = TRUE
)
sb_raw_highrisk


sb_corr_highrisk <- sum(
  merged$peri_or_sb == 1 & merged$SGA10_corr_interp == TRUE,
  na.rm = TRUE
)
sb_corr_highrisk


# build interpolation functions for multiple centiles
interp_3  <- approxfun(wrong_chart$WGA, wrong_chart$`3rd`,  rule = 2)
interp_10 <- approxfun(wrong_chart$WGA, wrong_chart$`10th`, rule = 2)
interp_50 <- approxfun(wrong_chart$WGA, wrong_chart$`50th`, rule = 2)

sb <- merged[merged$peri_or_sb == 1, ]

sb$P3  <- interp_3(sb$GA_decimal)
sb$P10 <- interp_10(sb$GA_decimal)
sb$P50 <- interp_50(sb$GA_decimal)

# classify which centile band each stillbirth falls into
sb$band <- with(sb, ifelse(thirtysix_1_estweight < P3, "<3rd",
                    ifelse(thirtysix_1_estweight < P10, "3rd–10th",
                    ifelse(thirtysix_1_estweight < P50, "10th–50th", "≥50th"))))

sb[, c("WGA","GA_decimal","thirtysix_1_estweight","band")]


dataset$severe_morbidity <- with(dataset,
  as.integer(
    peri_ga_36plus == 1 |
    nd_ga_36plus == 1 |
    ph_71 == 1 |
    apgar_7 == 1 |
    sepsis_confirmed_summary == 1 |
    respiratory_distress_summary == 1 |
    nnu_admissions == 1
  )
)


# components to include
components <- c("ph_71","apgar_7",
                "respiratory_distress_summary","nnu_admissions")

# coerce to numeric 0/1 (safe)
dataset[components] <- lapply(dataset[components], function(x) {
  as.numeric(as.character(x))
})

# create composite with the robust rule
dataset$severe_morbidity <- apply(dataset[components], 1, function(row) {
  if(all(is.na(row))) return(NA_integer_)
  if(any(row == 1, na.rm = TRUE)) return(1L)
  return(0L)
})

# diagnostics
table(dataset$severe_morbidity, useNA = "ifany")




merged$severe_morbidity <- dataset$severe_morbidity

comp_corr_highrisk <- sum(
  merged$severe_morbidity == 1 & merged$SGA10_corr_interp == TRUE,
  na.rm = TRUE
)
comp_corr_highrisk


comp_raw_highrisk <- sum(
  merged$severe_morbidity == 1 & merged$SGA10_interp == TRUE,
  na.rm = TRUE
)
comp_raw_highrisk


merged[c(4860, 11226, 8002, 13086, 8506 ), ]


# one-off numbers (from your wrong_chart)
WGA_target <- 36
EFW_target <- 2385

row36 <- wrong_chart[wrong_chart$WGA == WGA_target, ]
centiles <- c(3, 10, 50, 90, 97)
weights  <- as.numeric(row36[c("3rd","10th","50th","90th","97th")])

# function mapping weight -> centile by linear interpolation of the table at that GA
map_weight_to_centile <- approxfun(x = weights, y = centiles, rule = 2)

# get centile
centile_chart <- map_weight_to_centile(EFW_target)
centile_chart


# Equation based chart
wga <- 36  # weeks gestational age
efw_obs <- 2385  # observed estimated fetal weight in grams

# Compute mu (mean in log-scale)
mu <- 0.578 + 0.332 * wga - 0.00354 * wga^2

# Constant SD in log-scale
sd_val <- 0.12

# Z-score
z <- (log(efw_obs) - mu) / sd_val

# Percentile (multiply by 100 for %)
percentile <- pnorm(z) * 100

# Outputs
cat("Percentile:", percentile, "\n")

sum(merged$severe_morbidity == 1 & merged$SGA10_corr_interp == TRUE, na.rm = TRUE)
mean(merged$SGA10_corr_interp[merged$severe_morbidity == 1], na.rm = TRUE)



sum(merged$severe_morbidity == 1)

contingency_table <- matrix(c(60, 33,0, 1197), nrow = 2, byrow = TRUE)
contingency_table

mcnemar.test(contingency_table, correct = FALSE)

mcnemar.test(table(merged$SGA10_interp, merged$SGA10_corr_interp))





