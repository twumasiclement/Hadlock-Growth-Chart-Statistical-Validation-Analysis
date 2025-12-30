# hadlock_efw_utils.R
# Utilities to compute Hadlock (1991) "optimal" model percentiles and bootstrap CIs
# Follows Roberts et al. (2025) description:
#   MeanEFW = exp(0.578 + 0.332*WGA - 0.00354*WGA^2)
#   SD (on log scale) = 0.12
#   pth percentile = exp( ln(MeanEFW) + z_p * SD )

####Setting the plot size and resolution (300 dpi)
options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(fda);         # functional data analysis
library(Hotelling)    # Hotelling's T2 test (paired)
library(mvtnorm)      # for multivariate simulation if needed


# A function to compute the  mean EFW based on the optimal model by Hadlock (1991)
hadlock_mean_efw <- function(wga){
  # wga: numeric vector of week gestational ages (WGA)
  mu_log <- 0.578 + 0.332 * wga - 0.00354 * (wga^2)
  mean_efw <- exp(mu_log)
  return(mean_efw)
}

# Testing the function at WGA=24 weeks for example 
hadlock_mean_efw(wga=24)

# Not needed, just for exploration
out<- data.frame(WGA=10:40, EFW=hadlock_mean_efw(wga=10:40))
head(out)

# Reproducing the optimal regression coefficients
mod <- lm(log(EFW) ~ WGA + I(WGA^2), data = out)
summary(mod)
coef(mod)


 mu_log <- 0.578 + 0.332 * wga - 0.00354 * (wga^2)

# A function to compute the percentiles of EFW at varying values of WGA, assuming SD of log(EFW)=0.12
hadlock_percentile <- function(wga, p = 50, sd_log = 0.12){
  # Compute a single percentile (or vectorised) from the Hadlock model.
  # wga: numeric vector of WGA
  # p: percentile(s) in percent (e.g., 3, 10, 50, 90, 97) - can be vectorised
  # sd_log: standard deviation on log(EFW) scale (default 0.12)
  stopifnot(all(p > 0 & p < 100))
  mean_efw <- hadlock_mean_efw(wga)
  # z for percentile
  z <- qnorm(p/100, lower.tail=TRUE)
  # vectorisation: expand grid if lengths differ
  # If wga length > 1 and p length > 1, return matrix-like expanded result via outer:
  if(length(wga) > 1 && length(p) > 1){
    # returns a data.frame in long format. Below convenience wrapper handles this.
    m <- outer(log(mean_efw), z, FUN = function(m, z) exp(m + z * sd_log))
    colnames(m) <- paste0(p, "th")
    rownames(m) <- as.character(wga)
    return(m)
  } else {
    # simple vectorised calculation
    mu_log <- log(mean_efw)
    res <- exp(mu_log + z * sd_log)
    return(res)
  }
}

# Testing the function at WGA=10 weeks for example 
percentiles_p<- c(3, 10, 50, 90, 97)

print(hadlock_percentile(wga=40, p = percentiles_p, sd_log = 0.12))


print(ceiling(signif(hadlock_percentile(wga= 40, p = percentiles_p, sd_log = 0.12), 3)))

print(round(signif(hadlock_percentile(wga=40, p = percentiles_p, sd_log = 0.12), 3)))

print(round(hadlock_percentile(wga=40, p = percentiles_p, sd_log = 0.12)))


# A function to return WGA, percentile, model-based EFW
hadlock_percentiles_df <- function(wga = seq(24, 40, by = 0.5),
                                   percentiles = c(3, 10, 50, 90, 97),
                                   sd_log = 0.12,
                                   sigfig = 3){
  # Return a tidy data.frame with columns: WGA, percentile, efw_model
  wga <- as.numeric(wga)
  ps <- sort(unique(percentiles))
  mu_log <- 0.578 + 0.332 * wga - 0.00354 * (wga^2)
  
  out <- do.call(rbind, lapply(seq_along(wga), function(i){
    z <- qnorm(ps / 100)
    efw_vals <- exp(mu_log[i] + z * sd_log)
    data.frame(WGA = rep(wga[i], length(ps)),
               percentile = ps,
               efw_model = round(signif(efw_vals, sigfig)),   # round only at output
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  return(out)
}

head(hadlock_percentiles_df(wga = 10:40, percentiles = c(3,10,50,90,97)), n=10)

# A function to construct bootstrap confidence intervals of the mean EFW
bootstrap_percentile_CI <- function(wga,
                                    percentile = 50,
                                    n = 100,
                                    B = 2000,
                                    sd_log = 0.12,
                                    conf = 0.95,
                                    seed = NULL,
                                    return_bootstrap_samples = FALSE){
  # Parametric bootstrap to estimate CI for the sample percentile (for a given sample size n).
  # wga: single numeric WGA value (if vector provided, function will vectorise via sapply/lapply externally)
  # percentile: p in percent (single value)
  # n: sample size per bootstrap sample (integer)
  # B: number of bootstrap replicates
  # sd_log: sd on log scale (default 0.12)
  # conf: confidence level for interval (default 0.95)
  # seed: optional integer
  # return_bootstrap_samples: if TRUE will return bootstrap vector (makes result heavier)
  stopifnot(length(wga) == 1)
  if(!is.null(seed)) set.seed(seed)
  mu_log <- 0.578 + 0.332 * wga - 0.00354 * (wga^2)
  # simulate B bootstrap samples, each with n observations drawn from log-EFW distribution:
  # ln(EFW) ~ N(mu_log, sd_log^2), EFW = exp(ln(EFW))
  # For each bootstrap sample, compute the sample percentile (type=7 default)
  boot_pers <- numeric(B)
  for(b in seq_len(B)){
    samp <- rnorm(n = n, mean = mu_log, sd = sd_log)
    efw_samp <- exp(samp)
    # Percentile based on the theoretical formula of Hadlock and bootstrap samples (Not advisable to introduce variability)
    # boot_pers[b] <- hadlock_percentile(wga=wga, p = percentile, sd_log =sd(samp))
     boot_pers[b] <- as.numeric(quantile(efw_samp, probs = percentile/100, type = 7))
  }
   alpha_over_2 <- (1 - conf) / 2
   lower <- quantile(boot_pers, probs = alpha_over_2)
   upper <- quantile(boot_pers, probs = 1 - alpha_over_2)
   #  # # # #bootstrap normal approximation CI
   # boot_mean <- mean(boot_pers)
   # boot_sd   <- sd(boot_pers)
  
  # lower <- boot_mean - qnorm(1 - alpha_over_2) * (boot_sd / sqrt(B))
  # upper <- boot_mean + qnorm(1 -alpha_over_2) * (boot_sd / sqrt(B))
    
  # also compute model (Roberts/Hadlock) analytic percentile for reference
  model_pred <- as.numeric(exp(mu_log + qnorm(percentile/100) * sd_log))
  out <- list(WGA = wga,
              percentile = percentile,
              model_pred = model_pred,
              boot_mean = mean(boot_pers),
              CI_lower = as.numeric(lower),
              CI_upper = as.numeric(upper),
              B = B,
              n = n)
  if(return_bootstrap_samples) out$bootstrap_samples <- boot_pers
  return(out)
}


# out<- bootstrap_percentile_CI(wga=10, percentile = 50,
                                    # n = 100,
                                    # B = 2000,
                                    # sd_log = 0.12,
                                    # conf = 0.95,
                                    # seed = 123,
                                    # return_bootstrap_samples = TRUE)

# hist(out$bootstrap_samples, col="blue")

# A function to return  WGA, percentile, model-based EFW, bootstrap mean EFW, 95% Prediction Intervals, 
# bootstrap sample size (n), and  number of bootstrap replicates (B)
bootstrap_percentiles_grid <- function(wga = 10:40,
                                       percentiles = c(3, 10, 50, 90, 97),
                                       n = 100,
                                       B = 2000,
                                       sd_log = 0.12,
                                       conf = 0.95,
                                       sigfig = 3,
                                       seed = NULL,
                                       parallel = FALSE){
  # Run bootstrap_percentile_CI across a grid of WGA and percentiles.
  # Returns tidy data.frame: WGA, percentile, model_pred, boot_mean, CI_lower, CI_upper, n, B
  if(!is.null(seed)) set.seed(seed)
  wga <- as.numeric(wga)
  ps <- sort(unique(percentiles))
  results <- vector("list", length = length(wga) * length(ps))
  idx <- 1
  for(i in seq_along(wga)){
    for(j in seq_along(ps)){
      res <- bootstrap_percentile_CI(wga = wga[i],
                                     percentile = ps[j],
                                     n = n,
                                     B = B,
                                     sd_log = sd_log,
                                     conf = conf,
                                     seed = NULL,
                                     return_bootstrap_samples = FALSE)
      results[[idx]] <- data.frame(WGA = res$WGA,
                                   percentile = res$percentile,
                                   model_pred_EFW = round(signif(res$model_pred,sigfig)),
                                   boot_mean_EFW = round(signif(res$boot_mean, sigfig)),
                                   Pred.Int_lower = round(signif(res$CI_lower,sigfig),1),
                                   Pred.Int_upper = round(signif(res$CI_upper,sigfig), 1),
                                   n_boot_size = res$n,
                                   N_boot_replicates = res$B,
                                   stringsAsFactors = FALSE)
      idx <- idx + 1
    }
  }
  out_df <- do.call(rbind, results)
  rownames(out_df) <- NULL
  return(out_df)
}

# for reproducibility of results
wgas<- 10:40 # 10 to 40 WGAs
percentile_values<- c(3,10,50,90,97)
# model percentiles
model_df <- hadlock_percentiles_df(wga = wgas, percentiles = percentile_values)
head(model_df)

# bootstrap grid with moderate B for demo
EFW_boot_df <- bootstrap_percentiles_grid(wga = wgas, percentiles =percentile_values,
                                     n = 200, B = 1000, sigfig = 3,seed = 1)


# look at 10th percentile results
subset(EFW_boot_df, percentile == 10)


# look at the results for WGA=24 weeks across all percentile CIs
subset(EFW_boot_df,WGA== 24)

# Create Table 2 dataframe from Hadlock 1991 critical point values
Uncorrected_table2_df<- data.frame(
  WGA = 10:40,
  `3rd` = c(26,34,43,55,70,88,110,136,167,205,248,299,359,426,503,
            589,685,791,908,1034,1169,1313,1465,1622,1783,1946,2110,
            2271,2427,2576,2714),
  `10th` = c(29,37,48,61,77,97,121,150,185,227,275,331,398,471,556,
             652,758,876,1004,1145,1294,1453,1621,1794,1973,2154,2335,
             2513,2686,2851,3004),
  `50th` = c(35,45,58,73,93,117,146,181,223,273,331,399,478,568,670,
             785,913,1055,1210,1379,1559,1751,1953,2162,2377,2595,2813,
             3028,3236,3435,3619),
  `90th` = c(41,53,68,85,109,137,171,212,261,319,387,467,559,665,784,
             918,1068,1234,1416,1613,1824,2049,2285,2530,2781,3036,3291,
             3543,3786,4019,4234),
  `97th` = c(44,56,73,91,116,146,183,226,279,341,414,499,598,710,838,
             981,1141,1319,1513,1724,1949,2189,2441,2703,2971,3244,3516,
             3785,4045,4294,4524)
)

names(Uncorrected_table2_df)[2:6]<- c("3rd", "10th", "50th", "90th", "97th")

cat("Standard Hadlock Chart:", "\n")
head(Uncorrected_table2_df, n=10)


# Convert standard Hadlock (wide) to long format
Uncorrected_long_df <- Uncorrected_table2_df %>%
  pivot_longer(
    cols = c(`3rd`, `10th`, `50th`, `90th`, `97th`),
    names_to = "percentile",
    values_to = "EFW"
  ) %>%
  mutate( percentile = as.numeric(gsub("(st|nd|rd|th)", "", percentile)))

head(Uncorrected_long_df)

EFW_boot_df$Source <- "Corrected estimates based on Hadlock 1991 regression equation"
Uncorrected_long_df$Source <- "Standard Hadlock 1991 published estimates"


combined_df <- bind_rows(
  EFW_boot_df %>% dplyr::select(WGA, percentile, model_pred_EFW, Pred.Int_lower, Pred.Int_upper, Source) %>% 
    rename(EFW = model_pred_EFW),
  Uncorrected_long_df
)


percentile_colors <- c(
  "3"   = "#1f78b4",  # deep blue
  "10"  = "#33a02c",  # green
  "50"  = "#ff7f00",  # orange
  "90"  = "#e31a1c",  # red
  "97"  = "#6a3d9a"   # purple
)

percentile_fill <- percentile_colors  # same for ribbons


################ FIGURE 1 of the Paper ############## 

library(ggplot2)

ggplot() +
  # Corrected Hadlock ribbon (only for corrected data)
  geom_ribbon(
    data = combined_df %>% filter(Source == "Corrected estimates based on Hadlock 1991 regression equation"),
    aes(x = WGA, ymin = Pred.Int_lower, ymax = Pred.Int_upper, fill = factor(percentile)),
    alpha = 0.3, colour = NA
  ) +
  
  # Lines for both datasets
  geom_line(
    data = combined_df,
    aes(x = WGA, y = EFW, colour = factor(percentile), linetype = Source),
    size = .65
  ) +
  
  # Labels
  labs(
    y = "Estimated Fetal Weight (grams)",
    x = "Gestational age (weeks)",
    colour = "Percentile:",
    fill = "Percentile:",
    linetype = "Source:"
  ) +
  
  # Line types
  scale_linetype_manual(values = c("Corrected estimates based on Hadlock 1991 regression equation" = "solid",
                                   "Standard Hadlock 1991 published estimates" = "dashed")) +
  
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "top"
  )+
theme(
  legend.position = c(0.3, 0.7),  # 0.5 = center horizontally, 0.95 = near top
  legend.justification = "center",
  legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")  # bold legend labels
)+
guides(
  colour = guide_legend(nrow = 1, byrow = TRUE),
  fill = guide_legend(nrow = 1, byrow = TRUE),
  linetype = guide_legend(nrow = 2)
)+
scale_colour_manual(values = percentile_colors,
                    labels = c("3rd", "10th", "50th", "90th", "97th")) +
scale_fill_manual(values = percentile_fill,
                  labels = c("3rd", "10th", "50th", "90th", "97th"))

library(ggplot2)

ggplot(EFW_boot_df, aes(x = WGA, colour = factor(percentile), fill = factor(percentile))) +
  # Corrected Hadlock (solid lines)
  geom_line(aes(y = model_pred_EFW), size = 0.5) +
  
  # Bootstrap 95% prediction intervals
  geom_ribbon(aes(ymin = Pred.Int_lower, ymax = Pred.Int_upper), alpha = 0.2, colour = NA) +
  
  # Standard (original) Hadlock chart – dashed lines
  geom_line(
    data = Uncorrected_long_df,
    aes(x = WGA, y = EFW, colour = factor(percentile)),
    linetype = "dashed",
    size = 0.8
  ) +
  
  labs(
    y = "EFW (grams)",
    x = "Gestational age (weeks)",
    colour = "Percentile:",
    fill = "Percentile:"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "top"
  ) +
  scale_colour_brewer(
    palette = "Set1",
    labels = c("3rd", "10th", "50th", "90th", "97th")
  ) +
  scale_fill_brewer(
    palette = "Set1",
    labels = c("3rd", "10th", "50th", "90th", "97th")
  )

ggplot(EFW_boot_df, aes(x = WGA, colour = factor(percentile), fill = factor(percentile))) +
  geom_line(aes(y = model_pred_EFW), size = .95) +
  geom_ribbon(aes(ymin = Pred.Int_lower, ymax = Pred.Int_upper), alpha = 0.2, colour = NA) +
  labs(
    # title = "Hadlock percentiles with bootstrap 95% prediction intervals",
    y = "EFW (grams)",
    x = "Gestational age (weeks)",
    colour = "Percentile:",
    fill = "Percentile:"
  ) +
  theme_minimal() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "top"
  ) +
  scale_colour_brewer(palette = "Set1",
                      labels = c("3rd", "10th", "50th", "90th", "97th")) +
  scale_fill_brewer(palette = "Set1",
                    labels = c("3rd", "10th", "50th", "90th", "97th"))
  library(ggplot2)
   # pdat <- subset(EFW_boot_df, percentile == 10)
   ggplot(EFW_boot_df, aes(x = WGA)) +
    geom_line(aes(y = model_pred_EFW)) +
     geom_ribbon(aes(ymin = Pred.Int_lower, ymax = Pred.Int_upper), alpha = 0.2) +
    labs(title = "Hadlock 10th percentile with bootstrap 95% CI",
          y = "EFW (grams)", x = "Gestational age (weeks)")+
theme_minimal() + 
facet_wrap(~ percentile) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "top")

 
head(EFW_boot_df)

################ Table 1 of the Paper ############## 


table1_df <- EFW_boot_df %>%
  dplyr::mutate(
    percentile = paste0(percentile, "th"),
    value = sprintf("%.0f (%.1f–%.1f)", model_pred_EFW, Pred.Int_lower, Pred.Int_upper)
  ) %>%
  dplyr::select(WGA, percentile, value) %>%
  tidyr::pivot_wider(names_from = percentile, values_from = value) %>%
  dplyr::arrange(WGA)

table1_df 

setwd("/Users/clementtwumasi/Desktop/Research Collaborations_Theo/Paper_HadlockChart_Evaluation_2025")
write.csv(table1_df, "CorrectedHadlockChart_PredInts_df.csv")

library(knitr)
kable(table1_df, caption = "Estimated fetal weight percentiles by gestational age (Hadlock model with bootstrap 95% CI)")


# The corrected-Hadlock table
Corrected_table2_df  <- EFW_boot_df %>%
  dplyr::mutate(
    percentile = paste0(percentile, "th"),
     value = model_pred_EFW # format like Roberts table
  ) %>%
  dplyr::select(WGA, percentile, value) %>%
  tidyr::pivot_wider(names_from = percentile, values_from = value) %>%
  dplyr::arrange(WGA)

names(Corrected_table2_df)[2]<- "3rd"

setwd("/Users/clementtwumasi/Desktop/Research Collaborations_Theo/Paper_HadlockChart_Evaluation_2025")
write.csv(Corrected_table2_df, "CorrectedHadlockChart_WithoutPredInts_df.csv")

# Create Table 2 dataframe from Hadlock 1991 critical point values
Uncorrected_table2_df<- data.frame(
  WGA = 10:40,
  `3rd` = c(26,34,43,55,70,88,110,136,167,205,248,299,359,426,503,
            589,685,791,908,1034,1169,1313,1465,1622,1783,1946,2110,
            2271,2427,2576,2714),
  `10th` = c(29,37,48,61,77,97,121,150,185,227,275,331,398,471,556,
             652,758,876,1004,1145,1294,1453,1621,1794,1973,2154,2335,
             2513,2686,2851,3004),
  `50th` = c(35,45,58,73,93,117,146,181,223,273,331,399,478,568,670,
             785,913,1055,1210,1379,1559,1751,1953,2162,2377,2595,2813,
             3028,3236,3435,3619),
  `90th` = c(41,53,68,85,109,137,171,212,261,319,387,467,559,665,784,
             918,1068,1234,1416,1613,1824,2049,2285,2530,2781,3036,3291,
             3543,3786,4019,4234),
  `97th` = c(44,56,73,91,116,146,183,226,279,341,414,499,598,710,838,
             981,1141,1319,1513,1724,1949,2189,2441,2703,2971,3244,3516,
             3785,4045,4294,4524)
)

names(Uncorrected_table2_df)[2:6]<- c("3rd", "10th", "50th", "90th", "97th")

cat("Standard Hadlock Chart:", "\n")
head(Uncorrected_table2_df, n=10)


setwd("/Users/clementtwumasi/Desktop/Research Collaborations_Theo/Paper_HadlockChart_Evaluation_2025")
write.csv(Uncorrected_table2_df, "StandardHadlockChart_df.csv")

Corrected_table2_df<- as.data.frame(Corrected_table2_df)

cat("Corrected Hadlock Chart:", "\n")
head(Corrected_table2_df, n=10)

# convert Table2 wide -> long tidy format (WGA, percentile, efw)
table2_wide_to_long <- function(df_wide){
  df_long <- df_wide %>%
    pivot_longer(cols = -WGA, names_to = "percentile", values_to = "efw") %>%
    mutate(percentile = as.character(percentile))
  return(df_long)
}

# parametric simulator for a given WGA and percentile using the Hadlock "optimal" model
# returns B simulated sample percentile estimates (EFW in grams) using ln(EFW) ~ N(mu_log, sd_log^2)
hadlock_parametric_boot <- function(WGA, percentile = 50, sd_log = 0.12, n = 100, B = 2000, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  mu_log <- 0.578 + 0.332 * WGA - 0.00354 * (WGA^2)
  z <- qnorm(percentile / 100)
  # The analytic value (one draw) is:
  analytic <- exp(mu_log + z * sd_log)
  # simulate B bootstrap percentile estimates by generating n observations per replicate,
  # computing the sample percentile from each replicate
  res <- numeric(B)
  for(b in 1:B){
    samp_ln <- rnorm(n, mean = mu_log, sd = sd_log)
    samp <- exp(samp_ln)
    res[b] <- as.numeric(quantile(samp, probs = percentile/100, type = 7))
  }
  list(analytic = analytic, boot = res)
}


head(Uncorrected_table2_df)

# Function to compute bias table with bootstrap CI
compute_bias_table_bootCI <- function(Uncorrected_df, Corrected_df,
                            percentiles = c("3rd","10th","50th","90th","97th"),
                            sd_log = 0.12, n = 100, B = 2000, seed = 123, sigfig = 3){
  # ensure both have same WGA rows and ordering
  stopifnot(all(Uncorrected_df$WGA == Corrected_df$WGA))
  WGAs <- Uncorrected_df$WGA
  out_list <- list()
  set.seed(seed)
  for(i in seq_along(WGAs)){
    w <- WGAs[i]
    for(p in percentiles){
      unc_val <- as.numeric(Uncorrected_df[i, p])
      # use parametric bootstrap to get corrected's bootstrap draws
      pnum <- as.numeric(gsub("th|rd|st|nd","", p)) # try to parse percent
      hb <- hadlock_parametric_boot(WGA = w, percentile = pnum,
                                    sd_log = sd_log, n = n, B = B, seed = NULL)
      corrected_analytic <- round(signif(hb$analytic, sigfig))
      corrected_boot <- round(signif(hb$boot, sigfig))
      # compute bias bootstrap draws (uncorrected fixed, corrected variable)
      bias_boot <- corrected_boot-unc_val 
      bias_point <- corrected_analytic-unc_val
        
      ci_lower <- quantile(bias_boot, probs = 0.025)
      ci_upper <- quantile(bias_boot, probs = 0.975)
        
      out_list[[length(out_list) + 1]] <- data.frame(
        WGA = w,
        percentile = p,
        Uncorrected = unc_val,
        Corrected = corrected_analytic,
        bias = bias_point,
        CI_lower = as.numeric(ci_lower),
        CI_upper = as.numeric(ci_upper),
        stringsAsFactors = FALSE
      )
    }
  }
  out_df <- do.call(rbind, out_list)
  # adding a column indicating significance (CI excludes 0 -> significant)
  out_df <- out_df %>%
    mutate(significant = !(CI_lower <= 0 & CI_upper >= 0))
  rownames(out_df)<- NULL
  return(out_df)
}


bias_tbl <- compute_bias_table_bootCI(Uncorrected_df=Uncorrected_table2_df, 
                    Corrected_df=Corrected_table2_df, B = 1000, n = 200, seed = 1)
head(bias_tbl)

cat("Comparing the bias between the standard and corrected Hadlock chart at the 3rd percentile:", "\n")
subset(bias_tbl,percentile=="3rd")

cat("Comparing the bias between the standard and corrected Hadlock chart at the 10th percentile:", "\n")
subset(bias_tbl,percentile=="10th")

cat("Comparing the bias between the standard and corrected Hadlock chart at the 50th percentile:", "\n")
subset(bias_tbl,percentile=="50th")

cat("Comparing the bias between the standard and corrected Hadlock chart at the 90th percentile:", "\n")
subset(bias_tbl,percentile=="90th")

cat("Comparing the bias between the standard and corrected Hadlock chart at the 97th percentile:", "\n")
subset(bias_tbl,percentile=="97th")

# FDA analysis for one percentile
fda_bias_analysis <- function(Unc_df, Corr_df, percentile = "10th",
                              sd_log = 0.12, n = 200, B = 1000,
                              nbasis = 15, norder = 4, seed = 1){
  # returns a list with smoothed bias function, bootstrap CIs for bias and derivative
  set.seed(seed)
  # extract WGA and vectors
  WGAs <- Unc_df$WGA
  x <- WGAs
  unc <- as.numeric(Unc_df[[percentile]])
  # analytic corrected curve from corrected df
  corr_analytic <- as.numeric(Corr_df[[percentile]])
  bias_point <-  corr_analytic-unc
  
  # Bootstrap corrected curves (simulate corrected estimates at each WGA B times)
  Bmat <- matrix(NA, nrow = length(x), ncol = B)
  for(i in seq_along(x)){
    pnum <- as.numeric(gsub("th|rd|st|nd","", percentile))
    hb <- hadlock_parametric_boot(WGA = x[i], percentile = pnum,
                                  sd_log = sd_log, n = n, B = B, seed = NULL)
    Bmat[i, ] <- hb$boot
  }
  # B bootstrap bias curves (unc fixed - corrected_boot)
  bias_boot_mat <- matrix(NA, nrow = length(x), ncol = B)
  for(b in 1:B){
    bias_boot_mat[, b] <-  Bmat[, b]-unc
  }
  # smoothing basis
  rng <- range(x)
  nbasis <- nbasis
  basis <- create.bspline.basis(rangeval = rng, nbasis = nbasis, norder = norder)
  # smooth each bootstrap bias curve and also the analytic bias
  fd_boot <- vector("list", B)
  fd_vals <- matrix(NA, nrow = length(x), ncol = B)
  fd_deriv_vals <- matrix(NA, nrow = length(x), ncol = B)
  # create fine grid for evaluation (we'll evaluate at the original x)
  for(b in 1:B){
    yb <- bias_boot_mat[, b]
    # smooth with least squares smoothing onto basis using fdPar
    fdParobj <- fdPar(basis, Lfdobj = int2Lfd(2), lambda = 1e-6)
    smoothres <- smooth.basis(x, yb, fdParobj)
    fdobj <- smoothres$fd
    eval_y <- eval.fd(x, fdobj)
    eval_d1 <- eval.fd(x, fdobj, 1)
    fd_vals[, b] <- eval_y
    fd_deriv_vals[, b] <- eval_d1
  }
  # analytic smoothed bias (point)
  fdParobj <- fdPar(basis, Lfdobj = int2Lfd(2), lambda = 1e-6)
  smoothres_point <- smooth.basis(x, bias_point, fdParobj)
  fd_point <- smoothres_point$fd
  point_vals <- as.numeric(eval.fd(x, fd_point))
  point_d1 <- as.numeric(eval.fd(x, fd_point, 1))
  # compute bootstrap CIs (pointwise)
  bias_CI_lower <- apply(fd_vals, 1, quantile, probs = 0.025)
  bias_CI_upper <- apply(fd_vals, 1, quantile, probs = 0.975)
  deriv_CI_lower <- apply(fd_deriv_vals, 1, quantile, probs = 0.025)
  deriv_CI_upper <- apply(fd_deriv_vals, 1, quantile, probs = 0.975)
  out <- list(
    WGA = x,
    bias_point = point_vals,
    bias_CI_lower = bias_CI_lower,
    bias_CI_upper = bias_CI_upper,
    deriv_point = point_d1,
    deriv_CI_lower = deriv_CI_lower,
    deriv_CI_upper = deriv_CI_upper,
    fd_point = fd_point,
    fd_boot_vals = fd_vals,
    fd_boot_deriv_vals = fd_deriv_vals
  )
  return(out)
}


# FDA across all percentiles
percentile_levels<- c("3rd", "10th", "50th", "90th", "97th")
fda_res<- list()
for(p in seq_along(percentile_levels)){
    fda_res[[p]] <- fda_bias_analysis(Uncorrected_table2_df, Corrected_table2_df, percentile = percentile_levels[p],
                           sd_log = 0.12, n = 200, B = 1000, seed=1)
}

fda_df_plot <- NULL
for(p in seq_along(percentile_levels)){
    fda_df_plot[[p]] <- data.frame(WGA = fda_res[[p]]$WGA,
                       bias =fda_res[[p]]$bias_point,
                       lower = fda_res[[p]]$bias_CI_lower,
                      upper = fda_res[[p]]$bias_CI_upper,
                      percentile= percentile_levels[p])
    }

# Combined all percentile results
fda_df_plot_combined<- do.call("rbind",  fda_df_plot)

head(fda_df_plot_combined, n=10)

levels(as.factor(fda_df_plot_combined$percentile))

fda_df_plot_combined$percentile<- factor(fda_df_plot_combined$percentile, levels=percentile_levels)
levels(fda_df_plot_combined$percentile)

# Compute the overall bias
overall_bias <- fda_df_plot_combined %>%
  dplyr::group_by(WGA) %>%
  dplyr::summarize(
    median_bias = median(bias),
    lower_overall = median(lower),
    upper_overall = median(upper)
  )


overall_bias_facet <- overall_bias %>%
  mutate(percentile = "Overall (median)")  # This will allow it to appear in every facet


plot_data <- fda_df_plot_combined %>%
  bind_rows(overall_bias_facet %>%
              rename(bias = median_bias, lower = lower_overall, upper = upper_overall))

print(percentile_levels)
names(fda_res) <- percentile_levels


levels(as.factor(plot_data$percentile))

plot_data$percentile<- factor(plot_data$percentile, levels=c(percentile_levels, "Overall (median)"))

levels(as.factor(plot_data$percentile))

tail(plot_data)
percentile_bias_split<- split(plot_data, plot_data$percentile)


library(refund)

# Construct bias matrix: rows = curves (percentiles), columns = WGA points
bias_matrix <- rbind(
  percentile_bias_split$`3rd`[["bias"]],
  percentile_bias_split$`10th`[["bias"]],
  percentile_bias_split$`50th`[["bias"]],
  percentile_bias_split$`90th`[["bias"]],
  percentile_bias_split$`97th`[["bias"]],
  percentile_bias_split$`Overall (median)`[["bias"]]
)
# Now each row = a curve (percentile), columns = WGA

#  Covariates: one row per curve
covariates <- data.frame(
  percentile = factor(c("3rd","10th","50th","90th","97th","Overall (median)"))
)

covariates$percentile<- factor(covariates$percentile, levels=c("3rd","10th","50th","90th","97th","Overall (median)"))

# yind = WGA grid
yind <- 10:40  # length must match number of columns in bias_matrix

# Penalised flexible functional regression
# Implements additive regression for functional and scalar covariates and functional responses.
# Fit function-on-scalar regression
fit <- pffr(
  bias_matrix ~ percentile,
  yind = yind,
  data = covariates,
  algorithm = "gam"
)

# Summary and visualisation
summary(fit)
# plot(fit)

ggplot(plot_data, aes(x = WGA, y = bias)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(percentile)), alpha = 0.2) +
  geom_line(aes(color = factor(percentile)), size = .8) +
  labs(
    title = "",
    x = "Gestational Age (weeks)",
    y = "Functional Bias (Corrected - Standard Hadlock percentiles)",
    color = "Percentile:",
    fill = "Percentile:"
  ) +
  facet_wrap(~ percentile) + 
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )+
scale_colour_manual(values = percentile_colors,
                    labels = c("3rd", "10th", "50th", "90th", "97th","Overall (median)")) +
scale_fill_manual(values = percentile_fill,
                  labels = c("3rd", "10th", "50th", "90th", "97th","Overall (median)"))+
guides(
  colour = guide_legend(nrow = 1, byrow = TRUE),
  fill = guide_legend(nrow = 1, byrow = TRUE)
)
######################## FIGURE 2 ##############################

# percentile_colors <- c(
  # "3rd"   = "#1f78b4",  # deep blue
  # "10th"  = "#33a02c",  # green
  # "50th"  = "#ff7f00",  # orange
  # "90th"  = "#e31a1c",  # red
  # "97th"  = "#6a3d9a",   # purple
  # "Overall (median)"="black"
# )

percentile_fill <- rep("black", 6) # same for ribbons
percentile_colors<-percentile_fill


# Extract summary table from pffr
fit_sum <- summary(fit)$s.table
term_names <- rownames(fit_sum)

# Build results data frame
results_df <- data.frame(
  percentile = gsub("percentile|\\(yind\\)", "", term_names),
  edf = fit_sum[, "edf"],
  p_value = fit_sum[, "p-value"],
  stringsAsFactors = FALSE
)

# Replace empty percentile label (Intercept)
results_df$percentile[results_df$percentile == ""] <- "Intercept"

# Keep only percentiles used in plot
results_df <- subset(results_df, percentile %in% unique(plot_data$percentile))


# Create simplified label (edf + p-value only)
results_df$label_detailed <- paste0(
    ifelse(
  results_df$edf < 0.1, "Flat bias curve",
  ifelse(results_df$edf < 3, "Moderately nonlinear bias curve",
         "Highly nonlinear bias")
), "\n",
  "(edf = ", format(results_df$edf,digits=2),")" ,"\n",
  ifelse(
    results_df$p_value < 0.001, "p < 0.001***",
    ifelse(results_df$p_value < 0.05, "p < 0.05",
           paste0("p = ", sprintf("%.3f", results_df$p_value)))
  )
)




# Combine with the existing label_detailed
results_df$label_detailed <- paste0(
  results_df$label_detailed, "\n", results_df$bias_interpretation
)


# Merge labels with plotting data for geom_text
# Compute facet-specific label positions (top-right of each panel)
# Compute a per-facet label position (top centre)
plot_labels <- plot_data %>%
  group_by(percentile) %>%
  dplyr::summarise(
    x = mean(range(WGA, na.rm = TRUE)),  # mid x-position
    y = 150 # slightly above top
  ) %>%
  left_join(results_df[, c("percentile", "label_detailed")], by = "percentile")

plot_data$percentile<- factor(plot_data$percentile, levels=c(percentile_levels, "Overall (median)"))
plot_labels$percentile<- factor(plot_labels$percentile, levels=c(percentile_levels, "Overall (median)"))


# Create a data frame for the legend text
# Get all facet levels
facet_levels <- levels(plot_data$percentile)

# Create legend data for all facets
signif_legend <- data.frame(
  percentile = factor(facet_levels, levels = facet_levels),
  x = min(plot_data$WGA) + 1,    # horizontal position
  y = max(plot_data$bias, na.rm = TRUE) + 20,  # vertical position
  label = "p < 0.05 → Bias significant from 0\np ≥ 0.05 → Bias insignificant from 0"
)


ggplot(plot_data, aes(x = WGA, y = bias)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(percentile)), alpha = 0.2) +
  geom_line(aes(color = factor(percentile)), size = .8) +
  geom_text(
    data = plot_labels,
    aes(x = x, y = y, label = label_detailed, color = percentile),
    hjust = 0.55,
    vjust = 0,
    size = 3.1,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    data = signif_legend,
    aes(x = x, y = y, label = label),
    hjust = 0.005, vjust = 10,
    size = 2.8,
    fontface = "italic",
    color = "black",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Gestational Age (weeks)",
    y = "Functional Bias (Corrected − Standard Hadlock percentiles)",
    color = "Percentile:",
    fill = "Percentile:"
  ) +
  facet_wrap(~ percentile) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)  # add frame
  ) +
  scale_colour_manual(values = percentile_colors,
                      labels = names(percentile_colors)) +
  scale_fill_manual(values = percentile_fill,
                    labels = names(percentile_fill)) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    fill = guide_legend(nrow = 1, byrow = TRUE)
  )


# For derivative plot

fda_df_plot.fderiv <- NULL
for(p in seq_along(percentile_levels)){
    fda_df_plot.fderiv[[p]] <- data.frame(WGA = fda_res[[p]]$WGA,
                       deriv = fda_res[[p]]$deriv_point,
                        lower = fda_res[[p]]$deriv_CI_lower,
                       upper = fda_res[[p]]$deriv_CI_upper,
                      percentile= percentile_levels[p])
    }

# Combined all percentile results
fda_df_plot_deriv_combined<- do.call("rbind", fda_df_plot.fderiv)

fda_df_plot_deriv_combined$percentile<- factor(fda_df_plot_deriv_combined$percentile, levels=percentile_levels)
# levels(fda_df_plot_deriv_combined$percentile)

head(fda_df_plot_deriv_combined, n=10)


# Compute the overall bias
overall_d.bias <- fda_df_plot_deriv_combined %>%
  dplyr::group_by(WGA) %>%
  dplyr::summarize(
    median_d.bias = median(deriv),
    lower_overall = median(lower),
    upper_overall = median(upper)
  )


overall_d.bias_facet <- overall_d.bias %>%
  mutate(percentile = "Overall (median)")  # This will allow it to appear in every facet


plot_d.data <- fda_df_plot_deriv_combined %>%
  bind_rows(overall_d.bias_facet %>%
              rename(deriv = median_d.bias, lower = lower_overall, upper = upper_overall))
######################## FIGURE 3 ##############################

percentile_colors <- c(
  "3rd"   = "#1f78b4",  # deep blue
  "10th"  = "#33a02c",  # green
  "50th"  = "#ff7f00",  # orange
  "90th"  = "#e31a1c",  # red
  "97th"  = "#6a3d9a",   # purple
  "Overall (median)"="black"
)

percentile_fill <- percentile_colors  # same for ribbons

 ggplot(plot_d.data, aes(x = WGA, y = deriv)) +
   geom_line(aes(color = factor(percentile)), size = .8) +
   geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(percentile)), alpha = 0.2) +
  labs(title = "",
        y = "Functional rate of change of bias (grams/week)", x = "WGA",  color = "Percentile:",
    fill = "Percentile:")+
facet_wrap(~ percentile) + 
 theme_minimal()+
  facet_wrap(~ percentile) + 
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )+
scale_colour_manual(values = percentile_colors,
                    labels = c("3rd", "10th", "50th", "90th", "97th","Overall (median)")) +
scale_fill_manual(values = percentile_fill,
                  labels = c("3rd", "10th", "50th", "90th", "97th","Overall (median)"))+
guides(
  colour = guide_legend(nrow = 1, byrow = TRUE),
  fill = guide_legend(nrow = 1, byrow = TRUE)
)





# Uncorrected/standard percentiles
head(Uncorrected_table2_df)

# Corrected percentiles
head(Corrected_table2_df)

# Load libraries
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(readr)


# ---------------------------
# Loading the corrected Hadlock data
# ---------------------------
 corrected_Hadlock_df <- Corrected_table2_df
# names(corrected_Hadlock_df) <- c("WGA", "3rd", "10th", "50th", "90th", "97th")

# ---------------------------
# Define the published Hadlock reference (Uncorrected)
# ---------------------------
standard_Hadlock_df<- Uncorrected_table2_df 

standard_Hadlock_df <- standard_Hadlock_df%>%
  rename(
    p3  = `3rd`,
    p10 = `10th`,
    p50 = `50th`,
    p90 = `90th`,
    p97 = `97th`
  )


# ---------------------------
#  1. Corrected model: mean function + sd
# ---------------------------
# w=WGA
meanEFW <- function(w) {
  exp(0.578 + 0.332 * w - 0.00354 * w^2)
}

sdlog <- 0.12


mapped <- standard_Hadlock_df %>%
  mutate(
    mean_corr = meanEFW(WGA),
    corr_p_for_hadlock_3rd  = plnorm(p3,  meanlog = log(mean_corr), sdlog = sdlog),
    corr_p_for_hadlock_10th = plnorm(p10, meanlog = log(mean_corr), sdlog = sdlog),
    corr_p_for_hadlock_50th = plnorm(p50, meanlog = log(mean_corr), sdlog = sdlog),
    corr_p_for_hadlock_90th = plnorm(p90, meanlog = log(mean_corr), sdlog = sdlog),
    corr_p_for_hadlock_97th = plnorm(p97, meanlog = log(mean_corr), sdlog = sdlog)
  ) %>%
  # Estimated Percentile levels in proportions
  mutate(across(starts_with("corr_p_for_"), ~ . * 100)) %>%
  # Add rounded (whole-number) percentile columns
  mutate(
  across(
    starts_with("corr_p_for_"),
    list(rounded = ~ formatC(round(.), digits = 3, format = "fg")),
    .names = "{.col}_rounded"
  )
)







# Option II: Adjusting the standard Hadlock growth chart based on the correct critical point values

nondiscrepany_matrix<- corrected_Hadlock_df[, -1]==standard_Hadlock_df[, -1]

percentile_levels<- c("3rd", "10th", "50th", "90th", "97th")
percentile_levels_numeric<- c(3, 10, 50, 90, 97)

rounded_percentile_labels<- c('corr_p_for_hadlock_3rd_rounded','corr_p_for_hadlock_10th_rounded',
                              'corr_p_for_hadlock_50th_rounded','corr_p_for_hadlock_90th_rounded',
                              'corr_p_for_hadlock_97th_rounded')
mapped_rounded_percentile <- mapped[, c("WGA",rounded_percentile_labels)]

# Convert all columns to numeric
mapped_rounded_percentile[] <- lapply(mapped_rounded_percentile, as.numeric)

for(i in 1:dim(nondiscrepany_matrix)[1]){
    for(p in seq_along(percentile_levels)){
        if(nondiscrepany_matrix[i, p]==TRUE){
            # in case of truncation or rounding-off errors
          if(mapped_rounded_percentile[i, p+1]!=percentile_levels_numeric[p]){
          mapped_rounded_percentile[i, p+1]<- percentile_levels_numeric[p]
              
                 }
            }

        }

 }




# clean up for a publication table
# mapped_table <- mapped %>%
  # dplyr::select(WGA, ends_with("_rounded")) %>%
  # rename_with(~ gsub("corr_p_for_hadlock_|_rounded", "", .x))

mapped_table<- mapped_rounded_percentile

names(mapped_table)[1:6]<- c("Gestational age (weeks)","Corrected percentile corresponding to Hadlock 3rd",
    "Corrected percentile corresponding to Hadlock 10th", "Corrected percentile corresponding to Hadlock 50th",
        "Corrected percentile corresponding to Hadlock 90th", "Corrected percentile corresponding to Hadlock 97th")


# Corrected estimated percentile levels
mapped_table



write.csv(mapped_table, "Corrected.est.percentile.levels_Hadlock1991chart.csv")

# Sensitivity analysis (To check if critical point values at the reverse corrected estimated percentile levels
# corresponds to their corresponding standard Hadlock percentiles at say WGA=10 
cat("Standard/uncorrected 3rd or 0.03 percentile critical value at WGA=10:", "\n")
subset(mapped, WGA=="10")[["p3"]]

corrected_p<- subset(mapped, WGA=="10")$corr_p_for_hadlock_3rd
cat("Corrected percent level=", round(corrected_p, 0),"\n")

cat("Corresponding percentile critical value at this corrected percent level at WGA=10:", corrected_p,"", "gives:", "\n")
print(hadlock_percentile(wga=10, p = subset(mapped, WGA=="10")$corr_p_for_hadlock_3rd, sd_log = 0.12))



# Sensitivity analysis (To check if critical point values at the reverse corrected estimated percentile levels
# corresponds to their corresponding standard Hadlock percentiles at say WGA=24 
cat("Standard/uncorrected 3rd or 0.03 percentile critical value at WGA=24:", "\n")
subset(mapped, WGA=="24")[["p3"]]

corrected_p<- subset(mapped, WGA=="24")$corr_p_for_hadlock_3rd
cat("Corrected percent level=", round(corrected_p, 0),"\n")

cat("Corresponding percentile critical value at this corrected percent level at WGA=24:", corrected_p,"", "gives:", "\n")
print(hadlock_percentile(wga=24, p = subset(mapped, WGA=="24")$corr_p_for_hadlock_3rd, sd_log = 0.12))



# cat("Corrected Hadlock chart:", "\n")
# head(corrected_Hadlock_df, n=10)

# cat("Uncorrected/standard Hadlock chart:", "\n")
# head(standard_Hadlock_df, n=10)

# The true corrected Hadlock chart/optimal model percentile levels are perfectly aligned with the theoretical model.

# ---------------------------
# 1. Define corrected model parameters
# ---------------------------
meanEFW <- function(w) {
  exp(0.578 + 0.332 * w - 0.00354 * w^2)
}
sdlog <- 0.12

# ---------------------------
# 2. Generate corrected EFWs that exactly correspond to each percentile
# ---------------------------
# Define the function to compute qlnorm-based EFWs
qEFW <- function(p, w) {
  qlnorm(p, meanlog = log(meanEFW(w)), sdlog = sdlog)
}

# Build the corrected Hadlock table (Option 2)
corrected_exact_Hadlock_df <- tibble(
  WGA = main_corrected_Hadlock_df$WGA,
  p3  = qEFW(0.03, WGA),
  p10 = qEFW(0.10, WGA),
  p50 = qEFW(0.50, WGA),
  p90 = qEFW(0.90, WGA),
  p97 = qEFW(0.97, WGA)
)

# ---------------------------
# 3. (Optional) Verify that plnorm() of these values returns the target percentiles
# ---------------------------
validation <- corrected_exact_Hadlock_df %>%
  mutate(
    mean_corr = meanEFW(WGA),
    cdf_3rd  = plnorm(p3,  meanlog = log(mean_corr), sdlog = sdlog),
    cdf_10th = plnorm(p10, meanlog = log(mean_corr), sdlog = sdlog),
    cdf_50th = plnorm(p50, meanlog = log(mean_corr), sdlog = sdlog),
    cdf_90th = plnorm(p90, meanlog = log(mean_corr), sdlog = sdlog),
    cdf_97th = plnorm(p97, meanlog = log(mean_corr), sdlog = sdlog)
  )

# ---------------------------
# 4. Check validation (should be exactly 0.03, 0.10, 0.50, 0.90, 0.97)
# ---------------------------
# summary(dplyr::select(validation, starts_with("cdf_")))

corrected_exact_Hadlock_df <- tibble(
  WGA = main_corrected_Hadlock_df$WGA,
  p3  = qEFW(0.03, WGA),
  p10 = qEFW(0.10, WGA),
  p50 = qEFW(0.50, WGA),
  p90 = qEFW(0.90, WGA),
  p97 = qEFW(0.97, WGA)
) %>%
  mutate(across(-WGA, ~ round(., 0))) %>%  # round to nearest gram
  rename(
    `3rd`  = p3,
    `10th` = p10,
    `50th` = p50,
    `90th` = p90,
    `97th` = p97
  )

# corrected_exact_Hadlock_df

corrected_exact_p_Hadlock_df <- validation  %>%
   mutate(across(starts_with("cdf_"), ~ . * 100)) %>%
  # Add rounded (whole-number) percentile columns
  mutate(
  across(
    starts_with("cdf_"),
    list(rounded = ~ formatC(round(.), digits = 3, format = "fg")),
    .names = "{.col}_rounded"
  )
)


cat("The true corrected Hadlock chart/optimal model percentile levels are perfectly aligned with the theoretical model:
", "\n")
corrected_exact_p_Hadlock_df[, c("WGA",'cdf_3rd_rounded','cdf_10th_rounded','cdf_50th_rounded',
                                'cdf_90th_rounded','cdf_97th_rounded')]



head(standard_Hadlock_df)

head(corrected_Hadlock_df)

# Option II: Adjusting the standard Hadlock growth chart based on the correct critical point values

OptionII_Adjusted.standard_Hadlock_df<- standard_Hadlock_df

percentile_levels<- c("3rd", "10th", "50th", "90th", "97th")
percentile_levels_numeric<- c(3, 10, 50, 90, 97)

rounded_percentile_labels<- c('corr_p_for_hadlock_3rd_rounded','corr_p_for_hadlock_10th_rounded',
                              'corr_p_for_hadlock_50th_rounded','corr_p_for_hadlock_90th_rounded',
                              'corr_p_for_hadlock_97th_rounded')
mapped_rounded_percentile <- mapped[, c("WGA",rounded_percentile_labels)]

# Convert all columns to numeric
mapped_rounded_percentile[] <- lapply(mapped_rounded_percentile, as.numeric)

for(i in 1:dim(mapped_rounded_percentile)[1]){
    for(p in seq_along(percentile_levels)){
      

         # if there is a discrepancy in percentile levels, replace with the correct critical point value
        if(mapped_rounded_percentile[i, p+1]!=percentile_levels_numeric[p]){
           
              OptionII_Adjusted.standard_Hadlock_df[i, p+1]<- corrected_Hadlock_df[i, p+1]
              
                 
            }

        }

 }




cat("Option II: Adjusted Standard Hadlock Chart:")

names(OptionII_Adjusted.standard_Hadlock_df)[2:6]<- percentile_levels

write.csv(OptionII_Adjusted.standard_Hadlock_df, "OptionII_Adjusted.standard_Hadlock_df.csv")

OptionII_Adjusted.standard_Hadlock_df

####Setting the plot size and resolution (300 dpi)
options(repr.plot.width=8, repr.plot.height=8,repr.plot.res = 300) #Setting plot size

setwd("/Users/clementtwumasi/Desktop/Research Collaborations_Theo/Paper_HadlockChart_Evaluation_2025")

library("haven")

# A function to compute the  mean EFW based on the optimal model by Hadlock (1991)
hadlock_mean_efw <- function(wga){
  # wga: numeric vector of week gestational ages (WGA)
  mu_log <- 0.578 + 0.332 * wga - 0.00354 * (wga^2)
  mean_efw <- exp(mu_log)
  return(mean_efw)
}

# Testing the function at WGA=24 weeks for example 
hadlock_mean_efw(wga=24)

# Create Table 2 dataframe from Hadlock 1991 critical point values
Uncorrected_table2_df<- data.frame(
  WGA = 10:40,
  `3rd` = c(26,34,43,55,70,88,110,136,167,205,248,299,359,426,503,
            589,685,791,908,1034,1169,1313,1465,1622,1783,1946,2110,
            2271,2427,2576,2714),
  `10th` = c(29,37,48,61,77,97,121,150,185,227,275,331,398,471,556,
             652,758,876,1004,1145,1294,1453,1621,1794,1973,2154,2335,
             2513,2686,2851,3004),
  `50th` = c(35,45,58,73,93,117,146,181,223,273,331,399,478,568,670,
             785,913,1055,1210,1379,1559,1751,1953,2162,2377,2595,2813,
             3028,3236,3435,3619),
  `90th` = c(41,53,68,85,109,137,171,212,261,319,387,467,559,665,784,
             918,1068,1234,1416,1613,1824,2049,2285,2530,2781,3036,3291,
             3543,3786,4019,4234),
  `97th` = c(44,56,73,91,116,146,183,226,279,341,414,499,598,710,838,
             981,1141,1319,1513,1724,1949,2189,2441,2703,2971,3244,3516,
             3785,4045,4294,4524)
)

names(Uncorrected_table2_df)[2:6]<- c("3rd", "10th", "50th", "90th", "97th")

standard_Hadlock_chart<- Uncorrected_table2_df
cat("Standard Hadlock Chart:", "\n")
head(Uncorrected_table2_df, n=6)
tail(Uncorrected_table2_df, n=6)

clinical_data<- read_sav("Working clean dataset.sav")

dim(clinical_data)

# print(names(clinical_data))

variables_of_interest<- c("week", "birthweight")
clinical_data_main<- clinical_data[, variables_of_interest]

names(clinical_data_main)<- c("WGA", "true_birthweight")
dim(clinical_data_main)

head(clinical_data_main)

table(clinical_data_main$WGA)

# Extracting WGA of interest up to week 40 for direct comparison with the Hadlock chart
clinical_data_subset<- subset(clinical_data_main, WGA<=40)
dim(clinical_data_subset)
table(clinical_data_subset$WGA)

standard_Hadlock_chart_subset<- subset(standard_Hadlock_chart, WGA>=37 & WGA<=40)

standard_Hadlock_chart_subset

# Subsetting the data at each available observed WGA
Data_WGA_subset_available<- NULL
Data_WGA_subset_available[["WGA==37"]]<- subset(clinical_data_subset, WGA==37)
Data_WGA_subset_available[["WGA==38"]]<- subset(clinical_data_subset, WGA==38)
Data_WGA_subset_available[["WGA==39"]]<- subset(clinical_data_subset, WGA==39)
Data_WGA_subset_available[["WGA==40"]]<- subset(clinical_data_subset, WGA==40)

proportion_below_Hadlock10thCentile<- NULL

# proportion of the observed birth weight that falls below the Hadlock 10th centile
proportion_below_Hadlock10thCentile[["WGA==37"]]<- mean(Data_WGA_subset_available[["WGA==37"]]$true_birthweight
                                < standard_Hadlock_chart_subset[standard_Hadlock_chart_subset$WGA==37, "10th"],na.rm=TRUE)


proportion_below_Hadlock10thCentile[["WGA==38"]]<- mean(Data_WGA_subset_available[["WGA==38"]]$true_birthweight
                                < standard_Hadlock_chart_subset[standard_Hadlock_chart_subset$WGA==38, "10th"],na.rm=TRUE)

proportion_below_Hadlock10thCentile[["WGA==39"]]<- mean(Data_WGA_subset_available[["WGA==39"]]$true_birthweight
                                 < standard_Hadlock_chart_subset[standard_Hadlock_chart_subset$WGA==39, "10th"],na.rm=TRUE)

proportion_below_Hadlock10thCentile[["WGA==40"]]<- mean(Data_WGA_subset_available[["WGA==40"]]$true_birthweight
                                    < standard_Hadlock_chart_subset[standard_Hadlock_chart_subset$WGA==40, "10th"],na.rm=TRUE)

proportion_below_Hadlock10thCentile

# 1) Compute 95% CIs using the MLE standard error

# Extract proportions into a numeric vector
WGA <- c(37, 38, 39, 40)
p <- c(
  proportion_below_Hadlock10thCentile[["WGA==37"]],
  proportion_below_Hadlock10thCentile[["WGA==38"]],
  proportion_below_Hadlock10thCentile[["WGA==39"]],
  proportion_below_Hadlock10thCentile[["WGA==40"]]
)

# Sample sizes
n <- c(
  sum(!is.na(Data_WGA_subset_available[["WGA==37"]]$true_birthweight)),
  sum(!is.na(Data_WGA_subset_available[["WGA==38"]]$true_birthweight)),
  sum(!is.na(Data_WGA_subset_available[["WGA==39"]]$true_birthweight)),
  sum(!is.na(Data_WGA_subset_available[["WGA==40"]]$true_birthweight))
)

# Standard errors (MLE-based)
SE <- sqrt(p * (1 - p) / n)

# 95% confidence intervals
lower_CI <- p - 1.96 * SE
upper_CI <- p + 1.96 * SE

# Combine into a results data frame
results <- data.frame(
  WGA = WGA,
  proportion = p,
  n = n,
  SE = SE,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

print(results)


library(ggplot2)

ggplot(results, aes(x = WGA, y = proportion)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(
    label = sprintf("p = %.3f", proportion)
  ),
  nudge_x =0.2,    # moves text slightly to the right of the point
  size = 3
  ) +
  scale_x_continuous(breaks = results$WGA) +
  labs(
    x = "Observed WGA (weeks)",
    y = "Proportion of observed birth weight below Hadlock 10th centile (p)"
  ) +
  theme_minimal()+
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "top"
  )+
geom_hline(yintercept = 0.1, linetype = "dashed", color = "red", linewidth = 1)+
annotate("text", x = median(results$WGA), y = 0.1,
           label = "10th centile level (p = 0.1)",
           vjust = -0.5, color = "red")


# Estimating the prediction error (based on MAPE) between the predicted EFW from the model and the actual birth weights across all 4 observed WGA

# Combine all WGA data into one vector
observed_bw <- c(
  Data_WGA_subset_available[["WGA==37"]]$true_birthweight,
  Data_WGA_subset_available[["WGA==38"]]$true_birthweight,
  Data_WGA_subset_available[["WGA==39"]]$true_birthweight,
  Data_WGA_subset_available[["WGA==40"]]$true_birthweight
)

# Create matching WGA vector
wga_all <- c(
  rep(37, length(Data_WGA_subset_available[["WGA==37"]]$true_birthweight)),
  rep(38, length(Data_WGA_subset_available[["WGA==38"]]$true_birthweight)),
  rep(39, length(Data_WGA_subset_available[["WGA==39"]]$true_birthweight)),
  rep(40, length(Data_WGA_subset_available[["WGA==40"]]$true_birthweight))
)


predicted_efw <- hadlock_mean_efw(wga_all)


# Remove missing values safely
valid_idx <- !is.na(observed_bw) & !is.na(predicted_efw) & observed_bw > 0

obs <- observed_bw[valid_idx]
pred <- predicted_efw[valid_idx]

# Compute MAPE (%)
MAPE <- mean(abs((obs - pred) / obs)) * 100

# Print result
MAPE


mape_by_wga <- tapply(
  abs((observed_bw - predicted_efw) / observed_bw) * 100,
  wga_all,
  mean,
  na.rm = TRUE
)

mape_by_wga


observed_bw <- c(
  Data_WGA_subset_available[["WGA==37"]]$true_birthweight,
  Data_WGA_subset_available[["WGA==38"]]$true_birthweight,
  Data_WGA_subset_available[["WGA==39"]]$true_birthweight,
  Data_WGA_subset_available[["WGA==40"]]$true_birthweight
)

wga_all <- c(
  rep(37, length(Data_WGA_subset_available[["WGA==37"]]$true_birthweight)),
  rep(38, length(Data_WGA_subset_available[["WGA==38"]]$true_birthweight)),
  rep(39, length(Data_WGA_subset_available[["WGA==39"]]$true_birthweight)),
  rep(40, length(Data_WGA_subset_available[["WGA==40"]]$true_birthweight))
)

predicted_efw <- hadlock_mean_efw(wga_all)


valid_idx <- !is.na(observed_bw) & !is.na(predicted_efw) & observed_bw > 0


obs  <- observed_bw[valid_idx]
pred <- predicted_efw[valid_idx]
wga  <- wga_all[valid_idx]


df_mape <- data.frame(
  WGA = wga,
  Observed = obs,
  Predicted = pred
)


nrow(df_mape)

# 1) Scatter plot: Observed vs Predicted (with identity line)
library(ggplot2)

ggplot(df_mape, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    x = "Observed birth weight (g)",
    y = "Hadlock-predicted EFW (g)",
    title = "Observed vs Predicted Fetal Weight"
  ) +
  theme_minimal()

# 2) Bland–Altman plot (agreement plot)
df_mape$Mean <- (df_mape$Observed + df_mape$Predicted) / 2
df_mape$Diff <- df_mape$Observed - df_mape$Predicted

bias <- mean(df_mape$Diff)
loa_upper <- bias + 1.96 * sd(df_mape$Diff)
loa_lower <- bias - 1.96 * sd(df_mape$Diff)

ggplot(df_mape, aes(x = Mean, y = Diff)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = bias, linetype = "solid", color = "blue") +
  geom_hline(yintercept = loa_upper, linetype = "dashed", color = "red") +
  geom_hline(yintercept = loa_lower, linetype = "dashed", color = "red") +
  labs(
    x = "Mean of observed and predicted weight (g)",
    y = "Difference (Observed – Predicted, g)",
    title = "Bland–Altman Plot"
  ) +
  theme_minimal()
# 3) Distribution of percentage errors
df_mape$APE <- abs((df_mape$Observed - df_mape$Predicted) / df_mape$Observed) * 100

ggplot(df_mape, aes(x = APE)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  labs(
    x = "Absolute Percentage Error (%)",
    y = "Count",
    title = "Distribution of Absolute Percentage Errors"
  ) +
  theme_minimal()

# Make sure APE exists
df_mape$APE <- abs((df_mape$Observed - df_mape$Predicted) / df_mape$Observed) * 100

# MAPE by WGA
mape_by_wga <- aggregate(APE ~ WGA, data = df_mape, FUN = mean)

# Overall MAPE
overall_mape <- mean(df_mape$APE, na.rm = TRUE)


library(ggplot2)

ggplot(df_mape, aes(x = factor(WGA), y = APE)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  labs(
    x = "Observed WGA (weeks)",
    y = "Absolute Percentage Error (%)"
      # ,title = "MAPE by Gestational Age"
  ) +
  theme_minimal() +
  
  # Annotate WGA-specific MAPE values
  geom_text(
    data = mape_by_wga,
    aes(
      x = factor(WGA),
      y = max(df_mape$APE, na.rm = TRUE) * 0.80,
      label = sprintf("MAPE = %.2f%%", APE)
    ),
    size = 3.5
  ) +
  
  # Annotate overall MAPE
  annotate(
    "text",
    x = 2,
    y = max(df_mape$APE, na.rm = TRUE),
    label = sprintf("Overall MAPE = %.2f%%", overall_mape),
    hjust = 0,
    size = 3.5,
    fontface = "bold"
  )+
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "top"
  )





