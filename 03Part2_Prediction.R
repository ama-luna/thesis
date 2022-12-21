### Part 2: Prediction ###

## -------------------------------------------------------------------------------------------------
# Packages
library(pROC)
library(psfmi)
library(tidyverse)
library(ggnewscale)
library(ggh4x)
library(coxed) 
library(mcreplicate) 
library(parallel)
library(pbmcapply)
library(pbapply)
library(reshape2)
library(latex2exp)
library(showtext)

## -------------------------------------------------------------------------------------------------
# Performance measure calculation

fo_cog_data <- readRDS("fo_cog_data.rds")
po_tot_data <- readRDS("po_tot_data.rds")

fo_cog_data <- fo_cog_data %>% filter(!is.na(fo_cog))
po_tot_data <- po_tot_data %>% filter(!is.na(po_tot))

fo_cog_data <- fo_cog_data %>% arrange(PMC_ResearchID, TimeSinceDiag)
po_tot_data <- po_tot_data %>% arrange(PMC_ResearchID, TimeSinceDiag)

brier_score <- function(obs, pred) {
  mean((obs - pred)^2)
}

calibration_int_slope <- function(prob,y){
  
  dat <- data.frame(e = prob, o = y)
  dat$e[dat$e < 0.0000000001] = 0.0000000001
  dat$e[dat$e > 0.9999999999] = 0.9999999999
  dat$logit_func <- qlogis(dat$e)
  
  model <- glm(formula = o~logit_func, 
               family = binomial(link = "logit"), dat)
  
  slope <- as.numeric(model$coefficients[2])
  intercept <- as.numeric(model$coefficients[1])
  
  cox_results <- c(intercept, slope)
  names(cox_results) <- c("Intercept", "Slope")
  
  return(cox_results)
}

calculate_performance <- function(preds, obs){
  
  if (length(levels(as.factor(obs))) >= 2){
    brier <- mean((obs - preds)^2)
    brier_scaled <- 1 - (brier_score(obs, preds) / (mean(obs) * (1 - mean(obs))^2 + (1 - mean(obs)) * mean(obs)^2))
    auc <- as.numeric(roc(obs, preds, direction = "<", quiet = TRUE)$auc)
    slope <- calibration_int_slope(preds, obs)[2]
  } else {
    brier <- NA
    brier_scaled <- NA
    auc <- NA
    slope <- NA
  }
  
  performance_measures <- c(brier, brier_scaled, auc, slope)
  names(performance_measures) <- c("Brier", "Brier_scaled", "AUC", "Slope")
  
  return(performance_measures)
  
}

## -------------------------------------------------------------------------------------------------
# Create time window for risk scores
# FO_COG  

hist(fo_cog_data$TimeSinceDiag)

time_window_fo <- fo_cog_data %>% group_by(PMC_ResearchID) %>%
  filter(TimeSinceDiag  >= 8 & TimeSinceDiag <= 16)

time_window_fo <- time_window_fo %>% group_by(PMC_ResearchID) %>%  slice(which.min(abs(TimeSinceDiag - 12)))
hist(time_window_fo$TimeSinceDiag)

time_window_fo_p2_single_risk <- time_window_fo %>% distinct(PMC_ResearchID, .keep_all = TRUE)
time_window_fo_p2_single_risk$risk <- ifelse(time_window_fo_p2_single_risk$fo_cog <67, 0, 1)
time_window_fo_p2_single_risk$risk  <- factor(time_window_fo_p2_single_risk$risk , levels = c(0,1), labels = c("risk", "no_risk"))

table(time_window_fo_p2_single_risk$risk)

##
# PO_TOT

time_window_po <- po_tot_data %>% group_by(PMC_ResearchID) %>%
  filter(TimeSinceDiag  >= 8 & TimeSinceDiag <= 16)

time_window_po <- time_window_po %>% group_by(PMC_ResearchID) %>%  slice(which.min(abs(TimeSinceDiag - 12)))
hist(time_window_po$TimeSinceDiag)

time_window_po_p2_single_risk <- time_window_po %>% distinct(PMC_ResearchID, .keep_all = TRUE)
time_window_po_p2_single_risk$risk <- ifelse(time_window_po_p2_single_risk$po_tot < 67, 0, 1)
time_window_po_p2_single_risk$risk  <- factor(time_window_po_p2_single_risk$risk , levels = c(0,1), labels = c("risk", "no_risk"))

table(time_window_po_p2_single_risk$risk)

## -------------------------------------------------------------------------------------------------
# Create time window for delta
# FO_COG

time_window_fo <- fo_cog_data %>% group_by(PMC_ResearchID) %>%
  filter(TimeSinceDiag  >= 8 & TimeSinceDiag <= 22)

time_window_fo <- time_window_fo %>% group_by(PMC_ResearchID) %>% mutate(timediff = TimeSinceDiag - lag(TimeSinceDiag))
time_window_fo <- time_window_fo %>% group_by(PMC_ResearchID) %>% mutate(cumul_timediff = cumsum(ifelse(is.na(timediff), 0, timediff)) + timediff*0)
time_window_fo <- time_window_fo %>% group_by(PMC_ResearchID) %>%
  filter((is.na(cumul_timediff) | cumul_timediff <= 12))

time_window_fo <- time_window_fo%>% group_by(PMC_ResearchID) %>%
  filter(n() >= 2)

time_window_fo <- time_window_fo %>% group_by(PMC_ResearchID) %>% filter((row_number()==1 | row_number()==n()) & n() >= 2)
time_window_fo <- time_window_fo %>% ungroup()

delta <- time_window_fo %>% group_by(PMC_ResearchID) %>% summarise(delta = last(fo_cog) - first(fo_cog))
time_window_fo_p2 <- merge(time_window_fo, delta, by = "PMC_ResearchID")

time_window_fo_p2_single <- time_window_fo_p2 %>% distinct(PMC_ResearchID, .keep_all = TRUE)
time_window_fo_p2_single$delta_binary <- ifelse(time_window_fo_p2_single$delta < 0, 0, 1)
time_window_fo_p2_single$delta_binary <- factor(time_window_fo_p2_single$delta_binary, levels = c(0, 1), labels = c("worse", "better"))

table(time_window_fo_p2_single$delta_binary)

##
# PO_TOT

time_window_po <- po_tot_data %>% group_by(PMC_ResearchID) %>%
  filter(TimeSinceDiag  >= 8 & TimeSinceDiag <= 22)

time_window_po <- time_window_po %>% group_by(PMC_ResearchID) %>% mutate(timediff = TimeSinceDiag - lag(TimeSinceDiag))
time_window_po <- time_window_po %>% group_by(PMC_ResearchID) %>% mutate(cumul_timediff = cumsum(ifelse(is.na(timediff), 0, timediff)) + timediff*0)
time_window_po <- time_window_po %>% group_by(PMC_ResearchID) %>% 
  filter((is.na(cumul_timediff) | cumul_timediff <= 12))

time_window_po <- time_window_po%>% group_by(PMC_ResearchID) %>% 
  filter(n() >= 2)

time_window_po <- time_window_po %>% group_by(PMC_ResearchID) %>% filter((row_number()==1 | row_number()==n()) & n() >= 2)

delta <- time_window_po %>% group_by(PMC_ResearchID) %>% summarise(delta = last(po_tot) - first(po_tot))
time_window_po_p2 <- merge(time_window_po, delta, by = "PMC_ResearchID")

time_window_po_p2_single <- time_window_po_p2 %>% distinct(PMC_ResearchID, .keep_all = TRUE)
time_window_po_p2_single$delta_binary <- ifelse(time_window_po_p2_single$delta < 0, 0, 1)
time_window_po_p2_single$delta_binary <- factor(time_window_po_p2_single$delta_binary, levels = c(0, 1), labels = c("worse", "better"))

table(time_window_po_p2_single$delta_binary)

## -------------------------------------------------------------------------------------------------
# Repeated K-Fold Cross-Validation (RCV)

repeated_cv <- function(data, k, r){
  
  data_cv <- as.data.frame(data)
  
  # Rename outcome variable
  if ("risk" %in% colnames(data)){
    data_cv$y <- data_cv$risk
    data_cv$y <- ifelse(data_cv$y == "risk", 0, 1)
  } else if ("delta_binary" %in% colnames(data)){
    data_cv$y <- data_cv$delta_binary
    data_cv$y <- ifelse(data_cv$y == "worse", 0, 1)
  }
  
  #data_cv$y <- ifelse(as.numeric(data_cv$risk) ==1, 0, 1)
  groups <- unique(data_cv$y)
  num_groups <- length(groups) 
  data_cv$folds <- numeric(nrow(data_cv))
  results_brier <- matrix(nrow = k, ncol = r, dimnames = list(c(seq(1:k)), paste0(c("rep"),c(seq(1:r)))))
  results_brier_scaled <- matrix(nrow = k, ncol = r, dimnames = list(c(seq(1:k)), paste0(c("rep"),c(seq(1:r)))))
  results_auc <- matrix(nrow = k, ncol = r, dimnames = list(c(seq(1:k)), paste0(c("rep"),c(seq(1:r)))))
  results_slope <- matrix(nrow = k, ncol = r, dimnames = list(c(seq(1:k)), paste0(c("rep"),c(seq(1:r)))))
  baseline <- matrix(nrow = k, ncol = r, dimnames = list(c(seq(1:k)), paste0(c("rep"),c(seq(1:r)))))
  predictions <- matrix(nrow = nrow(data_cv), ncol = k)

  for(j in 1:r){
    
    for(g in 1:length(unique(data_cv$y))){
      data_cv$index_group <- data_cv$y==groups[g]
      n <- sum(data_cv$index_group)
      indexes <- which(data_cv$index_group==TRUE, arr.ind=TRUE)
      data_cv$folds[indexes] <- match(rep(1:k,length=n),sample(1:k))
      data_cv <- data_cv %>% select(-index_group)
    }
    
    for(i in 1:k){
      test_ind <- which(data_cv$folds==i,arr.ind=TRUE)
      train <- data_cv[-test_ind, ]
      test <- data_cv[test_ind, ]

      model_fit <- glm(y ~ overall_CT + AgeAtDiagnosis + cumul_duration + TimeSinceDiag, data = train, family = "binomial")
      fitted_results <- predict(model_fit, newdata = test, type = "response")
      results_brier[i,j] <-  calculate_performance(fitted_results, test$y)[1]
      results_brier_scaled[i,j] <- calculate_performance(fitted_results, test$y)[2]
      results_auc[i,j] <- calculate_performance(fitted_results, test$y)[3]
      results_slope[i,j] <- calculate_performance(fitted_results, test$y)[4]
      
    }
    
  }
  
  final_brier <- mean(apply(results_brier, 2, mean, na.rm = TRUE), na.rm = TRUE)
  final_brier_scaled <- mean(apply(results_brier_scaled, 2, mean, na.rm = TRUE), na.rm = TRUE)
  final_auc <- mean(apply(results_auc, 2, mean, na.rm = TRUE), na.rm = TRUE)
  final_slope <- mean(apply(results_slope, 2, mean, na.rm = TRUE), na.rm = TRUE)
  
  cv_error <- c(final_brier, final_brier_scaled, final_auc, 
                final_slope
  )
  names(cv_error) <- c("Brier", "Brier_scaled", "AUC", "Slope")
  
  return(list(cv_error))
  
}

# Confidence interval calculation
confint_validation <- function(method, data, b, k = k, r = r, f = f, boot_num = boot_num){
  
  results_boot <- pbreplicate(b, {
    
    if (method == "boot") {
      bootstrap_0632plus(f, data, boot_num)$error
    } else if (method == "cv") {
      repeated_cv(data, k, r)
    }
  })
  
  estimates <- do.call(rbind, results_boot)
  value <- apply(estimates, 2, mean)
  confints <- rbind(bca(estimates[,1]), bca(estimates[,2]), bca(estimates[,3]), bca(estimates[,4]))
  colnames(confints) <- c("2.5%", "97.5%")
  confints <- cbind(value, confints)
  confints <- as.data.frame(confints)
  confints <- rownames_to_column(confints, "variable")
  
  return(confints)
  
}

# Results RCV
fo_results_risk_cv_ci <- confint_validation("cv", time_window_fo_p2_single_risk, 1000, k = 10, r = 10)
po_results_risk_cv_ci <- confint_validation("cv", time_window_po_p2_single_risk, 1000, k = 10, r = 10)
fo_results_delta_cv_ci <- confint_validation("cv", time_window_fo_p2_single, 1000, k = 10, r = 10)
po_results_delta_cv_ci <- confint_validation("cv", time_window_po_p2_single, 1000, k = 10, r = 10)

## Plotting results
# Prepare the data

fo_risk_cv <- fo_results_risk_cv_ci
fo_risk_cv$Outcome <- "fo_risk_cv"

fo_delta_cv <- fo_results_delta_cv_ci
fo_delta_cv$Outcome <- "fo_delta_cv"

po_risk_cv <- po_results_risk_cv_ci
po_risk_cv$Outcome <- "po_risk_cv"

po_delta_cv <- po_results_delta_cv_ci
po_delta_cv$Outcome <- "po_delta_cv"

all_estimates <- as.data.frame(rbind(fo_delta_cv, po_delta_cv, fo_risk_cv, po_risk_cv))
col_order <- c( "Outcome", "variable", "value", "2.5%",
                "97.5%")
all_estimates <- all_estimates[, col_order]

# Create long data
all_estimates$Outcome <- factor(all_estimates$Outcome, levels = c("fo_delta_cv", "po_delta_cv", "fo_risk_cv", "po_risk_cv"), labels = c("Delta - CF", "Delta - QOL", "Risk - CF", "Risk - QOL"))
all_estimates$variable <- factor(all_estimates$variable, levels = c("Brier", "Brier_scaled", "AUC", "Slope"), labels = c("Brier", TeX("$Brier_{scaled}$"), "AUC", "Slope"))

names(all_estimates)[4:5] <- c("lower", "upper")
all_estimates[3:5] <- apply(all_estimates[3:5], 2, round, 3)

# Remove delta-cf for plotting
all_estimates <- all_estimates[5:nrow(all_estimates),]

#########

# 0.632+ bootstrap

bootstrap_0632plus <- function(f, data, boot_num) {
  
  # Rename outcome variable
  if ("risk" %in% colnames(data)){
    data$y <- data$risk
    data$y <- ifelse(data$y == "risk", 0, 1)
  } else if ("delta_binary" %in% colnames(data)){
    data$y <- data$delta_binary
    data$y <- ifelse(data$y == "worse", 0, 1)
  }
  
  # Fit models, get predictions for apparent training data
  model_fit <- glm(formula = f, data = data, family = "binomial")
  fitted_results <- predict(model_fit, newdata = data, type = "response")
  
  # Calculate apparent error
  boot_apparent <- calculate_performance(fitted_results, data$y)
  
  seq_y <- seq_along(data$y)
  results_brier <- rep(NA, length(data$y))
  results_brier_scaled <- rep(NA, length(data$y))
  results_auc <- rep(NA, length(data$y))
  results_cox_slope <- rep(NA, length(data$y))
  
  loo_boot_error_rates <- pblapply(seq_len(boot_num), function(b) {
    
    # Create bootstrap samples
    train <- sample(seq_y, replace = TRUE)
    test <- which(!(seq_y %in% train))
    
    # Fit models, get predictions and performance estimates
    model_fit <- glm(formula = f, data = data[train,], family = "binomial")
    fitted_results <- predict(model_fit, newdata = data[test,], type = "response")
    brier <- replace(results_brier, test,  calculate_performance(fitted_results, data[test,]$y)[1])
    b_scaled <-  replace(results_brier_scaled, test, calculate_performance(fitted_results, data[test,]$y)[2])
    auc <- replace(results_auc, test, calculate_performance(fitted_results, data[test,]$y)[3])
    cox_slope <-  replace(results_cox_slope, test, calculate_performance(fitted_results, data[test,]$y)[4])
    boot_results <- rbind(brier, b_scaled, auc, cox_slope)
  })
  
  # Get mean performance over bootstrap samples
  loo_boot_error_rates <- as.data.frame(do.call(rbind, loo_boot_error_rates))
  brier <- loo_boot_error_rates[grep("brier", rownames(loo_boot_error_rates)),]
  brier_scaled <- loo_boot_error_rates[grep("b_scaled", rownames(loo_boot_error_rates)),]
  auc <- loo_boot_error_rates[grep("auc", rownames(loo_boot_error_rates)),]
  slope <- loo_boot_error_rates[grep("cox_slope", rownames(loo_boot_error_rates)),]
  
  brier_mean <- colMeans(brier, na.rm = TRUE)
  auc_mean <- colMeans(auc, na.rm = TRUE)
  brier_scaled_mean <- colMeans(brier_scaled, na.rm = TRUE)
  slope_mean <- colMeans(slope, na.rm = TRUE)
  boot_loo <- c(mean(brier_mean, na.rm = TRUE), mean(brier_scaled_mean, na.rm = TRUE),
                mean(auc_mean, na.rm = TRUE), mean(slope_mean, na.rm = TRUE))
  names(boot_loo) <- c("Brier", "Brier_scaled", "AUC", "Slope")
  
  # Store estimates for each repeat
  brier_est <- as.vector(rowMeans(brier, na.rm = TRUE))
  brier_scaled_est <- as.vector(rowMeans(brier_scaled, na.rm = TRUE))
  auc_est <- as.vector(rowMeans(auc, na.rm = TRUE))
  slope_est <- as.vector(rowMeans(slope, na.rm = TRUE))
  estimates <- as.data.frame(cbind(brier_est, brier_scaled_est, auc_est, slope_est))
  colnames(estimates) <- c("Brier", "Brier_scaled", "AUC", "Slope")
  
  ### Calculate 0.632+ estimator
  
  # Calculate no-information rate
  obs <- data$y
  combos <- tidyr::crossing(fitted_results, obs) # Calculate performance for all combinations
  gamma_hat <- calculate_performance(combos$fitted_results, combos$obs)
  
  # Calculate overfitting rate
  R_hat <- (boot_loo - boot_apparent) / (gamma_hat - boot_apparent)
  
  # Calculate relative weight
  w_hat <- 0.632 / (1 - 0.368 * R_hat)
  
  # Calculate estimator
  boot_error <- (1 - w_hat) * boot_apparent + w_hat * boot_loo
  
  return(list(error = list(boot_error), estimates = estimates))
  
}

# FO results risk
data <- as.data.frame(time_window_fo_p2_single_risk %>% ungroup() %>% select(overall_CT, AgeAtDiagnosis, cumul_duration, TimeSinceDiag, risk))
x <- data %>% select(-risk)
f <- reformulate(colnames(x), response = "risk")

fo_risk_ci <- confint_validation("boot", data, 1000, f = f, boot_num = 100)

## PO results risk
data <- as.data.frame(time_window_po_p2_single_risk %>% ungroup() %>% select(overall_CT, AgeAtDiagnosis, cumul_duration, TimeSinceDiag, risk))
x <- data %>% select(-risk)
f <- reformulate(colnames(x), response = "risk")

po_risk_ci <- confint_validation("boot", data, 1000, f = f, boot_num = 100)

## FO results delta 
data <- as.data.frame(time_window_fo_p2_single %>% ungroup() %>% select(overall_CT, AgeAtDiagnosis, cumul_duration, TimeSinceDiag, delta_binary))
x <- data %>% select(-delta_binary)
f <- reformulate(colnames(x), response = "delta_binary")

fo_delta_ci <- confint_validation("boot", data, 1000, f = f, boot_num = 100)

## PO results delta 
data <- as.data.frame(time_window_po_p2_single %>% ungroup() %>% select(overall_CT, AgeAtDiagnosis, cumul_duration, TimeSinceDiag, delta_binary))
x <- data %>% select(-delta_binary)
f <- reformulate(colnames(x), response = "delta_binary")

po_delta_ci <- confint_validation("boot", data, 1000, f = f, boot_num = 100)

## Plotting results
# Prepare the data

fo_risk <- fo_risk_ci
fo_risk$Outcome <- "fo_risk"

fo_delta <- fo_delta_ci
fo_delta$Outcome <- "fo_delta"

po_risk <- po_risk_ci
po_risk$Outcome <- "po_risk"

po_delta <- po_delta_ci
po_delta$Outcome <- "po_delta"

all_estimates <- as.data.frame(rbind(fo_delta, po_delta, fo_risk, po_risk))
col_order <- c( "Outcome", "variable", "value", "2.5%",
                "97.5%")
all_estimates <- all_estimates[, col_order]

# Create long data
all_estimates$Outcome <- factor(all_estimates$Outcome, levels = c("fo_delta","po_delta", "fo_risk",  "po_risk"), labels = c("Delta - CF", "Delta - QOL", "Risk - CF", "Risk - QOL"))
all_estimates$variable <- factor(all_estimates$variable, levels = c("Brier", "Brier_scaled", "AUC", "Slope"), labels = c("Brier", TeX("$Brier_{scaled}$"), "AUC", "Slope"))
names(all_estimates)[4:5] <- c("lower", "upper")

all_estimates[3:5] <- apply(all_estimates[3:5], 2, round, 3)

# Remove delta-cf

all_estimates <- all_estimates[5:nrow(all_estimates),]

# Plotting
scales_y <- list(
  "Brier" = scale_y_continuous(limits = c(0.1, 0.3), breaks = seq(0.1, 0.3, 0.05)),
  "Brier[scaled]" = scale_y_continuous(limits = c(-0.2, 0.1), breaks = seq(-0.2, 0.1, 0.05)),
  "AUC" = scale_y_continuous(limits = c(0.4, 0.7), breaks = seq(0.4, 0.7, 0.05)),
  "Slope" = scale_y_continuous(limits = c(-0.1, 1.2), breaks = seq(-0.1, 1.2, 0.2))
)

plotting_results_dot <- function(font, title){
  plot1 <- ggplot(all_estimates, mapping = aes(Outcome, value, col = Outcome, fill  = Outcome)) + 
    geom_point(size = 3, stat = "identity") +
    geom_errorbar(aes(ymin=lower,ymax=upper), size = 1) +
    facet_wrap(~variable, scales="free_y", labeller = label_parsed) +
    labs(title = paste0(title, ": \n Estimates of performance measures"), x = "Outcome", y = "Value") +
    theme_minimal() +
    scale_fill_manual(values = alpha(c("#536CB1","#F78C00", "#6DBF1B", "#B266FF"), 0.5)) +
    scale_colour_manual(values = alpha(c("#536CB1","#F78C00", "#6DBF1B", "#B266FF"), 1)) +
    theme(text=element_text(family=font, face="plain", size=18),
          plot.title = element_text(hjust = .5, face = "bold", size = 22, margin = unit(c(0, 0, 10, 0), "mm")),
          plot.subtitle = element_text(hjust = .5, face = "italic", size = 24),
          legend.position = "bottom",
          legend.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5, linetype = "solid"),
          plot.margin=grid::unit(c(20,20,20,20), "mm"),
          panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                            colour = "#F5F5F5"),
          panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                            colour = "#F5F5F5"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = "lightgrey", fill = NA),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 90, vjust = 3, margin = unit(c(t = 0, r = 10, l = 0, b = 0), "mm")))
  plot2 <- plot1 +
    facetted_pos_scales(
      y = scales_y)
  
  plot2
}


