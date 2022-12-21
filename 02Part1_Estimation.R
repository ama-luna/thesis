### Part 1: Estimation ###

## -------------------------------------------------------------------------------------------------
# Packages
library(nlme)
library(tidyverse)
library(stargazer)
library(lmeresampler)
library(doParallel)
library(foreach)
library(ggeffects)
library(ggthemes)
library(lme4)
library(performance)
library(RColorBrewer)
library(gsubfn)

## -------------------------------------------------------------------------------------------------
# Assessing linear model, ignoring correlational structure of the data 

fo_cog_data <- readRDS("fo_cog_data.rds")
po_tot_data <- readRDS("po_tot_data.rds")

fo_cog_data <- fo_cog_data %>% filter(!is.na(fo_cog))
po_tot_data <- po_tot_data %>% filter(!is.na(po_tot))

fo_cog_data <- fo_cog_data %>% arrange(PMC_ResearchID, TimeSinceDiag)
po_tot_data <- po_tot_data %>% arrange(PMC_ResearchID, TimeSinceDiag)

model_lm <- lm(fo_cog ~ AgeAtDiagnosis + overall_CT  + total_anstat_cases*TimeSinceDiag,
               data = fo_cog_data, na.action = na.exclude)

plot(model_lm, which = c(1), col = 1, add.smooth = FALSE,
     caption = "")
plot(fo_cog_data$total_anstat_cases, resid(model_lm), xlab = "AnStat",
     ylab = "Residuals")
plot(fo_cog_data$TimeSinceDiag, resid(model_lm), xlab = "Time",
     ylab = "Residuals")
plot(fo_cog_data$overall_CT, resid(model_lm), xlab = "CT",
     ylab = "Residuals")
plot(fo_cog_data$AgeAtDiagnosis, resid(model_lm), xlab = "Age",
     ylab = "Residuals")

## -------------------------------------------------------------------------------------------------
# Comparison of models

model_lme <- lme(fo_cog ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 + TimeSinceDiag|PMC_ResearchID,
                 method = "REML", data = fo_cog_data, na.action = na.exclude
                 
)

model_lme2 <- lme(fo_cog ~ overall_CT*TimeSinceDiag + total_anstat_cases*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*total_anstat_cases + overall_CT*TimeSinceDiag*AgeAtDiagnosis + total_anstat_cases*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 |PMC_ResearchID,
                  method = "REML", data = fo_cog_data, na.action = na.exclude,
                  correlation = corCAR1(form = ~ 1|PMC_ResearchID),
                  weights = varExp(form = ~TimeSinceDiag)
                  
)

model_lme3 <- lme(fo_cog ~ overall_CT  + total_anstat_cases*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag, random = ~ 1 + TimeSinceDiag|PMC_ResearchID,
                  method = "REML", data = fo_cog_data, na.action = na.exclude,
                  correlation = corCAR1(form = ~ 1|PMC_ResearchID)
                  
)

model_lme4 <- lme(fo_cog ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 |PMC_ResearchID, 
                  method = "REML", data = fo_cog_data, na.action = na.exclude,
                  correlation = corCAR1(form = ~ 1|PMC_ResearchID)
                  
)

model_lme5 <- lme(fo_cog ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = list(PMC_ResearchID = pdDiag(~ TimeSinceDiag)),
                  method = "REML", data = fo_cog_data, na.action = na.exclude
                  
)

model_lme7 <- lme(fo_cog ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 |PMC_ResearchID,
                  method = "REML", data = fo_cog_data, na.action = na.exclude
                  
)

model_lme8 <- lme(fo_cog ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = list(PMC_ResearchID = pdDiag(~ TimeSinceDiag)),
                  method = "REML", data = fo_cog_data, na.action = na.exclude,
                  correlation = corExp(form = ~ 1|PMC_ResearchID),
                  weights = varExp(form = ~TimeSinceDiag)
                  
)

# Only showing one option of the several anova's that were run to compare models in terms of random structure
anova(model_lme, model_lme2)

#######

# Only selected number of models shown, same comparison was used as for fo_cog_data models directly above

model_lme <- lme(po_tot ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 + TimeSinceDiag|PMC_ResearchID,
                 method = "REML", data = po_tot_data
                 
)

model_lme2 <- lme(po_tot ~ overall_CT*TimeSinceDiag + total_anstat_cases*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*total_anstat_cases + overall_CT*TimeSinceDiag*AgeAtDiagnosis + total_anstat_cases*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 |PMC_ResearchID,
                  method = "ML", data = po_tot_data, na.action = na.exclude,
                  correlation = corCAR1(form = ~ 1|PMC_ResearchID),
                  weights = varExp(form = ~TimeSinceDiag)
                  
)

model_lme3 <- lme(po_tot ~ overall_CT + DiagnosisMainFinal + total_anstat_cases*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag, random = ~ 1 + TimeSinceDiag|PMC_ResearchID,
                  method = "REML", data = po_tot_data, na.action = na.exclude,
                  correlation = corCAR1(form = ~ 1|PMC_ResearchID)
                  
)

model_lme5 <- lme(po_tot ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = list(PMC_ResearchID = pdDiag(~ TimeSinceDiag)),
                  method = "REML", data = po_tot_data, na.action = na.exclude,
                  correlation = corCAR1(form = ~ 1|PMC_ResearchID)
                  
)

model_lme7 <- lme(po_tot ~ overall_CT*TimeSinceDiag + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + overall_CT*TimeSinceDiag*cumul_duration + overall_CT*TimeSinceDiag*AgeAtDiagnosis + cumul_duration*TimeSinceDiag*AgeAtDiagnosis, random = ~ 1 |PMC_ResearchID,
                  method = "REML", data = po_tot_data, na.action = na.exclude
                  
)

## -------------------------------------------------------------------------------------------------

# lme4-package comparison for random slopes
fm0 <- lmer(po_tot ~ overall_CT*TimeSinceDiag  + total_anstat_cases*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + (1+TimeSinceDiag|PMC_ResearchID), data=po_tot_data)
fm0 <- lmer(po_tot ~ overall_CT*TimeSinceDiag  + cumul_duration*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + (1+TimeSinceDiag|PMC_ResearchID), data=po_tot_data)

fm1 <- update(fm0, .~. - TimeSinceDiag)

# Check convergence and singularity for slopes model
check_convergence(fm0) 
check_singularity(fm0) 

# rePCA in lme4-package for slopes model
summary(rePCA(fm0))

fm0 <- lmer(fo_cog ~ overall_CT + total_anstat_cases* + AgeAtDiagnosis + TimeSinceDiag + (1+TimeSinceDiag|PMC_ResearchID), data=fo_cog_data)
fm0 <- lmer(fo_cog ~ overall_CT + cumul_duration + AgeAtDiagnosis + TimeSinceDiag + (1+TimeSinceDiag|PMC_ResearchID), data=fo_cog_data)

fm1 <- update(fm0, .~. - TimeSinceDiag)

# Check convergence and singularity for slopes model
check_convergence(fm0) 
check_singularity(fm0) 

# rePCA in lme4-package for slopes model
summary(rePCA(fm0))

# Example of method by Stram and Lee (1994), see thesis for details, using mixture of reference distributions
0.5 * ((1 - pchisq(4.360107, 1)) + (1 - pchisq(4.360107, 2)))
0.5 * ((1 - pchisq(4.246151, 1)) + (1 - pchisq(4.246151, 2)))

## -------------------------------------------------------------------------------------------------
# Check assumptions, only one model shown here for illustration
qqnorm(model_lme2, ~ranef(., level=1))
qqnorm(model_lme2$residuals)
plot(model_lme2)

## -------------------------------------------------------------------------------------------------
# Final models

model_freq1 <- lme(fo_cog ~ overall_CT + AgeAtDiagnosis + total_anstat_cases + TimeSinceDiag, random = ~ 1 |PMC_ResearchID, 
                   method = "ML", data = fo_cog_data, na.action = na.exclude,
                   correlation = corCAR1(form = ~ 1|PMC_ResearchID),
                   weights = varExp(form = ~TimeSinceDiag)
)

summary(model_freq1)

model_duration1 <- lme(fo_cog ~ overall_CT + AgeAtDiagnosis + cumul_duration + TimeSinceDiag, random = ~ 1 |PMC_ResearchID, 
                       method = "REML", data = fo_cog_data, na.action = na.exclude,
                       correlation = corCAR1(form = ~ 1|PMC_ResearchID),
                       weights = varExp(form = ~TimeSinceDiag)
                       
)

summary(model_duration1)

model_freq2 <- lme(po_tot ~ overall_CT*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + total_anstat_cases*TimeSinceDiag, random = ~ 1 |PMC_ResearchID,
                   method = "REML", data = po_tot_data, na.action = na.exclude,
                   correlation = corCAR1(form = ~ 1|PMC_ResearchID),
                   weights = varExp(form = ~TimeSinceDiag)
                   
)

summary(model_freq2)

model_duration2 <- lme(po_tot ~ overall_CT*TimeSinceDiag + AgeAtDiagnosis*TimeSinceDiag + cumul_duration*TimeSinceDiag, random = ~ 1 |PMC_ResearchID,
                       method = "REML", data = po_tot_data, na.action = na.exclude,
                       correlation = corCAR1(form = ~ 1|PMC_ResearchID),
                       weights = varExp(form = ~TimeSinceDiag)
                       
)

summary(model_duration2)

## -------------------------------------------------------------------------------------------------
# Create table
stargazer(model_freq1, model_duration1, model_freq2, model_duration2, 
          type = "latex", 
          out = "all_27-09.htm", 
          ci = FALSE, 
          dep.var.caption = "",
          digits = 2,
          notes.align = "r",
          keep.stat = "sigma2",
          dep.var.labels=c("Cognitive Fatigue", "Overall Health-Related QOL"),
          covariate.labels = c("Constant", "CT", "Age at Diagnosis", "Cumulative Frequency", "Cumulative Duration", 
                               "Time since Diagnosis", "Cumul.Freq x Time Since Diagnosis", 
                               "Cumul.Duration x Time Since Diagnosis", "Age at Diagnosis x Time Since Diagnosis"),
          intercept.top = TRUE, intercept.bottom = FALSE,
          float.env = "sidewaystable"
)

## -------------------------------------------------------------------------------------------------
# CI's based on parametric bootstrap, use parallel processing 

# Starting a cluster with 2 cores
no_cores <- 2
cl <- makeCluster(no_cores)
registerDoParallel(cores = no_cores)

# Run 500 bootstrap iterations on each core
freq1_boot <- foreach(
  B = rep(500, 2),
  .combine = combine_lmeresamp,
  .packages = c("lmeresampler", "nlme")
) %dopar% {
  bootstrap(model_freq1, type = "parametric", B = B)
}

freq1_conf <- confint(freq1_boot, type = "basic")
freq1_conf[3:4] <- round(freq1_conf[3:4], 2)
view(freq1_conf)

dur1_boot <- foreach(
  B = rep(500, 2),
  .combine = combine_lmeresamp,
  .packages = c("lmeresampler", "nlme")
) %dopar% {
  bootstrap(model_duration1, type = "parametric", B = B)
}

dur1_conf <- confint(dur1_boot, type = "basic")
dur1_conf[3:4] <- round(dur1_conf[3:4], 2)
view(dur1_conf)

##

freq2_boot <- foreach(
  B = rep(500, 2),
  .combine = combine_lmeresamp,
  .packages = c("lmeresampler", "nlme")
) %dopar% {
  bootstrap(model_freq2, type = "parametric", B = B)
}

freq2_conf <- confint(freq2_boot, type = "basic")
freq2_conf[3:4] <- round(freq2_conf[3:4], 2)
view(freq2_conf)

dur2_boot <- foreach(
  B = rep(500, 2),
  .combine = combine_lmeresamp,
  .packages = c("lmeresampler", "nlme")
) %dopar% {
  bootstrap(model_duration2, type = "parametric", B = B)
}

dur2_conf <- confint(dur2_boot, type = "basic")
dur2_conf[3:4] <- round(dur2_conf[3:4], 2)
view(dur2_conf)

stopCluster(cl)

## -------------------------------------------------------------------------------------------------
# Plot effects

freq1_df_cog <- ggpredict(model_freq1, terms = c("total_anstat_cases","TimeSinceDiag"), type = "fixed")
freq2_df_cog <- ggpredict(model_freq1, terms = c("AgeAtDiagnosis", "TimeSinceDiag"), type = "fixed")
freq3_df_cog <- ggpredict(model_freq1, terms = c("TimeSinceDiag", "overall_CT"), type = "fixed")

dur1_df_cog <- ggpredict(model_duration1, terms = c("cumul_duration", "TimeSinceDiag"), type = "fixed")
dur2_df_cog <- ggpredict(model_duration1, terms = c("AgeAtDiagnosis", "TimeSinceDiag"), type = "fixed")
dur3_df_cog <- ggpredict(model_duration1, terms = c("TimeSinceDiag", "overall_CT"), type = "fixed")

##

freq1_df <- ggpredict(model_freq2, terms = c("total_anstat_cases", "TimeSinceDiag"), type = "fixed")
freq2_df <- ggpredict(model_freq2, terms = c("AgeAtDiagnosis", "TimeSinceDiag"), type = "fixed")
freq3_df <- ggpredict(model_freq2, terms = c("TimeSinceDiag", "overall_CT"), type = "fixed")

dur1_df <- ggpredict(model_duration2, terms = c("cumul_duration", "TimeSinceDiag"), type = "fixed")
dur2_df <- ggpredict(model_duration2, terms = c("AgeAtDiagnosis", "TimeSinceDiag"), type = "fixed")
dur3_df <- ggpredict(model_duration2, terms = c("TimeSinceDiag", "overall_CT"), type = "fixed")

int_time <- function(data, font, pred, outcome, title){
  
  if (pred %in% c("Age at diagnosis", "Cumulative duration", "Cumulative frequency")){
    numeric_group <-  as.numeric(as.character(data$group))
  }
  
  plot_int_time <- ggplot() +
    {if (pred == "CT") {geom_line(data, mapping = aes(x, predicted, col = group), size  = 0.8)
    }else if (pred %in% c("Age at diagnosis", "Cumulative duration", "Cumulative frequency")){
      geom_line(data, mapping = aes(x, predicted, col = numeric_group, 
                                    group = numeric_group), size = 0.8)
    }} +
    {if (pred == "CT"){
      labs(x = "Time since diagnosis", y = outcome, title = title) 
    }else {
      labs(x = pred, y = outcome, title = title) 
    }}+
    theme_minimal() +
    {if (pred == "CT"){
      scale_color_manual(name = str_wrap(pred, 10), labels = c("CT", "No CT"), values = c("#536CB1", "#F78C00", "#6DBF1B"))
    } else if (pred %in% c("Age at diagnosis", "Cumulative duration", "Cumulative frequency")) {
      scale_color_gradient(name = str_wrap("Time since diagnosis", 10), low = "#536CB1", high = "#F78C00")
    }
    } +
    theme(text=element_text(family=font, face="plain", size=24),
          plot.title = element_text(hjust = .5, face = "bold", size = 28, margin = unit(c(0, 0, 10, 0), "mm")),
          plot.subtitle = element_text(hjust = .5, face = "italic", size = 24),
          legend.position = "right",
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
          axis.text.x = element_text(angle = 0),
          axis.title.x = element_text(angle = 0, vjust = -3, margin = unit(c(t = 10, r = 0, l = 0, b = 0), "mm")),
          axis.title.y = element_text(angle = 90, vjust = 3, margin = unit(c(t = 0, r = 10, l = 0, b = 0), "mm")))
  
  return(plot_int_time)
}

## -------------------------------------------------------------------------------------------------







