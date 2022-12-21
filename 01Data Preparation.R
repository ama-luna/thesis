### Data Preparation ###

## -------------------------------------------------------------------------------------------------
# Packages
library(foreign)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(readxl)
library(zoo)
library(ggplot2)
library(showtext)
library(reshape2)
library(ggthemes)
library(psych)

## -------------------------------------------------------------------------------------------------

# Import data
klik_data <- read.spss("klik_July.sav", to.data.frame=T)
klik_data2 <- read.spss("new_klik_July.sav", to.data.frame=T)
klik_data2 <- klik_data2 %>% select(PMC_ResearchID, DiagnosisDateFinal, DiagnosisMainFinal, DiagnosisSubFinal)

klik_data <- merge(klik_data, klik_data2, all.x = TRUE)

klik_data <- klik_data %>% select(-contains(c("pat", "PAT", "RISICO", "risico", "PROBLEMEN", "subscore", "ouders")))

################### Pre-processing the data #################

# Convert dates since the format was not preserved
klik_data$StartDatumEersteDBC <- as.Date(klik_data$StartDatumEersteDBC/86400, origin="1582-10-14")
klik_data$Entry_Date <- as.Date(klik_data$Entry_Date/86400, origin="1582-10-14")
klik_data$EersteDiagnoseDatum <- as.Date(klik_data$EersteDiagnoseDatum/86400, origin="1582-10-14")
klik_data$DiagnosisDateFinal <- as.numeric(as.character(klik_data$DiagnosisDateFinal))
klik_data$DiagnosisDateFinal <- as.Date(klik_data$DiagnosisDateFinal/86400, origin="1582-10-14")

# Remove extra spaces on strings
klik_data$Iccc_Main <- str_squish(klik_data$Iccc_Main)
klik_data$DBCDiagnoseCodeOmschrijving <- str_squish(klik_data$DBCDiagnoseCodeOmschrijving)
klik_data$AgeGroup <- str_squish(klik_data$AgeGroup)

klik_data <- klik_data %>% filter(AgeGroup == "<8 jaar")

## -------------------------------------------------------------------------------------------------

# Check number of cases to be included, based on date
ids <- unique(klik_data$PMC_ResearchID)

klik_filtered <- klik_data %>% filter(DiagnosisDateFinal >= "2019-01-01")
length(unique(klik_filtered$PMC_ResearchID))

# Remove rows with only NA
klik_filtered <- klik_filtered[rowSums(is.na(klik_filtered)) != ncol(klik_filtered), ]

## -------------------------------------------------------------------------------------------------

#### Anstat merge ####
anstat_data <- read.spss("anstat_July.sav", to.data.frame=T)

# Remove recovery rooms
anstat_data <- anstat_data %>%
  filter(!str_detect(kamer, "rec"))

anstat_data$begintijd <- as.POSIXlt(anstat_data$begintijd)
anstat_data$eindtijd <- as.POSIXlt(anstat_data$eindtijd)

library("timetk")

dates <- anstat_data %>%
  group_by(PMC_ResearchID) %>%
  summarise_by_time(
    .date_var = begintijd,
    .by       = "day", 
    value  = n()
  )

# Count rows per participant, since each row indicates a different session
anstat_cases <- dates %>% group_by(PMC_ResearchID) %>% summarise(total_anstat_cases = n())

# Merge datasets
klik_filtered <- merge(klik_filtered, anstat_cases, by = "PMC_ResearchID", all.x = TRUE)
anstat_new <- merge(anstat_data, anstat_cases, by = "PMC_ResearchID", all.x = TRUE)

# Count number of NAs
anstat_cases_count <- klik_filtered %>% group_by(PMC_ResearchID) %>% summarise(cases = sum(is.na(total_anstat_cases)))
sum(anstat_cases_count$cases > 0)

# Duration calculation
# Calculate difference in time between begintijd and eindtijd

data <- anstat_data
data$duration <- NA
  for (i in unique(data$PMC_ResearchID)){
    for (j in 1:nrow(data[data$PMC_ResearchID == i,])){
      data[data$PMC_ResearchID == i,]$duration[j] <- abs(as.numeric(difftime(data[data$PMC_ResearchID == i,]$eindtijd[j], data[data$PMC_ResearchID == i,]$begintijd[j], units = "hours")))
    }
  }

# Cumulative duration
cumulative_duration <- data %>% group_by(PMC_ResearchID) %>% summarise(cumul_duration = sum(duration))

klik_filtered <- merge(klik_filtered, cumulative_duration, by = "PMC_ResearchID", all.x = TRUE)

# Average duration
library(lubridate)

new_data_anstat <- merge(data, anstat_cases, by = "PMC_ResearchID", all.x = TRUE)

new_new_anstat_data <- new_data_anstat
new_new_new_anstat_data <- new_new_anstat_data %>% group_by(PMC_ResearchID) %>% summarise(avg_duration = mean(duration))

avg_duration <- new_new_new_anstat_data %>% select(PMC_ResearchID, avg_duration)

klik_filtered <- merge(klik_filtered, avg_duration, by = "PMC_ResearchID", all.x = TRUE)

# new_data_anstat <- merge(data, anstat_cases, by = "PMC_ResearchID", all.x = TRUE)
new_data_anstat <- merge(new_data_anstat, avg_duration, by = "PMC_ResearchID", all.x = TRUE)
new_data_anstat <- merge(new_data_anstat, cumulative_duration, by = "PMC_ResearchID", all.x = TRUE)

## -------------------------------------------------------------------------------------------------

##### Fo_cog ######

# Count number of cases of fo_cog
fo_cog_na <- klik_filtered %>% group_by(PMC_ResearchID) %>% summarize(fo_cog = sum(!is.na(fo_cog)))

# Number of cases in total
length(fo_cog_na$fo_cog)
sum(fo_cog_na$fo_cog >= 2)

# Only select those participants with at least two occurrences of the questionnaire data
ids_fo_cog <- unique(klik_filtered$PMC_ResearchID)[fo_cog_na$fo_cog >= 2]
klik_data_fo_cog <- klik_filtered[klik_filtered$PMC_ResearchID %in% ids_fo_cog,]

length(unique(klik_data_fo_cog$PMC_ResearchID))
anstat_cases_count <- klik_data_fo_cog %>% group_by(PMC_ResearchID) %>% summarise(cases = sum(is.na(total_anstat_cases)))
sum(anstat_cases_count$cases > 0)


###### Po_tot ########

# Count number of missings in po_tot per participant
po_tot_na <- klik_filtered %>% group_by(PMC_ResearchID) %>% summarize(po_tot = sum(!is.na(po_tot)))

# Number of cases in total
sum(po_tot_na$po_tot)
sum(po_tot_na$po_tot >= 2)

# Only select those participants with at least two occurrences of the questionnaire data
ids_po_tot <- unique(klik_filtered$PMC_ResearchID)[po_tot_na$po_tot >= 2]
klik_data_po_tot <- klik_filtered[klik_filtered$PMC_ResearchID %in% ids_po_tot,]

length(unique(klik_data_po_tot$PMC_ResearchID))
anstat_cases_count <- klik_data_po_tot %>% group_by(PMC_ResearchID) %>% summarise(cases = sum(is.na(total_anstat_cases)))
sum(anstat_cases_count$cases > 0)

########

fo_cog_full <- klik_data_fo_cog
po_tot_full <- klik_data_po_tot

## -------------------------------------------------------------------------------------------------

# Check number of rows with only NA values
fo_cog_unique <- fo_cog_full %>% distinct(PMC_ResearchID, .keep_all = TRUE)

# Count 
sum(rowSums(is.na(fo_cog_unique)) == ncol(fo_cog_unique))

# Remove rows with only NA values
fo_cog_full_filtered <- fo_cog_full[rowSums(is.na(fo_cog_full)) != ncol(fo_cog_full), ]

######

po_tot_unique <- po_tot_full %>% distinct(PMC_ResearchID, .keep_all = TRUE)

sum(rowSums(is.na(po_tot_unique)) == ncol(po_tot_unique))

po_tot_full_filtered <- po_tot_full[rowSums(is.na(po_tot_full)) != ncol(po_tot_full), ]

## -------------------------------------------------------------------------------------------------

fo_cog_full <- fo_cog_full_filtered
po_tot_full <- po_tot_full_filtered

# Prepare dates 
fo_cog_full$Geboortemaandjaar <- as.yearmon(fo_cog_full$Geboortemaandjaar, "%Y-%m")
fo_cog_full$Geboortemaandjaar <- as.Date(fo_cog_full$Geboortemaandjaar)

po_tot_full$Geboortemaandjaar <- as.yearmon(po_tot_full$Geboortemaandjaar, "%Y-%m")
po_tot_full$Geboortemaandjaar <- as.Date(po_tot_full$Geboortemaandjaar)

# Patients with wrongly specified diagnosis dates
fo_cog_full[fo_cog_full$PMC_ResearchID == "PMCID932AAL",]$DiagnosisDateFinal <- rep("2019-09-01", length(fo_cog_full[fo_cog_full$PMC_ResearchID == "PMCID932AAL",]$DiagnosisDateFinal))
po_tot_full[po_tot_full$PMC_ResearchID == "PMCID932AAL",]$DiagnosisDateFinal <- rep("2019-09-01", length(po_tot_full[po_tot_full$PMC_ResearchID == "PMCID932AAL",]$DiagnosisDateFinal))

fo_cog_full[fo_cog_full$PMC_ResearchID == "PMCID290AAA",]$DiagnosisDateFinal <- rep("2016-11-01", length(fo_cog_full[fo_cog_full$PMC_ResearchID == "PMCID290AAA",]$DiagnosisDateFinal))
po_tot_full[po_tot_full$PMC_ResearchID == "PMCID290AAA",]$DiagnosisDateFinal <- rep("2016-11-01", length(po_tot_full[po_tot_full$PMC_ResearchID == "PMCID290AAA",]$DiagnosisDateFinal))

fo_cog_full$DiagnosisDateFinal <- as.yearmon(fo_cog_full$DiagnosisDateFinal, "%Y-%m")
fo_cog_full$DiagnosisDateFinal <- as.Date(fo_cog_full$DiagnosisDateFinal)

po_tot_full$DiagnosisDateFinal <- as.yearmon(po_tot_full$DiagnosisDateFinal, "%Y-%m")
po_tot_full$DiagnosisDateFinal <- as.Date(po_tot_full$DiagnosisDateFinal)

fo_cog_full$Entry_Date <- as.yearmon(fo_cog_full$Entry_Date, "%Y-%m")
fo_cog_full$Entry_Date <- as.Date(fo_cog_full$Entry_Date)

po_tot_full$Entry_Date <- as.yearmon(po_tot_full$Entry_Date, "%Y-%m")
po_tot_full$Entry_Date <- as.Date(po_tot_full$Entry_Date)

# Function to calculate difference in time in months, which is not possible with difftime, source: https://stackoverflow.com/questions/25495833/extending-difftime-to-include-months-and-years
timeDiff <- function(dob, dod) 
{
  s <- sapply(list("months"), function(x) {
    (length(seq(dob, dod, x)) - 1)
  })
  setNames(s, c("months"))
}

# Calculate age at diagnosis
add_ageatdiag <- function(data){
  
  AgeAtDiagnosis <- list()
  for (i in unique(data$PMC_ResearchID)){
    AgeAtDiagnosis[[i]] <- age_calc(data[data$PMC_ResearchID == i,]$Geboortemaandjaar[2], data[data$PMC_ResearchID == i,]$DiagnosisDateFinal[2])
  }
  
  library(tibble)
  AgeAtDiagnosis <- as.data.frame(as.matrix(AgeAtDiagnosis))
  AgeAtDiagnosis <- rownames_to_column(AgeAtDiagnosis)
  
  AgeAtDiagnosis$V1 <- as.numeric(AgeAtDiagnosis$V1 )
  colnames(AgeAtDiagnosis) <- c("PMC_ResearchID", "AgeAtDiagnosis")
  AgeAtDiagnosis$AgeAtDiagnosis <- AgeAtDiagnosis$AgeAtDiagnosis/12 
  
  new_data <- merge(data, AgeAtDiagnosis, by = "PMC_ResearchID", all.x = TRUE)
  
  return(new_data)
  
}

library(eeptools)

# Calculate time since diagnosis
add_timesincediag <- function(data){
  
  data$TimeSinceDiag <- NA
  for (i in unique(data$PMC_ResearchID)){
    for (j in 1:nrow(data[data$PMC_ResearchID == i,])){
      
      if (data[data$PMC_ResearchID == i,]$DiagnosisDateFinal[j] < data[data$PMC_ResearchID == i,]$Entry_Date[j]){
        
        data[data$PMC_ResearchID == i,]$TimeSinceDiag[j] <- round(age_calc(data[data$PMC_ResearchID == i,]$DiagnosisDateFinal[j], data[data$PMC_ResearchID == i,]$Entry_Date[j]), 1)
        
      } else if (data[data$PMC_ResearchID == i,]$DiagnosisDateFinal[j] > data[data$PMC_ResearchID == i,]$Entry_Date[j]){
        
        data[data$PMC_ResearchID == i,]$TimeSinceDiag[j] <- round(age_calc(data[data$PMC_ResearchID == i,]$Entry_Date[j], data[data$PMC_ResearchID == i,]$DiagnosisDateFinal[j]), 1)
        
      }
      
      else if (data[data$PMC_ResearchID == i,]$DiagnosisDateFinal[j] == data[data$PMC_ResearchID == i,]$Entry_Date[j]){
        
        data[data$PMC_ResearchID == i,]$TimeSinceDiag[j] <- 0
        
      }
    }
  }
  
  return(data)
}

new_fo_cog <- add_ageatdiag(fo_cog_full)
new_po_tot <- add_ageatdiag(po_tot_full)

new_fo_cog <- add_timesincediag(new_fo_cog)
new_po_tot <- add_timesincediag(new_po_tot)

# Fill empty values that are not actually missing 

fo_cog_data <- new_fo_cog
po_tot_data <- new_po_tot

fo_cog_data <- fo_cog_data %>% select(PMC_ResearchID, Geslacht, AgeAtDiagnosis, Entry_Date, Geboortemaandjaar, CT, RT, SCT, fo_cog, po_tot, DiagnosisDateFinal, "DiagnosisDateFinal", "DiagnosisMainFinal", "DiagnosisSubFinal", 
                                      "total_anstat_cases", "cumul_duration", "avg_duration", "AgeAtDiagnosis", 
                                      "TimeSinceDiag", "overall_CT")

po_tot_data <- po_tot_data %>% select(PMC_ResearchID, Geslacht, AgeAtDiagnosis, Entry_Date, Geboortemaandjaar, CT, RT, SCT, fo_cog, po_tot, DiagnosisDateFinal, "DiagnosisDateFinal", "DiagnosisMainFinal", "DiagnosisSubFinal", 
                                      "total_anstat_cases", "cumul_duration", "avg_duration", "AgeAtDiagnosis", 
                                      "TimeSinceDiag", "overall_CT")


fo_cog_data <- fo_cog_data %>% distinct(.keep_all = T)
po_tot_data <- po_tot_data %>% distinct(.keep_all = T)

# Check double rows (double entries)
fo_cog_data %>%
  group_by(Entry_Date, PMC_ResearchID) %>%
  select(PMC_ResearchID, fo_cog, po_tot, Entry_Date,DiagnosisDateFinal,  TimeSinceDiag) %>%
  filter(n() > 1) %>%
  ungroup

View(po_tot_data %>%
       group_by(Entry_Date, PMC_ResearchID) %>%
       select(PMC_ResearchID, fo_cog, po_tot, Entry_Date,DiagnosisDateFinal,  TimeSinceDiag) %>%
       filter(n() > 1) %>%
       ungroup
)

# Fill missings, only keep one row
my_fun <- function(x) x[!is.na(x)]

f <- function(x) {
  x <- na.omit(x)
  if (length(x) > 0) last(x) else NA
}

fo_cog_data <- fo_cog_data %>%  group_by_(.dots = as.list(names(.)[!(names(.) %in% c("fo_cog", "po_tot"))])) %>% summarise_all(funs(f))
po_tot_data <- po_tot_data %>%  group_by_(.dots = as.list(names(.)[!(names(.) %in% c("fo_cog", "po_tot"))])) %>% summarise_all(funs(f))

fo_cog_data %>%
  group_by(Entry_Date, PMC_ResearchID) %>%
  select(PMC_ResearchID, fo_cog, po_tot, Entry_Date,DiagnosisDateFinal,  TimeSinceDiag) %>%
  filter(n() > 1) %>%
  ungroup

new_fo_cog <- fo_cog_data
new_po_tot <- po_tot_data

## -------------------------------------------------------------------------------------------------

### Until here, filtering on date of diagnosis, age at diagnosis was calculated, time variables were calculated

# Checking diagnoses
new_po_tot %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
  group_by(DiagnosisSubFinal) %>% summarise(Frequency = n())%>%
  mutate(Percentage = Frequency/sum(Frequency)*100) %>% arrange(-Frequency) 

new_fo_cog %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
  group_by(DiagnosisSubFinal) %>% summarise(Frequency = n())%>%
  mutate(Percentage = Frequency/sum(Frequency)*100) %>% arrange(-Frequency) 

# Finding IDs in both datasets
'%!in%' <- function(x,y)!('%in%'(x,y))
unique(po_tot_full$PMC_ResearchID)[which(unique(fo_cog_full$PMC_ResearchID) %!in% unique(po_tot_full$PMC_ResearchID))]

## -------------------------------------------------------------------------------------------------

# Descriptives

new_fo_cog <- new_fo_cog %>% arrange(PMC_ResearchID, Entry_Date)
new_po_tot <- new_po_tot %>% arrange(PMC_ResearchID, Entry_Date)

descriptives <- function(data){
  
  # Gender descriptives
  geslacht_descr <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    group_by(Geslacht) %>% summarise(Frequency = n())%>%
    mutate(Percentage = Frequency/sum(Frequency)*100)
  
  p1 <- geslacht_descr %>% 
    ggplot(aes(x=Geslacht, y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + labs(title = "Gender")
  
  # Age descriptives
  age_descr <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    summarise(Mean = mean(AgeAtDiagnosis, na.rm=T),
              SD = sd(AgeAtDiagnosis, na.rm=T))
  
  age_descr2 <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    group_by(AgeAtDiagnosis) %>% summarise(Frequency = n())%>%
    mutate(Percentage = Frequency/sum(Frequency)*100)
  
  p1.2 <- age_descr2 %>% 
    ggplot(aes(x=AgeAtDiagnosis, y= Frequency)) +
    geom_bar(stat="identity", binwidth = .3, fill="#163973") + theme_minimal() + labs(title = "Age at Diagnosis")
  
  # Type of diagnose descriptives
  diagnose_descr <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    group_by(DBCDiagnoseCodeOmschrijving) %>% summarise(Frequency = n())%>%
    mutate(Percentage = Frequency/sum(Frequency)*100) %>% arrange(-Frequency)
  
  p2 <- diagnose_descr %>% 
    ggplot(aes(x=reorder(DBCDiagnoseCodeOmschrijving, -Frequency), y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + 
    labs(title = "Diagnosis", x ="Diagnosis") + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
  
  # CT descriptives
  ct_descr <- data %>% select(PMC_ResearchID, CT) %>%
    arrange(PMC_ResearchID, CT) %>% 
    group_by(PMC_ResearchID) %>% 
    summarise_all(last) %>%ungroup() %>% group_by(CT) %>% 
    summarise(Frequency = n())
  
  p3 <- ct_descr %>% 
    ggplot(aes(x=CT, y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + labs(title = str_wrap("CT exposure before questionnaire", 30))
  
  # SCT descriptives
  sct_descr <- data %>% select(PMC_ResearchID, SCT) %>%
    arrange(PMC_ResearchID, SCT) %>% 
    group_by(PMC_ResearchID) %>% 
    summarise_all(last) %>%ungroup() %>% group_by(SCT) %>% 
    summarise(Frequency = n())
  
  p4 <- sct_descr %>% 
    ggplot(aes(x=SCT, y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + labs(title = str_wrap("Stem cell treatment exposure before questionnaire", 30))
  
  # RT descriptives
  rt_descr <- data %>% select(PMC_ResearchID, RT) %>%
    arrange(PMC_ResearchID, RT) %>% 
    group_by(PMC_ResearchID) %>% 
    summarise_all(last) %>%ungroup() %>% group_by(RT) %>% 
    summarise(Frequency = n())
  
  p5 <- rt_descr %>% 
    ggplot(aes(x=RT, y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + labs(title = str_wrap("Radiation therapy exposure before questionnaire", 30))
  
  # Type of tumor descriptives
  tumor_descr <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    group_by(Iccc_Main) %>% summarise(Frequency = n())%>%
    mutate(Percentage = Frequency/sum(Frequency)*100) %>% arrange(-Frequency)
  
  diag_descr <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    group_by(DiagnosisMainFinal, DiagnosisSubFinal) %>% summarise(Frequency = n())%>%
    mutate(Percentage = Frequency/sum(Frequency)*100) %>% arrange(-Frequency) 
  
  p6 <- tumor_descr %>% 
    ggplot(aes(x=reorder(Iccc_Main, -Frequency), y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + 
    labs(title = "Tumor Type", x ="Tumor type") + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
  
  # Anstat
  total_anstat_dat <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE)
  total_anstat <- sum(is.na(total_anstat_dat$total_anstat_cases))
  anstat_atleast2 <- sum(total_anstat_dat$total_anstat_cases[!is.na(total_anstat_dat$total_anstat_cases)] >= 2)
  
  anstat_descr <- data %>% distinct(PMC_ResearchID, .keep_all = TRUE) %>% 
    group_by(total_anstat_cases) %>% summarise(Frequency = n())%>%
    mutate(Percentage = Frequency/sum(Frequency)*100) %>% arrange(-Frequency)
  
  p7 <- anstat_descr %>% 
    ggplot(aes(x=reorder(total_anstat_cases, -Frequency), y= Frequency)) +
    geom_col(fill="#163973") + theme_minimal() + 
    labs(title = "Number of AnStat cases", x ="Number of cases") + 
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 30))
  
  all_plots <- ggarrange(p1, p1.2, p3, p4, p5, p7)
  type_plots <- ggarrange(p2, p6)
  
  return(list(geslacht_descr, age_descr, diag_descr, ct_descr, sct_descr, rt_descr, all_plots, type_plots, total_anstat))
  
}

fo_cog_final <- descriptives(new_fo_cog)
po_tot_final <- descriptives(new_po_tot)

## -------------------------------------------------------------------------------------------------

new_fo_cog %>%  arrange(PMC_ResearchID, AgeAtQuestionnaire) %>% group_by(PMC_ResearchID) %>% 
  summarise_all(last) %>%
  summarise(mean_age = mean(AgeAtQuestionnaire)) %>% 
  summarise(mean_age_overall = mean(mean_age)) 

new_po_tot %>% arrange(PMC_ResearchID, AgeAtQuestionnaire) %>% group_by(PMC_ResearchID) %>% 
  summarise_all(last) %>%
  summarise(mean_age = mean(AgeAtQuestionnaire)) %>% 
  summarise(mean_age_overall = mean(mean_age))

## -------------------------------------------------------------------------------------------------

# Correlation between duration vs. frequency
cor(new_fo_cog$cumul_duration, new_fo_cog$total_anstat_cases) 
cor(new_po_tot$cumul_duration, new_po_tot$total_anstat_cases) 

plot_data <- new_po_tot %>%
  arrange(PMC_ResearchID) %>% 
  group_by(PMC_ResearchID) %>% 
  summarise_all(last)

ggplot(data = plot_data, aes(x=avg_duration, y=po_tot)) + geom_point()

## -------------------------------------------------------------------------------------------------

# Descriptives 

# Avg duration fo_cog
new_fo_cog %>% group_by(PMC_ResearchID) %>%
  summarise(mean_val = mean(avg_duration)) %>% 
  summarise(mean_val_overall = mean(mean_val),
            sd = sd(mean_val))

# Cumul duration fo_cog
new_fo_cog %>% group_by(PMC_ResearchID) %>%
  summarise(mean_val = mean(cumul_duration)) %>% 
  summarise(mean_val_overall = mean(mean_val),
            sd = sd(mean_val))

# Frequency fo_cog
new_fo_cog %>% group_by(PMC_ResearchID) %>%
  summarise(mean_val = mean(total_anstat_cases)) %>% 
  summarise(mean_val_overall = mean(mean_val),
            sd = sd(mean_val))

###

# Avg duration po_tot
new_po_tot %>% group_by(PMC_ResearchID) %>%
  summarise(mean_val = mean(avg_duration)) %>% 
  summarise(mean_val_overall = mean(mean_val),
            sd = sd(mean_val))

# Cumul duration po_tot
new_po_tot %>% group_by(PMC_ResearchID) %>%
  summarise(mean_val = mean(cumul_duration)) %>% 
  summarise(mean_val_overall = mean(mean_val),
            sd = sd(mean_val))

# Frequency po_tot
new_po_tot %>% group_by(PMC_ResearchID) %>%
  summarise(mean_val = mean(total_anstat_cases)) %>% 
  summarise(mean_val_overall = mean(mean_val),
            sd = sd(mean_val))

## -------------------------------------------------------------------------------------------------

# Checking distributions, outliers of duration, frequency etc.

# Extract relevant PMCIDs
selected_data <- new_data_anstat[new_data_anstat$PMC_ResearchID %in% new_fo_cog$PMC_ResearchID | new_data_anstat$PMC_ResearchID %in% new_po_tot$PMC_ResearchID,]

# Create distinct data for histograms
distinct_anstat_data <- new_data_anstat %>% distinct(PMC_ResearchID, .keep_all = TRUE)
hist(distinct_anstat_data$total_anstat_cases)
hist(distinct_anstat_data$cumul_duration)
hist(distinct_anstat_data$avg_duration)
hist(new_fo_cog$AgeAtDiagnosis)
hist(new_po_tot$AgeAtDiagnosis)

hist(new_fo_cog$fo_cog)
hist(new_po_tot$po_tot)

test <- new_fo_cog %>% group_by(PMC_ResearchID) %>% 
  summarise_all(first)

hist(test$fo_cog)

# Summary of variables
summary(new_data_anstat$cumul_duration)
summary(new_data_anstat$avg_duration)
summary(new_data_anstat$total_anstat_cases)
summary(new_data_anstat$duration)

# Arrange data on max. and min. to check whether outliers are actually weird values
View(selected_data %>% arrange(-total_anstat_cases))
View(selected_data %>% arrange(-duration))
View(selected_data %>% arrange(-cumul_duration))
View(selected_data %>% arrange(-avg_duration))

## -------------------------------------------------------------------------------------------------

# Checking CRT

table(new_fo_cog$DiagnosisMainFinal, new_fo_cog$RT)

new_fo_cog %>% filter(DiagnosisMainFinal == "CNS tumor") %>% select(PMC_ResearchID, RT) %>%
  arrange(PMC_ResearchID, RT) %>% 
  group_by(PMC_ResearchID) %>% 
  summarise_all(last) %>%ungroup() %>% group_by(RT) %>% 
  summarise(Frequency = n())

new_po_tot %>% filter(DiagnosisMainFinal == "CNS tumor") %>% select(PMC_ResearchID, RT) %>%
  arrange(PMC_ResearchID, RT) %>% 
  group_by(PMC_ResearchID) %>% 
  summarise_all(last) %>%ungroup() %>% group_by(RT) %>% 
  summarise(Frequency = n())

## -------------------------------------------------------------------------------------------------

# Time points visualization

input_data <- fo_cog_data %>% ungroup()

plotting_time_points <- function(input_data, outcome, font){
  
  pmc_data <- input_data %>% select(PMC_ResearchID, Entry_Date, DiagnosisMainFinal)
  
  # Cast the data to have individuals as columns and dates as row
  cast_data <- dcast(pmc_data, Entry_Date~PMC_ResearchID,value.var="PMC_ResearchID")
  rownames(cast_data) <- cast_data$Entry_Date
  cast_data <- cast_data[,-1]
  cast_data[is.na(cast_data)] <- 0
  cast_data[cast_data != 0] <- 1
  
  obs_data <- apply(cast_data,2,function(x){
    seen_dates <- which(x == 1 | x == 2)
    x[seen_dates] <- "seen"
    x[setdiff(seq(seen_dates[1],tail(seen_dates,1)),seen_dates)] <- "not_seen"
    x[setdiff(seq(seen_dates[1],tail(seen_dates,1)),seen_dates)] <- NA
    other <- which(x == 0)
    x[other] <- NA
    x
  })
  
  # Melt the data and add a first_seen column that has the year of first sighting
  data<-as.data.frame(melt(obs_data))
  head(data)
  colnames(data) <- c("obs_date","ind_id","obs")
  
  data<-data[!is.na(data$obs),]
  data$obs_date<-as.Date(data$obs_date)
  for (i in unique(data$ind_id)){
    data[data$ind_id==i,4] <- format(data[data$ind_id==i,1][1],"%Y")
  }
  colnames(data) <- c("obs_date","PMC_ResearchID","obs","first_seen")
  
  diag_info <- pmc_data %>% select(PMC_ResearchID, DiagnosisMainFinal)
  
  data <- merge(data, diag_info, by = "PMC_ResearchID")
  data <- data[with(data, order(obs_date, first_seen)),]
  data$PMC_ResearchID <- factor(data$PMC_ResearchID, levels = unique(data$PMC_ResearchID))
  
  if (outcome == "Cognitive Fatigue"){
    plot_time_points <- ggplot(data,aes(obs_date, 
                                        reorder(PMC_ResearchID, first_seen),colour=first_seen)) +
      facet_wrap(~DiagnosisMainFinal) +
      geom_point(aes(shape=obs),size=2.5) +
      scale_x_date() +
      geom_line(size = 1, alpha = 0.3) +
      theme_few() +
      guides(shape="none") +
      labs(y = "ID", x = "Observation date", fill = "Year of Start", 
           title = paste0("Overview of Time Points: ", outcome), 
           subtitle = "Time points of entries for each participant by main diagnosis group") +
      scale_fill_manual(name= "Year of First Entry", values = alpha(c("#536CB1","#F78C00", "#6DBF1B", "#B266FF"), 0.7)) +
      scale_colour_manual(name= "Year of First Entry", values = alpha(c("#536CB1","#F78C00", "#6DBF1B", "#B266FF"), 1)) +
      theme(text=element_text(family=font, face="plain", size=24),
            plot.title = element_text(hjust = .5, face = "bold", size = 28),
            plot.subtitle = element_text(hjust = .5, face = "italic", size = 24),
            legend.position = "right",
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            legend.background=element_blank(),
            plot.margin=grid::unit(c(20,20,20,20), "mm"),
            panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                              colour = "#F5F5F5"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = "white"),
            panel.border = element_rect(color = "lightgrey"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 0),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(angle = 0, vjust = -3),
            axis.title.y = element_text(angle = 90, vjust = 3))
  }
  
  else {
    plot_time_points <- ggplot(data,aes(obs_date, 
                                        reorder(PMC_ResearchID, first_seen),colour=first_seen)) +
      facet_wrap(~DiagnosisMainFinal) +
      geom_point(aes(shape=obs),size=3) +
      scale_x_date() +
      geom_line(size = 1.5, alpha = 0.3) +
      theme_few() +
      guides(shape="none") +
      labs(y = "ID", x = "Observation date", fill = "Year of Start", 
           title = paste0("Overview of Time Points: ", outcome), 
           subtitle = "Time points of entries for each participant by main diagnosis group") +
      scale_fill_manual(name= "Year of First Entry", values = alpha(c("#536CB1","#F78C00", "#6DBF1B", "#B266FF"), 0.7)) +
      scale_colour_manual(name= "Year of First Entry", values = alpha(c("#536CB1","#F78C00", "#6DBF1B", "#B266FF"), 1)) +
      theme(text=element_text(family="Palatino_plot", face="plain", size=28),
            plot.title = element_text(hjust = .5, face = "bold", size = 32),
            plot.subtitle = element_text(hjust = .5, face = "italic", size = 28),
            legend.position = "right",
            plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white",
                                            colour = "white",
                                            size = 0.5, linetype = "solid"),
            legend.background=element_blank(),
            plot.margin=grid::unit(c(20,20,20,20), "mm"),
            #t, r, b, l
            panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                              colour = "#F5F5F5"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(fill = "white"),
            panel.border = element_rect(color = "lightgrey"),
            axis.ticks.x = element_blank(),
            axis.text.x = element_text(angle = 0),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(angle = 0, vjust = -3),
            axis.title.y = element_text(angle = 90, vjust = 3))
  }
  
  return(plot_time_points)
}

## -------------------------------------------------------------------------------------------------

# Visualize correlations

plotting_correlations <- function(data, oucome){
  
  data <- data %>% select(PMC_ResearchID, AgeAtDiagnosis, TimeSinceDiag, total_anstat_cases, avg_duration, cumul_duration)
  data <- data %>% group_by(PMC_ResearchID) %>% 
    summarise_all(last) 
  
  cors <- cor(data[,-1])
  
  colnames(cors) <- c("Age at diagnosis", "Time since diagnosis", "Cumulative frequency",
                      "Average duration", "Cumulative duration")
  rownames(cors) <- c("Age at diagnosis", "Time since diagnosis", "Cumulative frequency",
                      "Average duration", "Cumulative duration")
  
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  lower_tri <- get_lower_tri(cors)
  melted_cormat <- melt(lower_tri)
  
  corr_plot <- ggplot(data = subset(melted_cormat, !is.na(value)), aes(x=Var1, y=Var2, fill=value, label=value)) + 
    geom_tile(color = "white") +
    geom_text(aes(label = format(round(value, 3), nsmall = 3)), 
              family="Palatino_plot", size = 18*0.36) +
    labs(title = paste0("Correlations of Variables: ", oucome), subtitle = "Correlations between numeric variables") +
    scale_fill_gradient2(low = "#F78C00", high = "#0242B1", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation", 
                         na.value="transparent") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = .5, face = "bold", size = 24),
          plot.subtitle = element_text(hjust = .5, face = "italic", size = 20),
          plot.margin=grid::unit(c(10,10,10,10), "mm"),
          text=element_text(family="Palatino_plot", face="plain", size=20))+
    coord_fixed()
  
  return(corr_plot)
  
  
}

