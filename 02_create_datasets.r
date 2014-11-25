

# Copyright 2014 University of Washington
# 
# This file is part of BRFSS-SAE.
# 
# BRFSS-SAE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BRFSS-SAE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with BRFSS-SAE.  If not, see <http://www.gnu.org/licenses/>.

########################################################################################################################
## Description: Prep all datasets needed to fit each model
##              (1) Full data for final model
##              (2) Sampled down data for validation models
##              (3) Imputed datasets when there is a correction being used
########################################################################################################################

library(mvtnorm)
library(boot)
library(plyr)
library(data.table) 

rm(list=ls())
set.seed(98121)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
setwd(paste(root, "Project/NHIS/sae", sep=""))

## Get settings
outcome <- commandArgs()[3]

try(detach(settings), silent=T)
get_settings <- read.csv("05_small_area_models/model_settings.csv", skip=1, stringsAsFactors=F)[,-1]
settings <- lapply(1:nrow(get_settings), function(x) strsplit(get_settings[x, outcome], split=";")[[1]])
settings <- lapply(settings, function(x) if(sum(is.na(as.numeric(x)))==0) as.numeric(x) else x)
names(settings) <- get_settings[,1]
attach(settings)
rm(get_settings, settings)
 
## Load BRFSS data; keep only years in range 
load("03_survey_data/output/brfss_microdata.rdata")  
data <- data.table(data[data$year %in% year1:year2 & !is.na(data$age) & data$age >= min(age_groups) & data$age < max(age_groups),])
if (!is.na(skip_years[1])) data <- data[!year %in% skip_years,]
  
## Format phone variable
data[ownll == 1 & owncp == 1 & year > 2010, phone := "dual_phone"]
data[ownll == 1 & owncp == 0 & year > 2010, phone := "landline_only"]
data[ownll == 0 & owncp == 1 & year > 2010, phone := "cell_only"] 
data[, phone := factor(phone, levels=c("dual_phone", "landline_only", "cell_only"))] 

## Drop if sex, age, race, marital status, or phone is missing
data <- data[!is.na(sex) & !is.na(age),] 
if (1 %in% model_race) data <- data[!is.na(race),]
if (1 %in% model_marriage) data <- data[!is.na(marital),]
if (1 %in% model_edu) data <- data[!is.na(edu),]
if (correct_cells == 1) data <- data[!is.na(phone) | year <= 2010,]

## Recode missing fips
data[fips %% 1000 == 999, fips := NA]

## Format and (if necessary) correct outcome variable
# alcohol outcomes 
if (grepl("^any_drinking", outcome)) { 
  if (correction == 0) data[, outcome.0 := data$anydrnk]
  if (correction == 1) stop("*** No correction available for drinking ***") 
} 

if (grepl("^heavy_drinking", outcome)) { 
  if (correction == 0) { 
    data[, bingedrnk := pmin(30, bingedrnk)/30]
    data[, AF := pmax(0, alcday - bingedrnk)]
    data[, AQ := avedrnk]
    data[, BF := bingedrnk]
    data[, BQ := maxdrnks]
    data[, csmp := pmax(alcday*avedrnk, AF*AQ + BF*BQ)]
    data[, outcome.0 := as.numeric(csmp > ifelse(sex == 1, 2, 1))]
  } 
  if (correction == 1) stop("*** No correction available for heavy drinking ***") 
} 

if (grepl("^binge_drinking", outcome)) { 
  if (correction == 0) data[, outcome.0 := as.numeric(bingedrnk >= 1)]
  if (correction == 1) stop("*** No correction available for binge drinking ***") 
} 

if (grepl("^excess_drinking", outcome)) { 
  if (correction == 0) {
    data[, outcome.0 := as.numeric((bingedrnk >= 1) | (heavydrnk == 1))]
    data[is.na(bingedrnk) | is.na(heavydrnk), outcome.0 := NA]
  } 
  if (correction > 0) stop("*** No correction available for excessive drinking ***") 
} 
  
# self-reported health outcomes
if (grepl("^poor_hlth", outcome)) { 
  if (correction == 0) data[, outcome.0 := as.numeric(genhlth == "Fair" | genhlth == "Poor")]
  if (correction == 1) stop("*** No correction available for self-reported health status")
} 

if (grepl("^poor_phys_days", outcome)) { 
  if (correction == 0) data[, outcome.0 := as.numeric(physhlth >= 14)]
  if (correction == 1) stop("*** No correction available for self-reported unhealthy physical days")
} 

if (grepl("^poor_ment_days", outcome)) { 
  if (correction == 0) data[, outcome.0 := as.numeric(menthlth >= 14)]
  if (correction == 1) stop("*** No correction available for self-reported unhealthy mental days")
}      

# diabetes outcomes
if (grepl("^diabetes", outcome)) {
  # make sure correction is set properly
  if (grepl("^diabetes_diagnosed", outcome) & correction == 1) stop("*** No correction available for diagnosed diabetes")
  if (grepl("^diabetes_undiagnosed", outcome) & correction == 0) stop("*** Can't run undiagnosed diabetes without a correction")
  if (grepl("^diabetes_controlled", outcome) & correction == 0) stop("*** Can't run controlled diabetes without a correction")
  if (grepl("^diabetes_uncontrolled", outcome) & correction == 0) stop("*** Can't run uncontrolled diabetes without a correction")
  
  # drop if diabetes is missing -- we need this to properly do the correction
  data <- data[!is.na(diabetes),]
  
  # apply correction -- we do this even when we're not going to use it to make sure we drop the same data for all diabetes outcomes
  load("04_correction_models/output/diabetes_self_report_correction_model.rdata")
  formatted_data <- format_brfss(data)
  data[, diagnosed := formatted_data$diagnosed]
  for (ii in if (imputations == 0) 0 else 1:imputations) {
    data[, pred := correct_diabetes(formatted_data, fit_diabetes_undiag, fit_diabetes_diag)]
    if (grepl("^diabetes_diagnosed", outcome)) data[[paste0("outcome.", ii)]] <- data$diagnosed
    if (grepl("^diabetes_undiagnosed", outcome)) data[[paste0("outcome.", ii)]] <- as.numeric(data$diagnosed == 0 & data$pred == 1)
    if (grepl("^diabetes_controlled", outcome)) data[[paste0("outcome.", ii)]] <- as.numeric(data$diagnosed == 1 & data$pred == 0)  
    if (grepl("^diabetes_uncontrolled", outcome)) data[[paste0("outcome.", ii)]] <- as.numeric(data$diagnosed == 1 & data$pred == 1)
  }
  data <- data[!is.na(pred),]
  rm(formatted_data, format_brfss, fit_diabetes_diag, fit_diabetes_undiag, correct_diabetes); gc()
}

# cholesterol outcomes
if (grepl("^cholesterol", outcome)) {
  # make sure correction is set properly
  if (grepl("^cholesterol_diagnosed", outcome) & correction == 1) stop("*** No correction available for diagnosed cholesterol")
  if (grepl("^cholesterol_undiagnosed", outcome) & correction == 0) stop("*** Can't run undiagnosed cholesterol without a correction")
  if (grepl("^cholesterol_controlled", outcome) & correction == 0) stop("*** Can't run controlled cholesterol without a correction")
  if (grepl("^cholesterol_uncontrolled", outcome) & correction == 0) stop("*** Can't run uncontrolled cholesterol without a correction")
  
  # drop if cholesterol check is missing -- we need this to properly do the correction
  data <- data[!is.na(cholcheck),]
  
  # apply correction -- we do this even when we're not going to use it to make sure we drop the same data for all cholesterol outcomes
  load("04_correction_models/output/cholesterol_self_report_correction_model.rdata")
  formatted_data <- format_brfss(data)
  data[, diagnosed := formatted_data$diagnosed]
  for (ii in if (imputations == 0) 0 else 1:imputations) {
    data[, pred := correct_cholesterol(formatted_data, fit_cholesterol_undiag, fit_cholesterol_diag)]
    if (grepl("^cholesterol_diagnosed", outcome)) data[[paste0("outcome.", ii)]] <- data$diagnosed
    if (grepl("^cholesterol_undiagnosed", outcome)) data[[paste0("outcome.", ii)]] <- as.numeric(data$diagnosed == 0 & data$pred == 1)
    if (grepl("^cholesterol_controlled", outcome)) data[[paste0("outcome.", ii)]] <- as.numeric(data$diagnosed == 1 & data$pred == 0)  
    if (grepl("^cholesterol_uncontrolled", outcome)) data[[paste0("outcome.", ii)]] <- as.numeric(data$diagnosed == 1 & data$pred == 1)
  }
  data <- data[!is.na(pred),]
  rm(formatted_data, format_brfss, fit_cholesterol_diag, fit_cholesterol_undiag, correct_cholesterol); gc()
}

## Drop if outcome variable is missing  
for (var in grep("^outcome", names(data), value=T)) data <- data[!is.na(data[[var]]),]
  
## Generate age groups 
for (aa in 1:(length(age_groups)-1)) data[age >= age_groups[aa] & age < age_groups[aa+1], age := (aa-1)]
  
## Keep only necessary variables
setnames(data, "marital", "marriage")
levels(data$edu) <- gsub(" ", "_", levels(data$edu))
data <- data[, c("state", "fips", "year", "sample", "design_wt", "sex", "age", "race", "marriage", "edu", "phone", grep("outcome", names(data), value=T)), with=F]
setkeyv(data, c("fips", "year", "sex", "age", "race", "marriage", "edu", "phone"))
unlist(lapply(data, function(x) mean(is.na(x))))

## Calculate raw estimates and identify the gold standard 
if (!is.na(validation_type)) { 

# limit to validation years, collapse across imputations, and collapse across if performing a pooled validation
  raw <- data[!is.na(fips) & year %in% vyear1:vyear2 & !is.na(design_wt), ]
  if (imputations > 0) raw$outcome.0 <- apply(raw[, grep("outcome", names(raw), value=T), with=F], 1, mean, na.rm=T)
  raw <- na.omit(raw[, list(fips, year, design_wt, sex, age, race, marriage, edu, outcome.0)])
  if (validation_type == "pooled") raw[, year := median(vyear1:vyear2)]
 
# get counts for each county-year-sex and select the validation set
  raw[, fips_year_sex := paste(fips, year, sex, sep="_")]
  count <- table(raw$fips_year_sex)
  val.set <- names(count[count > k])
  raw <- raw[fips_year_sex %in% val.set,]
  
# get sex-age-race-marital-edu-specific estimates for each county   
  raw <- raw[, list(gs = weighted.mean(outcome.0, design_wt)), by='fips,year,fips_year_sex,sex,age,race,marriage,edu']
  
# weight the estimates by marital status, then education, then race to collapse to the sex-age level for each county 
  load(paste("05_small_area_models/inter/", outcome, "/pop_weights.rdata", sep=""))
  raw <- merge(raw, marriage_weights, by=c("fips", "year", "sex", "age", "marriage"), all.x=T)
  raw <- raw[, list(gs = weighted.mean(gs, wt, na.rm=T)), by='fips,year,fips_year_sex,sex,age,race,edu']
  raw <- merge(raw, edu_weights, by=c("fips", "year", "sex", "age", "edu"), all.x=T)
  raw <- raw[, list(gs = weighted.mean(gs, wt, na.rm=T)), by='fips,year,fips_year_sex,sex,age,race']
  raw <- merge(raw, race_weights, by=c("fips", "year", "sex", "age", "race"), all.x=T)
  raw <- raw[, list(gs = weighted.mean(gs, wt, na.rm=T)), by='fips,year,fips_year_sex,sex,age']
  
# age-standardize to collapse to the sex level for each county 
  load("02_covariates/inter/census_age_sex_standard.rdata")
  age.dis <- standard[standard$sex != "both" & standard$age >= age_groups[1] & standard$age < age_groups[length(age_groups)], c("sex", "age", "w")]
  for (aa in 1:(length(age_groups)-1)) age.dis$age[age.dis$age >= age_groups[aa] & age.dis$age < age_groups[aa+1]] <- aa - 1  
  age.dis <- aggregate(age.dis["w"], age.dis[c("sex", "age")], sum)
  age.dis$sex <- ifelse(age.dis$sex == "male", 1, 2) 
  raw <- merge(raw, age.dis, by=c("sex","age"))
  raw <- raw[, list(gs = weighted.mean(gs, w)), by='sex,fips,year,fips_year_sex']
  
# merge sample size count onto validation set information
  count <- data.frame(count)
  names(count) <- c("fips_year_sex", "count")
  raw <- merge(raw, count, by="fips_year_sex", all.x=T)
  val.set <- raw[, list(sex, fips, year, count, gs)]
  rm(raw, count); gc()

## Sample down county-years in the validation set 
  sampling <- sort(sampling[sampling != 9999])
  drop <- array(0, dim=c(length(sampling), iterations, nrow(data)))

# loop over county-sexes selected in the validation set, and then over year 
# NOTE: we sample down in all years that will be used in fitting validation models, not just validation years, to be more conservative: we don't want  
#       a county-year in the validation set to have an adjacent year with a large sample, as that will make prediction much easier
  all_cs <- unique(paste(val.set$fips, val.set$sex, sep="_"))
  for (cs in all_cs) { 
    cat(paste("\nvalidation county", which(all_cs == cs), "of", length(all_cs), "- ")); flush.console()
    fips <- as.numeric(strsplit(cs, split="_")[[1]][1])
    sex <- as.numeric(strsplit(cs, split="_")[[1]][2])
    
    for (year in year1:year2) { 
      cat(paste(year, " ")); flush.console()
    
# identify all observations to be sampled from 
      all <- which(data$fips == fips & data$sex == sex & data$year == year)

# loop over sampling levels and iterations and choose which observations to keep in each iteration 
      for (s in 1:length(sampling)) { 
        for (i in 1:iterations) { 
          drop[s, i, all[-sample(x=1:length(all), size=min(sampling[s], length(all)), replace=F)]] <- 1 # drop everything BUT the sampled observations
        } 
      }
    }
  }

# if we're not doing a validation, just create an empty object
} else { 
  val.set <- drop <- NA
} 
data[, design_wt := NULL]
  
## Save prepped data
save(data, val.set, drop, file=paste("05_small_area_models/inter/", outcome, "/model_data.rdata", sep=""))
