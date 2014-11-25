

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
## Description: Prep population weights for collapsing across education status, 
##              marital status, race, fips, and states. 
########################################################################################################################

library(data.table)
library(reshape)

rm(list=ls())
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
thisdate <- format(Sys.time(), format = "%Y_%m_%d")
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

## Write function for calculating weights but omitting anything that is overly small so we don't need to carry so many digits 
round_weights <- function(x, min=0.00001) { 
  x <- x/sum(x) 
  x[x < min] <- 0 
  x <- x/sum(x) 
  x <- round(x, nchar(1/min)+1)
  return(x)
} 

## Load population files (by race, by marital status, by education status) 
load("02_covariates/output/population.rdata")

## Calculate weights for collapsing from education-fips to fips (within age/sex) 
edu_pop <- collapse_pop(edu_pop, "edu", c(year1, year2), age_groups)
edu_weights <- edu_pop[, list(edu = edu, wt = round_weights(pop)), by='fips,year,age,sex']
temp <- edu_pop[, list(pop = sum(pop)), by='fips,year,sex,edu'][, list(edu = edu, wt2 = round_weights(pop)), by='fips,year,sex'] # this is the same as above but pools across age; we'll use this when there's no population in a given age group to calculate weights with
edu_weights <- merge(edu_weights, temp, by=c("fips", "year", "sex", "edu"))
edu_weights <- edu_weights[, list(fips, year, sex, edu = gsub(" ", "_", edu), age = as.numeric(as.factor(age))-1, wt = ifelse(is.na(wt), wt2, wt))]
rm(temp, edu_pop); gc()

## Calculate weights for collapsing from marital-fips to fips (within age/sex) 
marital_pop <- collapse_pop(marital_pop, "marital", c(year1, year2), age_groups)
setnames(marital_pop, "marital", "marriage")
marriage_weights <- marital_pop[, list(marriage = marriage, wt = round_weights(pop)), by='fips,year,age,sex']
temp <- marital_pop[, list(pop = sum(pop)), by='fips,year,sex,marriage'][, list(marriage = marriage, wt2 = round_weights(pop)), by='fips,year,sex'] # this is the same as above but pools across age; we'll use this when there's no population in a given age group to calculate weights with
marriage_weights <- merge(marriage_weights, temp, by=c("fips", "year", "sex", "marriage"))
marriage_weights <- marriage_weights[, list(fips, year, sex, marriage, age = as.numeric(as.factor(age))-1, wt = ifelse(is.na(wt), wt2, wt))]  
rm(temp, marital_pop); gc() 
  
## Calculate weights for collapsing from race-fips to fips (within age/sex) 
race_pop$race <- factor(race_pop$race, levels=c("white", "black", "hisp", "asian", "native"), labels=c("white", "black", "hisp", "other", "native"))
race_pop <- collapse_pop(race_pop, "race", c(year1, year2), age_groups)
race_weights <- race_pop[, list(race = race, wt = round_weights(pop)), by='fips,year,age,sex']
temp <- race_pop[, list(pop = sum(pop)), by='fips,year,sex,race'][, list(race = race, wt2 = round_weights(pop)), by='fips,year,sex'] # this is the same as above but pools across age; we'll use this when there's no population in a given age group to calculate weights with
race_weights <- merge(race_weights, temp, by=c("fips", "year", "sex", "race"))
race_weights <- race_weights[, list(fips, year, sex, race, age = as.numeric(as.factor(age))-1, wt = ifelse(is.na(wt), wt2, wt))]
rm(temp); gc()
  
## Calculate weights for collapsing from fips to state (within age/sex)
race_pop$age <- as.numeric(as.factor(race_pop$age))-1
race_pop <- race_pop[, list(state = floor(fips/1000), pop = sum(pop)), by='fips,year,age,sex']
fips_weights <- race_pop[, list(fips = fips, wt = round_weights(pop)), by='state,year,age,sex']

## Calculate weights for collapsing from state to national (within age/sex)   
race_pop <- race_pop[, list(pop = sum(pop)), by='state,year,age,sex'] 
state_weights <- race_pop[, list(state = state, wt = round_weights(pop)), by='year,age,sex']  
rm(race_pop); gc()

## Calculate weights for collapsing over phone usage categories 
# load phone usage estimates, remove no phone group
load("02_covariates/inter/CDC_NHSR_phone_usage.rdata")
data <- data.table(melt(data[, c("m.fips", "year", "landline_only", "cell_only", "dual_phone")], id.var=c("m.fips", "year")))
setnames(data, c("m.fips", "variable"), c("fips", "phone"))

# calculate weights
phone_weights <- data[, list(phone = phone, wt = round_weights(value)), by='fips,year']
rm(data); gc()

## Save all sets of weights
for (var in c("fips", "year", "sex", "age")) edu_weights[[var]] <- as.integer(edu_weights[[var]])
setkeyv(edu_weights, c("fips", "year", "sex", "age", "edu"))
for (var in c("fips", "year", "sex", "age")) marriage_weights[[var]] <- as.integer(marriage_weights[[var]])
setkeyv(marriage_weights, c("fips", "year", "sex", "age", "marriage"))
for (var in c("fips", "year", "sex", "age")) race_weights[[var]] <- as.integer(race_weights[[var]])
setkeyv(race_weights, c("fips", "year", "sex", "age", "race"))
for (var in c("fips", "state", "year", "sex", "age")) fips_weights[[var]] <- as.integer(fips_weights[[var]])
setkeyv(fips_weights, c("state", "fips", "year", "sex", "age"))
for (var in c("state", "year", "sex", "age")) state_weights[[var]] <- as.integer(state_weights[[var]])
setkeyv(state_weights, c("state", "year", "sex", "age"))
phone_weights$fips <- as.integer(phone_weights$fips)
setkeyv(phone_weights, c("fips", "year", "phone"))

save(edu_weights, marriage_weights, race_weights, fips_weights, state_weights, phone_weights, file=paste("05_small_area_models/inter/", outcome, "/pop_weights.rdata", sep=""))
