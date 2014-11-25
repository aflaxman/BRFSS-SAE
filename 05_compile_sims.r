

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
## Description: Compile model fits and simulations for all divisions and both sexes
##  (1) compile all sub-model fits and save
##  (2) compile simulations from all sub-models
##  (3) aggregate simulations to the state and national level
##  (4) age-standardize simulations
##  (5) apply the cell phone correction (if necessary)
##  (6) aggregate the simulations to both sexes combined
##  (7) collapse simulations to get CIs for levels, change, rate of change, and
##      ranks 
##  (8) calculate validation metrics (if running a validation) 
##  (9) save all output
########################################################################################################################

library(INLA)
library(data.table)

rm(list=ls())
set.seed(98121)
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
setwd(paste(root, "Project/NHIS/sae", sep=""))

## Get settings
outcome <- commandArgs()[3]
model <- as.numeric(commandArgs()[4])
race <- as.numeric(commandArgs()[5])
marriage <- as.numeric(commandArgs()[6])
edu <- as.numeric(commandArgs()[7])
imp <- as.numeric(commandArgs()[8])
samp <- as.numeric(commandArgs()[9])
iter <- as.numeric(commandArgs()[10])

try(detach(settings), silent=T)
get_settings <- read.csv("05_small_area_models/model_settings.csv", skip=1, stringsAsFactors=F)[,-1]
settings <- lapply(1:nrow(get_settings), function(x) strsplit(get_settings[x, outcome], split=";")[[1]])
settings <- lapply(settings, function(x) if(sum(is.na(as.numeric(x)))==0) as.numeric(x) else x)
names(settings) <- get_settings[,1]
attach(settings)
rm(get_settings, settings)  

# figure out if we need to get uncorrected estimates for the validation 
if (offset == 0) get_uncorrected <- 0 
if (offset == 1) if (sum(vyear1:vyear2 %in% offset_years) == 0) get_uncorrected <- 1 else get_uncorrected <- 0 


### Compile models -----------------------------------------------------------------------------------------------------

load(file=paste("05_small_area_models/inter/", outcome, "/models.rdata", sep=""))
if (correct_cells == 1 & samp != 9999) sub_models <- sub_models[sub_models$svy_sample == "l",]
rm(models); gc()

cat(paste("\nCompile model fits -", Sys.time())); flush.console() 
fit <- lapply(split(sub_models, 1:nrow(sub_models)), function(sm) { 
  load(paste("/clustertmp/counties/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu,  
             "/fit_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", sm$division, "_svy_sample_", sm$svy_sample, "_sex_", sm$sex, ".rdata", sep=""))
  fit <- fit[c(".args", "summary.fixed", "summary.random", "summary.hyperpar")]
  fit$.args$data <- NULL
  fit$.args$Ntrials <- NULL
  fit$.args$offset <- NULL
  fit
}) 
names(fit) <- sapply(split(sub_models, 1:nrow(sub_models)), function(sm) paste(c("sm_", sm), collapse=""))
save(fit, file=paste("05_small_area_models/inter/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu,  
                     "/fit_imp_", imp, "_samp_", samp, "_iter_", iter, ".rdata", sep=""))
rm(fit); gc()

### Process simulations ------------------------------------------------------------------------------------------------

## Compile simulations
cat(paste("\nCompile simulations -", Sys.time())); flush.console() 
sims <- lapply(split(sub_models, 1:nrow(sub_models)), function(sm) { 
  load(paste("/clustertmp/counties/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu,  
             "/sims_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", sm$division, "_svy_sample_", sm$svy_sample, "_sex_", sm$sex, ".rdata", sep=""))  
  sims[, svy_sample := sm$svy_sample]
  sims[, sex := sm$sex]
  sims
})
sims <- rbindlist(sims)

# pull out simulations for uncorrected estimates (only applicable if there is an internal correction and the validation period is in years not included in the offset)
if (get_uncorrected) { 
  uncorrected_sims <- sims[offset == 0,]
  uncorrected_sims[, offset := NULL]
  setkeyv(uncorrected_sims, c("svy_sample", "sex", "fips", "year", "age"))
  sims <- sims[offset == 1,]
}
sims[, offset := NULL] 
setkeyv(sims, c("svy_sample", "sex", "fips", "year", "age"))
simvars <- grep("^s[[:digit:]]", names(sims), value=T)

## Aggregate simulations to the state and national level (by age/sex/year/sample)
if (samp == 9999) { 
  cat(paste("\nAggregate simulations to state/national level -", Sys.time())); flush.console() 
  load(paste("05_small_area_models/inter/", outcome, "/pop_weights.rdata", sep="")) 
  setkeyv(fips_weights, c("fips", "year", "sex", "age"))
  wt <- fips_weights[sims[, list(fips, year, sex, age)],]$wt
  sims_state <- copy(sims)
  sims_state[, fips := floor(fips/1000)] 
  sims_state[, (simvars) := lapply(.SD, function(x) x * wt), .SDcols = simvars]
  sims_state[, (simvars) := lapply(.SD, function(x) sum(x)), by='fips,year,age,sex,svy_sample', .SDcols = simvars]
  setkeyv(sims_state, c("svy_sample", "sex", "fips", "year", "age"))
  sims_state <- unique(sims_state) 

  setnames(state_weights, "state", "fips")
  setkeyv(state_weights, c("fips", "year", "sex", "age"))
  wt <- state_weights[sims_state[, list(fips, year, sex, age)],]$wt
  sims_natl <- copy(sims_state)
  sims_natl[, fips := 0]
  sims_natl[, (simvars) := lapply(.SD, function(x) x * wt), .SDcols = simvars]
  sims_natl[, (simvars) := lapply(.SD, function(x) sum(x)), by='fips,year,age,sex,svy_sample', .SDcols = simvars]
  setkeyv(sims_natl, c("svy_sample", "sex", "fips", "year", "age"))
  sims_natl <- unique(sims_natl) 
  
  sims <- rbind(sims_natl, sims_state, sims, use.names=T)    
  rm(list=c(grep("_weights", ls(), value=T), "sims_state", "sims_natl", "wt")); gc()
} 

## Age-standardize to get county/state/national-level simulations
cat(paste("\nAge-standardize the simulations -", Sys.time())); flush.console() 
load("02_covariates/inter/census_age_sex_standard.rdata")
standard <- data.table(standard[standard$sex != "both" & standard$age >= age_groups[1] & standard$age < age_groups[length(age_groups)],])
standard[, sex := ifelse(sex == "male", 1, 2)]
standard[, age := as.numeric(as.character(cut(age, breaks=age_groups, labels=0:(length(age_groups)-2), right=F)))]
standard <- standard[, list(wt = sum(w)), by='age,sex']
standard[, wt := wt/sum(wt), by='sex']
setkeyv(standard, c("age", "sex"))

wt <- standard[sims[, list(age, sex)],]$wt
sims[, (simvars) := lapply(.SD, function(x) x * wt), .SDcols = simvars]
sims[, (simvars) := lapply(.SD, function(x) sum(x)), by='fips,year,sex,svy_sample', .SDcols = simvars]
sims[, age := NULL]
setkeyv(sims, c("svy_sample", "sex", "fips", "year"))
sims <- unique(sims) 

if (get_uncorrected) { 
  wt <- standard[uncorrected_sims[, list(age, sex)],]$wt
  uncorrected_sims[, (simvars) := lapply(.SD, function(x) x * wt), .SDcols = simvars]
  uncorrected_sims[, (simvars) := lapply(.SD, function(x) sum(x)), by='fips,year,sex,svy_sample', .SDcols = simvars]
  uncorrected_sims[, age := NULL]
  setkeyv(uncorrected_sims, c("svy_sample", "sex", "fips", "year"))
  uncorrected_sims <- unique(uncorrected_sims)
}  

## Apply cell phone correction (if necessary)
if (correct_cells == 1 & samp == 9999) { 
  cat(paste("\nApply the cell phone correction -", Sys.time())); flush.console() 
  cell_begin_year <- 2000
  sims_l <- sims[svy_sample == "l",]
  setkeyv(sims_l, c("fips", "sex", "year"))
  sims_b <- sims[svy_sample == "b",]
  setkeyv(sims_b, c("fips", "sex", "year"))

# calculate difference between estimates from combined sample and estimates from landline-only sample in 2011
  sims_c <- cbind(sims_b[year == 2011, list(fips, sex)], sims_b[year == 2011, simvars, with=F] - sims_l[year == 2011, simvars, with=F])
  
# linearly interpolate this difference back the last year in which the cell phone effect should be negligible, and then add it on to the landline only estimates
  sims_c <- sims_c[, lapply(simvars, function(ii) approx(x=c(cell_begin_year, 2011), y=c(0, .SD[[ii]]), xout=max(cell_begin_year,year1):2010)$y), by='fips,sex']
  setnames(sims_c, paste("V", 1:length(simvars), sep=""), simvars)
  sims_c <- cbind(sims_l[year >= cell_begin_year & year <= 2010, c("fips", "sex", "year"), with=F], 
                  sims_l[year >= cell_begin_year & year <= 2010, simvars, with=F] + sims_c[, simvars, with=F])
  sims_c$svy_sample <- "b" 
  
# in a small number of cases, this correction will make estimates for some simulations negative; replace these with the smallest observed non-negative value among the other predictions for that sim
  ng <- sims_c[, lapply(simvars, function(x) sum(.SD[[x]] < 0))] 
  for (ii in which(ng > 0)) sims_c[[simvars[ii]]] <- ifelse(sims_c[[simvars[ii]]] < 0, min(abs(sims_c[[simvars[ii]]])), sims_c[[simvars[ii]]])
  
# combine predictions back together
  sims_l2 <- sims[svy_sample == "l" & year < cell_begin_year,]
  sims_l2[, svy_sample := "b"]
  sims <- rbind(sims_l, sims_l2, sims_c, sims_b, use.names=T)
  setkeyv(sims, c("fips", "year", "sex", "svy_sample"))      
  rm(sims_l, sims_l2, sims_b, sims_c, ng, ii); gc()
} 

## Aggregate simulations to get both sexes combined
if (samp == 9999) { 
  cat(paste("\nCombine sex-specific simulations to get simulations for both sexes combined -", Sys.time())); flush.console() 
  wt <- sex_standard[sims$sex]
  sims_both <- copy(sims)
  sims_both[, (simvars) := lapply(.SD, function(x) x * wt), .SDcols = simvars]
  sims_both[, (simvars) := lapply(.SD, function(x) sum(x)), by='fips,year,svy_sample', .SDcols = simvars]
  sims_both <- sims_both[sex == 1,]
  sims_both[, sex := 3]
  sims <- rbind(sims, sims_both)
  rm(sims_both); gc()
}

## Save compiled simulations
setkeyv(sims, c("fips", "year", "sex", "svy_sample"))
save(sims, file=paste("/clustertmp/counties/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, "/sims_imp_", imp, "_samp_", samp, "_iter_", iter, "_compiled.rdata", sep=""))

### Calculate CI -------------------------------------------------------------------------------------------------------

## Write function for getting all summary statistics from the simulations
sumsims <- function(x, sig_change=F) {
  if (sig_change) output <- c(mean(x), quantile(x, c(0.025, 0.975)), var(x), as.numeric(mean(x<0) < 0.05 | mean(x<0) > 0.95))
  if (!sig_change) output <- c(mean(x), quantile(x, c(0.025, 0.975)), var(x))
  lapply(as.numeric(output), function(y) y)
}  

## Calculate confidence intervals for yearly estimates
cat(paste("\nCalculate CI for yearly estimates -", Sys.time())); flush.console() 
est <- sims[, sumsims(unlist(.SD)), by='fips,sex,year,svy_sample']
setnames(est, paste("V", 1:4, sep=""), c("prev", "prev.lb", "prev.ub", "prev.var"))  
setkeyv(est, c("fips", "year", "sex", "svy_sample"))

if (get_uncorrected) {   
  uncorrected_est <- uncorrected_sims[, sumsims(unlist(.SD)), by='fips,sex,year,svy_sample']
  setnames(uncorrected_est, paste("V", 1:4, sep=""), c("prev", "prev.lb", "prev.ub", "prev.var"))  
  setkeyv(uncorrected_est, c("fips", "year", "sex", "svy_sample"))
} else { 
  uncorrected_est <- NULL 
}

## Calculate yearly estimate ranks + confidence intervals
if (samp == 9999) {
  cat(paste("\nCalculate CI for yearly estimate ranks -", Sys.time())); flush.console()
  rank_est <- sims[fips > 100,]
  rank_est[, (simvars) := lapply(.SD, function(x) rank(x)), by='sex,year,svy_sample', .SDcols = simvars]
  rank_est <- rank_est[, sumsims(unlist(.SD)), by='fips,sex,year,svy_sample']
  setnames(rank_est, paste("V", 1:4, sep=""), c("rank", "rank.lb", "rank.ub", "rank.var"))  
  setkeyv(rank_est, c("fips", "year", "sex", "svy_sample"))
  rank_est <- est[rank_est,]
  rank_est[, rank := rank(prev), by='sex,year,svy_sample']
  rank_est[, rank.lb := pmin(rank, floor(rank.lb))]
  rank_est[, rank.ub := pmax(rank, ceiling(rank.ub))]
  rank_est <- rank_est[, list(fips, year, sex, svy_sample, rank, rank.lb, rank.ub, rank.var)]
  setkeyv(rank_est, c("fips", "year", "sex", "svy_sample"))
} else { 
  rank_est <- NULL 
}  
  
## Calculate absolute changes + confidence intervals 
if (samp == 9999) {
  cat(paste("\nCalculate CI for differences -", Sys.time())); flush.console()
  diff <- lapply(diff_years, function(x) { 
                   cat(paste("\n", x)); flush.console()
                   y1 <- as.numeric(substr(x, 1, 4))
                   y2 <- as.numeric(substr(x, 6, 9))
                   diff <- cbind(sims[year == y2, list(fips, sex, svy_sample)], sims[year == y2, simvars, with=F] - sims[year == y1, simvars, with=F]) 
                   diff <- diff[, sumsims(unlist(.SD), T), by='fips,sex,svy_sample']
                   setnames(diff, paste("V", 1:5, sep=""), c("diff", "diff.lb", "diff.ub", "diff.var", "ot_sig_change"))
                   diff[, years := x]
                 })
  diff <- rbindlist(diff)
  setkeyv(diff, c("fips", "years", "sex", "svy_sample"))
} else { 
  diff <- NULL 
}

## Calculate absolute change ranks + confidence intervals
if (samp == 9999) {
  cat(paste("\nCalculate CI for differences rank -", Sys.time())); flush.console()
  rank_diff <- lapply(diff_years, function(x) { 
                   cat(paste("\n", x)); flush.console()
                   y1 <- as.numeric(substr(x, 1, 4))
                   y2 <- as.numeric(substr(x, 6, 9))
                   rank_diff <- cbind(sims[year == y2, list(fips, sex, svy_sample)], sims[year == y2, simvars, with=F] - sims[year == y1, simvars, with=F]) 
                   rank_diff <- rank_diff[fips > 100,]
                   rank_diff[, (simvars) := lapply(.SD, function(x) rank(x)), by='sex,svy_sample', .SDcols = simvars]
                   rank_diff <- rank_diff[, sumsims(unlist(.SD), F), by='fips,sex,svy_sample']
                   setnames(rank_diff, paste("V", 1:4, sep=""), c("rank", "rank.lb", "rank.ub", "rank.var"))
                   rank_diff[, years := x]
                   setkeyv(rank_diff, c("fips", "years", "sex", "svy_sample"))
                   rank_diff <- diff[rank_diff,]
                   rank_diff[, rank := rank(diff), by='sex,years,svy_sample']
                   rank_diff[, rank.lb := pmin(rank, floor(rank.lb))]
                   rank_diff[, rank.ub := pmax(rank, ceiling(rank.ub))]
                   rank_diff <- rank_diff[, list(fips, years, sex, svy_sample, rank, rank.lb, rank.ub, rank.var)]                     
                 })
  rank_diff <- rbindlist(rank_diff)
  setkeyv(rank_diff, c("fips", "years", "sex", "svy_sample"))
} else { 
  rank_diff <- NULL 
}

## Calculate percent change + confidence intervals 
if (samp == 9999) {
  cat(paste("\nCalculate CI for percent change -", Sys.time())); flush.console()
  pct_diff <- lapply(diff_years, function(x) { 
    cat(paste("\n", x)); flush.console()
    y1 <- as.numeric(substr(x, 1, 4))
    y2 <- as.numeric(substr(x, 6, 9))
    pct_diff <- cbind(sims[year == y2, list(fips, sex, svy_sample)], (sims[year == y2, simvars, with=F] - sims[year == y1, simvars, with=F])/(sims[year == y1, simvars, with=F])) 
    pct_diff <- pct_diff[, sumsims(unlist(.SD), T), by='fips,sex,svy_sample']
    setnames(pct_diff, paste("V", 1:5, sep=""), c("pct_diff", "pct_diff.lb", "pct_diff.ub", "pct_diff.var", "ot_sig_change"))
    setkeyv(pct_diff, c("fips", "sex", "svy_sample"))
    pct_diff[, pct_diff := (est$prev[est$year == y2] - est$prev[est$year == y1])/(est$prev[est$year == y1])]
    pct_diff[, years := x]
  })
  pct_diff <- rbindlist(pct_diff)
  setkeyv(pct_diff, c("fips", "years", "sex", "svy_sample"))
} else { 
  pct_diff <- NULL 
}

## Calculate percent change ranks + confidence intervals 
if (samp == 9999) {
  cat(paste("\nCalculate CI for percent change rank -", Sys.time())); flush.console()
  rank_pct_diff <- lapply(diff_years, function(x) { 
    cat(paste("\n", x)); flush.console()
    y1 <- as.numeric(substr(x, 1, 4))
    y2 <- as.numeric(substr(x, 6, 9))
    rank_pct_diff <- cbind(sims[year == y2, list(fips, sex, svy_sample)], (sims[year == y2, simvars, with=F] - sims[year == y1, simvars, with=F])/(sims[year == y1, simvars, with=F])) 
    rank_pct_diff <- rank_pct_diff[fips > 100,]
    rank_pct_diff[, (simvars) := lapply(.SD, function(x) rank(x)), by='sex,svy_sample', .SDcols = simvars]
    rank_pct_diff <- rank_pct_diff[, sumsims(unlist(.SD), F), by='fips,sex,svy_sample']
    setnames(rank_pct_diff, paste("V", 1:4, sep=""), c("rank", "rank.lb", "rank.ub", "rank.var"))
    rank_pct_diff[, years := x]
    setkeyv(rank_pct_diff, c("fips", "years", "sex", "svy_sample"))
    rank_pct_diff <- pct_diff[rank_pct_diff,]
    rank_pct_diff[, rank := rank(pct_diff), by='sex,years,svy_sample']
    rank_pct_diff[, rank.lb := pmin(rank, floor(rank.lb))]
    rank_pct_diff[, rank.ub := pmax(rank, ceiling(rank.ub))]
    rank_pct_diff <- rank_pct_diff[, list(fips, years, sex, svy_sample, rank, rank.lb, rank.ub, rank.var)]                          
  })
  rank_pct_diff <- rbindlist(rank_pct_diff)
  setkeyv(rank_pct_diff, c("fips", "years", "sex", "svy_sample"))
} else { 
  rank_pct_diff <- NULL 
}

## Calculate annualized rates of change + confidence intervals 
if (samp == 9999) {
  cat(paste("\nCalculate CI for annualized rates of change -", Sys.time())); flush.console()
  arc <- lapply(diff_years, function(x) { 
                   cat(paste("\n", x)); flush.console()
                   y1 <- as.numeric(substr(x, 1, 4))
                   y2 <- as.numeric(substr(x, 6, 9))
                   arc <- cbind(sims[year == y2, list(fips, sex, svy_sample)], log(sims[year == y2, simvars, with=F]/sims[year == y1, simvars, with=F])/(y2-y1)) 
                   arc <- arc[, sumsims(unlist(.SD), T), by='fips,sex,svy_sample']
                   setnames(arc, paste("V", 1:5, sep=""), c("arc", "arc.lb", "arc.ub", "arc.var", "ot_sig_change"))
                   setkeyv(arc, c("fips", "sex", "svy_sample"))
                   arc[, arc := log(est$prev[est$year == y2]/est$prev[est$year == y1])/(y2-y1)]
                   arc[, years := x]
                 })
  arc <- rbindlist(arc)
  setkeyv(arc, c("fips", "years", "sex", "svy_sample"))
} else { 
  arc <- NULL 
}

## Calculate annualized rates of change ranks + confidence intervals 
if (samp == 9999) {
  cat(paste("\nCalculate CI for annualized rates of change rank -", Sys.time())); flush.console()
  rank_arc <- lapply(diff_years, function(x) { 
                       cat(paste("\n", x)); flush.console()
                       y1 <- as.numeric(substr(x, 1, 4))
                       y2 <- as.numeric(substr(x, 6, 9))
                       rank_arc <- cbind(sims[year == y2, list(fips, sex, svy_sample)], log(sims[year == y2, simvars, with=F]/sims[year == y1, simvars, with=F])/(y2-y1)) 
                       rank_arc <- rank_arc[fips > 100,]
                       rank_arc[, (simvars) := lapply(.SD, function(x) rank(x)), by='sex,svy_sample', .SDcols = simvars]
                       rank_arc <- rank_arc[, sumsims(unlist(.SD), F), by='fips,sex,svy_sample']
                       setnames(rank_arc, paste("V", 1:4, sep=""), c("rank", "rank.lb", "rank.ub", "rank.var"))
                       rank_arc[, years := x]
                       setkeyv(rank_arc, c("fips", "years", "sex", "svy_sample"))
                       rank_arc <- arc[rank_arc,]
                       rank_arc[, rank := rank(arc), by='sex,years,svy_sample']
                       rank_arc[, rank.lb := pmin(rank, floor(rank.lb))]
                       rank_arc[, rank.ub := pmax(rank, ceiling(rank.ub))]
                       rank_arc <- rank_arc[, list(fips, years, sex, svy_sample, rank, rank.lb, rank.ub, rank.var)]                          
                     })
  rank_arc <- rbindlist(rank_arc)
  setkeyv(rank_arc, c("fips", "years", "sex", "svy_sample"))
} else { 
  rank_arc <- NULL 
}
  
### Calculate error and coverage ---------------------------------------------------------------------------------------

if (!is.na(validation_type)) {
  load(paste("05_small_area_models/inter/", outcome, "/model_data.rdata", sep=""))
  val.set <- data.frame(val.set)[, c("fips", "sex", "year", "gs")]

  if (get_uncorrected) error <- data.frame(uncorrected_est) else error <- data.frame(est)
  if (correct_cells == 1) error <- error[(error$svy_sample == "l" & error$year <= 2010) | (error$svy_sample == "b" & error$year > 2010), c("fips", "sex", "year", "prev", "prev.lb", "prev.ub")]
  else error <- error[, c("fips", "sex", "year", "prev", "prev.lb", "prev.ub")]
  error <- merge(val.set, error, by=c("fips", "sex", "year"), all.x=T)
  
  val_output <- lapply(1:2, function(sex) { 
                         error <- error[error$sex == sex,]
                         rmse <- sqrt(mean((error$gs - error$prev)^2, na.rm=T))
                         cor <- (2*cor(error$gs, error$prev, use="pairwise.complete.obs")*sd(error$prev, na.rm=T)*sd(error$gs, na.rm=T))/(var(error$gs, na.rm=T) + var(error$prev, na.rm=T) + (mean(error$gs, na.rm=T) - mean(error$prev, na.rm=T))^2)
                         mre <- median(abs(error$gs - error$prev)/error$gs, na.rm=T)
                         coverage <- mean(error$prev.lb <= error$gs & error$gs <= error$prev.ub, na.rm=T)
                         c(sex=sex, rmse=rmse, cor=cor, mre=mre, coverage=coverage)
                       }) 
  val_output <- data.frame(do.call("rbind", val_output))
} else { 
  val_output <- NA
}  

## Save all output
save(est, uncorrected_est, diff, pct_diff, arc, rank_est, rank_diff, rank_pct_diff, rank_arc, val_output, file=paste("05_small_area_models/inter/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, "/preds_imp_", imp, "_samp_", samp, "_iter_", iter, ".rdata", sep=""))
