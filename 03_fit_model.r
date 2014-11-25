

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
## Description: Fit the division/sex-specific small area models
##              (1) Prep the model data
##              (2) Fit the initial state-level model to get empirical priors (if this option is turned on) 
##              (3) Prep the dataset needed for prediction from the county-level model 
##              (4) Fit the county-level model and save the model fit
########################################################################################################################

library(INLA)
library(data.table)
library(dummies)

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
svy_sample <- commandArgs()[11]
this_sex <- sex <- as.numeric(commandArgs()[12])
division <- as.numeric(commandArgs()[13])

try(detach(settings), silent=T)
get_settings <- read.csv("05_small_area_models/model_settings.csv", skip=1, stringsAsFactors=F)[,-1]
settings <- lapply(1:nrow(get_settings), function(x) strsplit(get_settings[x, outcome], split=";")[[1]])
settings <- lapply(settings, function(x) if(sum(is.na(as.numeric(x)))==0) as.numeric(x) else x)
names(settings) <- get_settings[,1]
attach(settings)
rm(get_settings, settings)

phone <- as.numeric(correct_cells == 1 & svy_sample == "b")

divisions <- read.csv("01_counties/census_regions_divisions.csv")[, c("state_fips", "division2_fips")]
include_states <- divisions$state_fips[divisions$division2_fips == division]
rm(divisions)

### Prep Data ----------------------------------------------------------------------------------------------------------

# load data and choose correct outcome variable
load(paste("05_small_area_models/inter/", outcome, "/model_data.rdata", sep=""))
setnames(data, paste("outcome", imp, sep="."), "outcome")

# sample down, if necessary
if (samp != 9999) data <- data[drop[which(sampling[sampling!=9999]==samp), iter, ] == 0,]
 
# keep only the sex and division we're currently fitting
data <- data[sex == this_sex & state %in% include_states, list(state, fips, year, age, race, marriage, edu, phone, outcome)]
if (!is.na(validation_type)) val.set <- val.set[floor(fips/1000) %in% include_states & sex == this_sex, list(fips, year)]
 
# choose specified survey sample 
if (svy_sample == "l") data <- data[year <= 2010 | phone %in% c("landline_only", "dual_phone"),]
if (svy_sample == "b" & correct_cells == 1) data <- data[year >= 2011,] 

# dummy code the demographic variables 
if (race == 1) data <- cbind(data, dummy("race", data, sep="_")[,-1])
data[, race := NULL]
if (edu == 1) data <- cbind(data, dummy("edu", data, sep="_")[,-1])
data[, edu := NULL]
if (marriage == 1) data <- cbind(data, dummy("marriage", data, sep="_")[,-1])
data[, marriage := NULL]
if (phone == 1) data <- cbind(data, dummy("phone", data, sep="_")[,-1])
data[, phone := NULL]

# identify years for offset
if (offset == 1) {
  data[, offset := as.numeric(year %in% offset_years)]
} else { 
  data[, offset := 1]
} 

### Fit the initial state-level model ----------------------------------------------------------------------------------
if (initial_state_model) { 

# collapse data to the state level 
  state_data <- copy(data)
  setkeyv(state_data, c("state", "year", "age"))
  state_data <- state_data[, c(list(Y = sum(outcome), N = length(outcome)), lapply(.SD[, -1, with=F], mean)), by='state,year,age,offset']

# create variables for random effects
  setkeyv(state_data, c("state", "year", "age", "offset"))
  state_data$state_iid <- as.numeric(as.factor(state_data$state))
  state_data$time_rw <- state_data$time_ln <- as.numeric(as.factor(state_data$year))
  state_data$int <- as.numeric(factor(paste(state_data$state, state_data$year, sep="_"), levels=unique(paste(state_data$state, state_data$year, sep="_"))))

# construct the model formula 
  form <- "Y ~ 1 + factor(age)" 
  if (offset == 1 & length(unique(data$offset)) > 1) form <- paste(form, "+ offset")
  for (var in c("race", "marriage", "edu", "phone")) { 
    if (get(var) == 0) next
    form <- paste(form, "+", paste(grep(var, names(data), value=T), collapse=" + "))
    if (var %in% interact_age) form <- paste(form, "+", paste(paste(grep(var, names(data), value=T), ":factor(age)", sep=""), collapse=" + "))
  }
  form <- paste(form, "+ time_ln + f(time_rw, model=\"rw1\") + f(state_iid, model=\"iid\") + f(int, model=\"iid\")") 
  form

# fit the model 
  fit <- inla(as.formula(form), data=state_data, family="binomial", Ntrials=N, 
              control.inla=list(strategy="gaussian"),
              control.fixed=list(prec.intercept=0.001),
              num.threads=6)

# extract fixed effects predictions to use as priors 
  prior_mean <- c(fit$summary.fixed[, "mean"], default=0)
  prior_precision <- c(1/fit$summary.fixed[, "sd"]^2, default=0.001)

} else {
  prior_mean <- 0
  prior_precision <- 0.001
  
}

### Prep for prediction from county-level model ------------------------------------------------------------------------

# collapse down the data
data <- data[!is.na(fips),]
setkeyv(data, c("fips", "year", "age", "offset"))
data <- data[, c(list(Y = sum(outcome), N = length(outcome)), lapply(.SD[, -1, with=F], mean)), by='fips,year,age,offset']

# add line for predictions based on population fractions for all demographic variables in the model 
load(paste("05_small_area_models/inter/", outcome, "/pop_weights.rdata", sep=""))
race_weights <- reshape(race_weights[sex == this_sex & floor(fips/1000) %in% include_states & year %in% min(data$year):max(data$year), list(fips, year, race, age, wt)], direction="wide", idvar=c("fips", "year", "age"), timevar="race", drop="sex")
setattr(race_weights, "names", gsub("wt.", "race_", names(race_weights)))
edu_weights <- reshape(edu_weights[sex == this_sex & floor(fips/1000) %in% include_states & year %in% min(data$year):max(data$year), list(fips, year, edu = gsub(" ", "_", edu), age, wt)], direction="wide", idvar=c("fips", "year", "age"), timevar="edu", drop="sex")
setattr(edu_weights, "names", gsub("wt.", "edu_", names(edu_weights)))
marriage_weights <- reshape(marriage_weights[sex == this_sex & floor(fips/1000) %in% include_states & year %in% min(data$year):max(data$year), list(fips, year, marriage, age, wt)], direction="wide", idvar=c("fips", "year", "age"), timevar="marriage", drop="sex")
setattr(marriage_weights, "names", gsub("wt.", "marriage_", names(marriage_weights)))
pred_square <- Reduce(function(x, y) merge(x, y, by=c("fips", "year", "age")), list(race_weights, edu_weights, marriage_weights))
if (phone == 1) { 
  phone_weights <- reshape(phone_weights[sex == this_sex & floor(fips/1000) %in% include_states & year %in% min(data$year):max(data$year), list(fips, year, phone, wt)], direction="wide", idvar=c("fips", "year"), timevar="phone")
  setattr(phone_weights, "names", gsub("wt.", "phone_", names(phone_weights)))
  pred_square <- merge(pred_square, phone_weights, by=c("fips", "year"))
} 
pred_square <- pred_square[, names(pred_square) %in% names(data), with=F]

# if running a validation model, subset to only those county-years which we need to predict for
if (samp != 9999) {
  miss <- merge(unique(pred_square[, list(fips, year)]), data[, list(N = sum(N)), by='fips,year'], by=c("fips", "year"), all=T)
  miss <- miss[is.na(N), list(fips, year)] # these are county-years that aren't already in the dataset -- we need to add these to estimate the interaction properly 
  pred_square <- pred_square[rbind(miss, val.set),] # keep county-years in the validation set (so we can get predictions) and county-years not already in the dataset (so we can estimate the interaction properly)
}

# if we need non-corrected estimates for the validation of a model with an offset, add in these lines
pred_square[, offset := 1]
if (offset == 1 & !is.na(validation_type)) {
  if (sum(vyear1:vyear2 %in% offset_years) == 0 & sum(!val.set$year %in% pred_square$year) == 0) { # if validation years are not offset years we need to predict uncorrected estimates for the validation set to use for calculating validation metrics
    setkeyv(pred_square, c("fips", "year", "age"))
    temp <- pred_square[val.set,] 
    temp[, offset := 0]
    pred_square <- rbind(pred_square, temp)
    rm(temp)
  }
} 
  
# combine the real data with the prediction data
data$type <- "data"
pred_square$type <- "pred"
data <- merge(data, pred_square, by=names(pred_square), all=T)
rm(pred_square, race_weights, edu_weights, marriage_weights, phone_weights, fips_weights, state_weights, drop, val.set, miss); gc()

### Fit the county-level model -----------------------------------------------------------------------------------------

# merge on covariates, if necessary 
if (model %in% 3:4) { 
  load("02_covariates/output/covariates.rdata")
  covar <- covar[[sex]]
  for (var in names(covar)[!names(covar) %in% c("fips", "year")]) covar[,var] <- (covar[,var] - mean(covar[,var]))/sd(covar[,var])
  covar <- covar[covar$year %in% min(data$year):max(data$year), c("fips", "year", area_vars)]
  data <- merge(data, covar, by=c("fips", "year"), all.x=T)
  rm(covar)
} 

# create variables for random effects
setkeyv(data, c("fips", "year", "age", "type", "offset"))
data$cnty_iid <- data$cnty_icar <- as.numeric(as.factor(data$fips))
data$time_rw <- data$time_ln <- as.numeric(as.factor(data$year))
data$int <- as.numeric(factor(paste(data$fips, data$year, sep="_"), levels=unique(paste(data$fips, data$year, sep="_"))))

# prep the structure matrix for the space effects (IID or ICAR) 
load("01_counties/output/county_adjacencies.rdata")  
neighbors <- neighbors[neighbors$fips %in% data$fips & neighbors$neighbor %in% data$fips,]
all_fips <- sort(unique(neighbors$fips))
if (model %in% c(2,4)) { # ICAR
  graph <- matrix(0, nrow=length(all_fips), ncol=length(all_fips))
  for (ff in 1:length(all_fips)) graph[ff, match(neighbors$neighbor[neighbors$fips == all_fips[ff]], all_fips)] <- 1
  S_cnty <- diag(apply(graph, 1, sum)) - graph
} else { # IID
  S_cnty <- diag(rep(1, length(all_fips))) 
}

# prep the structure matrix for the time effects (RW)    
t <- length(unique(data$year))    
S_time <- matrix(0, nrow=t, ncol=t)
for (ii in 1:t) {
  S_time[ii,ii] <- ifelse(ii == 1 | ii == t, 1, 2)
  if (ii < t) {
    S_time[ii,ii+1] <- -1
    S_time[ii+1,ii] <- -1
  }
} 

# prep the structure matrix for the interaction (type II or IV depending on space effect)  
S_int <- kronecker(S_cnty, S_time) 

# set up the linear constraints for the interaction (note: it takes too long to get eigenvectors from the interaction matrix directly, but we can use properties of kronecker products to get there faster (see theorem 13.2 in http://www.siam.org/books/textbooks/OT91sample.pdf)
e1 <- eigen(S_cnty, symmetric=T)
e2 <- eigen(S_time, symmetric=T)
vectors <- kronecker(e1$vectors, e2$vectors)
values <- kronecker(e1$values, e2$values)
A <- t(vectors[, which(abs(values)<1e-10)])
e <- rep(0, nrow(A))

# construct the model formula 
if (outcome == "binge_drinking" & sex == 1) offset <- 0 # we don't want to include the binge-drinking offset for males, since the definition of binge drinking changed for females only
if (outcome == "excess_drinking" & sex == 1) offset <- 0 # same as for binge-drinking

form <- "Y ~ 1 + factor(age)" 
if (offset == 1 & length(unique(data$offset)) > 1) form <- paste(form, "+ offset")
for (var in c("race", "marriage", "edu", "phone")) { 
  if (get(var) == 0) next
  form <- paste(form, "+", paste(grep(var, names(data), value=T), collapse=" + "))
  if (var %in% interact_age) form <- paste(form, "+", paste(paste(grep(var, names(data), value=T), ":factor(age)", sep=""), collapse=" + "))
}
if (model %in% 3:4) form <- paste(form, "+", paste(area_vars, collapse=" + "))
if (model %in% c(1, 3)) form <- paste(form, "+ f(cnty_iid, model=\"iid\", hyper=list(prec=list(param=c(1,0.01))))")
if (model %in% c(2, 4)) form <- paste(form, "+ f(cnty_iid, model=\"iid\", hyper=list(prec=list(param=c(1,0.01)))) + f(cnty_icar, model=\"besag\", graph=graph, hyper=list(prec=list(param=c(1,0.01))), adjust.for.con.comp=F)")  
form <- paste(form, "+ time_ln + f(time_rw, model=\"rw1\", hyper=list(prec=list(param=c(1,0.01)))) + f(int, model=\"generic0\", Cmatrix=S_int, extraconstr=list(A=A,e=e), hyper=list(prec=list(param=c(1,0.01))), diagonal=1e-4)") 
form
  
# fit the model 
cat(paste("***************** Fit model:", Sys.time(), " *****************\n")); flush.console()                                           
fit <- try(inla(as.formula(form), data=data, family="binomial", Ntrials=N, 
                control.inla=list(strategy="gaussian"),
                control.fixed=list(prec.intercept=0.001, mean=prior_mean, prec=prior_precision),
                control.compute=list(config=TRUE), 
                num.threads=6, verbose=T))

if (class(fit) == "try-error" | length(fit) < 10) { 
  cat(paste("***************** Direct fit failed; fit initial model:", Sys.time(), " *****************\n")); flush.console()  
  initial <- inla(as.formula(form), data=data, family="binomial", Ntrials=N,
                  control.inla=list(strategy="gaussian", int.strategy="eb", diagonal=100),
                  control.fixed=list(prec.intercept=0.001, mean=prior_mean, prec=prior_precision),
                  num.threads=6, verbose=T)  
  
  cat(paste("***************** Initial model succeeded; fit final model:", Sys.time(), " *****************\n")); flush.console()  
  fit <- inla(as.formula(form), data=data, family="binomial", Ntrials=N, 
              control.inla=list(strategy="gaussian"),
              control.fixed=list(prec.intercept=0.001, mean=prior_mean, prec=prior_precision),
              control.mode=list(result=initial, restart=TRUE), 
              control.compute=list(config=TRUE), 
              num.threads=6, verbose=T)        
} 

# save the model fit
save(fit, file=paste("/clustertmp/counties/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu,  
                     "/fit_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", division, "_svy_sample_", svy_sample, "_sex_", sex, ".rdata", sep=""))
