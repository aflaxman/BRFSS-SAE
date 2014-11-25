

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
## Description: Generate county-year level simulations from each division/sex-specific model
########################################################################################################################

library(INLA)
library(data.table)
library(boot)

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
sex <- as.numeric(commandArgs()[12])
division <- as.numeric(commandArgs()[13])

## Load model
load(file=paste("/clustertmp/counties/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu,
                "/fit_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", division, "_svy_sample_", svy_sample, "_sex_", sex, ".rdata", sep=""))

## Generate simulations (break up into chunks of 100 to reduce memory requirements 
pred_ids <- which(fit$.args$data$type == "pred")
sims <- lapply(1:10, function(x) {
  cat(paste(Sys.time(), x, "\n")); flush.console()
  sapply(inla.posterior.sample(100, fit), function(y) y$latent[pred_ids,])
})
sims <- do.call("cbind", sims)
colnames(sims) <- paste("s", 1:ncol(sims), sep="")
sims <- data.table(cbind(fit$.args$data[type == "pred", list(fips, year, age, offset)], inv.logit(sims)), key=c("fips", "year", "age", "offset"))

## Save simulations
save(sims, file=paste("/clustertmp/counties/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu,
                      "/sims_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", division, "_svy_sample_", svy_sample, "_sex_", sex, ".rdata", sep=""))
