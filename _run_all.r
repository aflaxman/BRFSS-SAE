

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
## Description: Submit all jobs to run specified models to the cluster
########################################################################################################################

rm(list=ls())
root <- ifelse(Sys.info()[1]=="Windows", "J:/", "/home/j/")
setwd(paste(root, "Project/NHIS/sae", sep=""))

### Get settings -------------------------------------------------------------------------------------------------------

outcome <- commandArgs()[3]

try(detach(settings), silent=T)
get_settings <- read.csv("05_small_area_models/model_settings.csv", skip=1, stringsAsFactors=F)[,-1]
settings <- lapply(1:nrow(get_settings), function(x) strsplit(get_settings[x, outcome], split=";")[[1]])
settings <- lapply(settings, function(x) if(sum(is.na(as.numeric(x)))==0) as.numeric(x) else x)
names(settings) <- get_settings[,1]
attach(settings)
rm(get_settings)

if (FALSE|is.na(sampling)) sampling <- 9999 # turn on to skip validation 
replace <- F # Set to T to submit all jobs; Set to F to skip jobs that are currently running or where the output file exists
 
## Create inter and output and temp folders for this outcome
dir.create(paste("05_small_area_models/inter/", outcome, sep=""), showWarnings=F)
dir.create(paste("05_small_area_models/output/", outcome, sep=""), showWarnings=F)
dir.create(paste("/clustertmp/counties/", outcome, sep=""), showWarnings=F)

save(settings, file=paste("05_small_area_models/output/", outcome, "/settings.rdata", sep=""))
rm(settings) 
 
## Define all models and sub-models 
models <- expand.grid(family = sort(model_families), race = model_race, edu = model_edu, marriage = model_marriage, stringsAsFactors=F)
sub_models <- expand.grid(division = 1:8, sex = 1:2, svy_sample = if (correct_cells) c("l", "b") else "b", stringsAsFactors=F)
save(models, sub_models, file=paste("05_small_area_models/inter/", outcome, "/models.rdata", sep=""))

### Submit jobs --------------------------------------------------------------------------------------------------------

## Define qsub function 
qsub <- function(name, code, arguments=NULL, hold=NULL, slots=1, test=(root == "J:/")) {
  if(!is.null(arguments)) arguments <- paste(paste("\"", arguments, "\"", sep=""), collapse=" ")
  if(!is.null(hold) & length(hold)>1) hold <- paste(hold, collapse=",")
  x <- paste("/usr/local/bin/SGE/bin/lx24-amd64/qsub -cwd",
             if (slots>1) paste("-pe multi_slot", slots),
             "-N", name,
             if (!is.null(hold)) paste("-hold_jid", hold),
             "r_shell.sh" , code,
             if (!is.null(arguments)) arguments,
             sep=" ")
  if (test) print(x) else system(x)
}

count <- 0 

## Find out if any jobs have already been submitted
running <- system("qstat -r | grep jobname | awk {'print $3'}", intern=T)
running <- grep(outcome, running, value=T)
running

## Submit code to generate population weights
jobname_pop <- paste("pop_wts", outcome, sep="_") 
if (replace | (!file.exists(paste("05_small_area_models/inter/", outcome, "/pop_weights.rdata", sep="")) & !jobname_pop %in% running)) { 
  qsub(jobname_pop, "05_small_area_models/code/rw_icar/01_create_pop_weights.r", outcome)
  count <- count + 1
} 
  
## Submit code to prep data 
jobname_data <- paste("data_prep", outcome, sep="_") 
if (replace | (!file.exists(paste("05_small_area_models/inter/", outcome, "/model_data.rdata", sep="")) & !jobname_data %in% running)) {
  qsub(jobname_data, "05_small_area_models/code/rw_icar/02_create_datasets.r", outcome, jobname_pop)
  count <- count + 1
} 

## Loop through models 
all_jobname_comb <- NULL
for (mm in split(models, 1:nrow(models))) { 

## Create model-specific directories 
  inter_dir <- paste("05_small_area_models/inter/", outcome, "/model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, sep="")
  dir.create(inter_dir, showWarnings=F)     
  clustertmp_dir <- paste("/clustertmp/counties/", outcome, "/model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, sep="")
  dir.create(clustertmp_dir, showWarnings=F)
  
## Loop through sampling levels, iterations, and imputations 
  all_jobname_comp <- NULL 
  for (samp in sort(unique(c(sampling, 9999)))) { 
    for (iter in if(samp == 9999) 0 else 1:iterations) { 
      for (imp in if(imputations == 0) 0 else 1:imputations) {   
      
## Loop through submodels 
        all_jobname_sim <- NULL
        for (sm in split(sub_models, 1:nrow(sub_models)))  { 
          if (correct_cells == 1 & sm$svy_sample == "b" & samp != 9999) next
          
## Submit code to fit models and to generate simulations
          jobname_fit <- paste("fit", paste(mm, collapse=""), paste(sm, collapse=""), samp, iter, imp, outcome, sep="_")
          if (replace | (!file.exists(paste(clustertmp_dir, "/fit_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", sm$division, "_svy_sample_", sm$svy_sample, "_sex_", sm$sex, ".rdata", sep="")) & !jobname_fit %in% running)) { 
            qsub(jobname_fit, "05_small_area_models/code/rw_icar/03_fit_model.r", arguments=c(outcome, mm$family, mm$race, mm$marriage, mm$edu, imp, samp, iter, sm$svy_sample, sm$sex, sm$division), hold=jobname_data, slots=10)
            count <- count + 1  
          } 
             
          jobname_sim <- paste("sim", paste(mm, collapse=""), paste(sm, collapse=""), samp, iter, imp, outcome, sep="_")
          if (replace | (!file.exists(paste(clustertmp_dir, "/sims_imp_", imp, "_samp_", samp, "_iter_", iter, "_division_", sm$division, "_svy_sample_", sm$svy_sample, "_sex_", sm$sex, ".rdata", sep="")) & !jobname_sim %in% running)) {  
            qsub(jobname_sim, "05_small_area_models/code/rw_icar/04_gen_sims.r", arguments=c(outcome, mm$family, mm$race, mm$marriage, mm$edu, imp, samp, iter, sm$svy_sample, sm$sex, sm$division), hold=jobname_fit, slots=ifelse(samp==9999, 10, 5))
            count <- count + 1 
          } 
          all_jobname_sim <- c(all_jobname_sim, jobname_sim)
        } # close sub-model loop
        
## Submit code to compile all submodels 
        jobname_comp <- paste("comp", paste(mm, collapse=""), samp, iter, imp, outcome, sep="_")
        if (replace | (!file.exists(paste("05_small_area_models/inter/", outcome, "/model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, "/preds_imp_", imp, "_samp_", samp, "_iter_", iter, ".rdata", sep="")) & !jobname_comp %in% running)) { 
          qsub(jobname_comp, "05_small_area_models/code/rw_icar/05_compile_sims.r", arguments=c(outcome, mm$family, mm$race, mm$marriage, mm$edu, imp, samp, iter), hold=all_jobname_sim, slots=ifelse(samp==9999, 10, 2))
          count <- count + 1
        } 
        all_jobname_comp <- c(all_jobname_comp, jobname_comp)
        
      } # close imputations loop
    } # close iterations loop 
  } # close sampling levels loop 
  
## Submit code to combine all sampling levels, iterations, and imputations for this model  
  jobname_comb <- paste("comb", paste(mm, collapse=""), outcome, sep="_")
  if (replace | (!file.exists(paste("05_small_area_models/output/", outcome, "/predictions_model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, ".rdata", sep="")) & !jobname_comb %in% running)) { 
    qsub(jobname_comb, "05_small_area_models/code/rw_icar/06_combine_results.r", arguments=c(outcome, mm$family, mm$race, mm$marriage, mm$edu), hold=all_jobname_comp)
    count <- count + 1
  } 
  all_jobname_comb <- c(all_jobname_comb, jobname_comb)

## Submit code to plot predictions and model parameters for this model 
  jobname_plotpred <- paste("map", paste(mm, collapse=""), outcome, sep="_")
  if (replace | (!file.exists(paste("05_small_area_models/output/", outcome, "/graphs_and_maps_model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, ".pdf", sep="")) & !jobname_plotpred %in% running)) { 
    qsub(jobname_plotpred, "05_small_area_models/code/rw_icar/07_plot_preds.r", arguments=c(outcome, mm$family, mm$race, mm$marriage, mm$edu), hold=jobname_comb)
    count <- count + 1
  } 
  
  jobname_plotmod <- paste("plot", paste(mm, collapse=""), outcome, sep="_")
  if (replace | (!file.exists(paste("05_small_area_models/output/", outcome, "/model_parameters_model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, ".pdf", sep="")) & !jobname_plotmod %in% running)) { 
    qsub(jobname_plotmod, "05_small_area_models/code/rw_icar/08_plot_models.r", arguments=c(outcome, mm$family, mm$race, mm$marriage, mm$edu), hold=jobname_comb)
    count <- count + 1
  }       
  
} # close model loop 

## Submit code to compile and plot validation results across models 
jobname_val <- paste("val", outcome, sep="_")
if (replace | (!file.exists(paste("05_small_area_models/output/", outcome, "/validation_results.pdf", sep="")) & !jobname_val %in% running)) { 
  qsub(jobname_val, "05_small_area_models/code/rw_icar/09_plot_validation.r", arguments=outcome, hold=all_jobname_comb)
  count <- count + 1
} 

## Count jobs submitted
cat(paste("\n\n******", count, "jobs submitted *******")); flush.console()
