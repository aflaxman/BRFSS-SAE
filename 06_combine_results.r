

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
## Description: Compile predictions across imputations (if necessary) and sex. 
##              Compile validation results across all iterations and sampling levels
########################################################################################################################

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

try(detach(settings), silent=T)
get_settings <- read.csv("05_small_area_models/model_settings.csv", skip=1, stringsAsFactors=F)[,-1]
settings <- lapply(1:nrow(get_settings), function(x) strsplit(get_settings[x, outcome], split=";")[[1]])
settings <- lapply(settings, function(x) if(sum(is.na(as.numeric(x)))==0) as.numeric(x) else x)
names(settings) <- get_settings[,1]
attach(settings)
rm(get_settings, settings)
  
## Load estimates from full models
all_est <- all_diff <- all_pct_diff <- all_arc <- all_rank_est <- all_rank_diff <- all_rank_pct_diff <- all_rank_arc <- NULL 
for (imp in if (imputations == 0) 0 else 1:imputations) { 
  load(paste("05_small_area_models/inter/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, "/preds_imp_", imp, "_samp_9999_iter_0.rdata", sep="")) 
  all_est[[length(all_est) + 1]] <- cbind(imp = imp, est) 
  all_diff[[length(all_diff) + 1]] <- cbind(imp = imp, diff)
  all_pct_diff[[length(all_pct_diff) + 1]] <- cbind(imp = imp, pct_diff)
  all_arc[[length(all_arc) + 1]] <- cbind(imp = imp, arc) 
  all_rank_est[[length(all_rank_est) + 1]] <- cbind(imp = imp, rank_est)  
  all_rank_diff[[length(all_rank_diff) + 1]] <- cbind(imp = imp, rank_diff)
  all_rank_pct_diff[[length(all_rank_pct_diff) + 1]] <- cbind(imp = imp, rank_pct_diff)
  all_rank_arc[[length(all_rank_arc) + 1]] <- cbind(imp = imp, rank_arc)
  rm(est, uncorrected_est, diff, pct_diff, arc, rank_est, rank_diff, rank_pct_diff, rank_arc); gc()  
} 
all <- list(all_est, all_diff, all_pct_diff, all_arc, all_rank_est, all_rank_diff, all_rank_pct_diff, all_rank_arc) 
all <- lapply(all, function(x) data.frame(rbindlist(x)))
rm(all_est, all_diff, all_pct_diff, all_arc, all_rank_est, all_rank_diff, all_rank_pct_diff, all_rank_arc); gc()
  
## Combine across imputations, if necessary, and format
if (imputations == 0) { 
  format_final <- function(data, idvars, estvar, ot_sig=F) {
    if (ot_sig) { 
      data <- data[do.call("order", data[idvars]), c(idvars, paste(estvar, c("", ".var", ".lb", ".ub"), sep=""), "ot_sig_change")]
      names(data) <- c(idvars, c("direct", "var", "lb", "ub", "ot_sig_change"))
    } else { 
      data <- data[do.call("order", data[idvars]), c(idvars, paste(estvar, c("", ".var", ".lb", ".ub"), sep=""))]
      names(data) <- c(idvars, c("direct", "var", "lb", "ub"))      
    } 
    rownames(data) <- 1:nrow(data)
    data
  }
  
  est <- format_final(all[[1]], c("fips", "sex", "year", "svy_sample"), "prev", F)
  diff <- format_final(all[[2]], c("fips", "sex", "years", "svy_sample"), "diff", T)
  pct_diff <- format_final(all[[3]], c("fips", "sex", "years", "svy_sample"), "pct_diff", T)
  arc <- format_final(all[[4]], c("fips", "sex", "years", "svy_sample"), "arc", T)
  rank_est <- format_final(all[[5]], c("fips", "sex", "year", "svy_sample"), "rank", F)
  rank_diff <- format_final(all[[6]], c("fips", "sex", "years", "svy_sample"), "rank", F)
  rank_pct_diff <- format_final(all[[7]], c("fips", "sex", "years", "svy_sample"), "rank", F)
  rank_arc <- format_final(all[[8]], c("fips", "sex", "years", "svy_sample"), "rank", F)    
  rm(all); gc() 

} else { 
  combine_imputations <- function(data, idvars, estvar, ot_sig=F) { 
    data$ot_sig_change <- NULL
    data <- reshape(data, direction="wide", idvar=idvars, timevar="imp", v.names=paste(estvar, c("", ".var"), sep=""), drop=paste(estvar, c(".lb", ".ub"), sep=""))
    data$mean_est <- apply(data[, paste(estvar, 1:imputations, sep=".")], 1, mean)
    data$mean_var <- apply(data[, paste(estvar, "var", 1:imputations, sep=".")], 1, mean)
    data$between_var <- apply(data[, c("mean_est", paste(estvar, 1:imputations, sep="."))], 1, function(x) sum((x[-1]-x[1])^2)/(imputations-1))
    data$total_var <- (1+1/imputations) * data$between_var + data$mean_var
    data$df <- round((imputations-1) * (data$total_var/((1+1/imputations)*data$between_var))^2)
    data$lwr <- apply(data[,c("mean_est", "df", "total_var")], 1, function(x) x[1] + qt(0.025, df=x[2]) * sqrt(x[3]))
    data$upr <- apply(data[,c("mean_est", "df", "total_var")], 1, function(x) x[1] + qt(0.975, df=x[2]) * sqrt(x[3])) 
    if (ot_sig) { 
      lwr_1t <- apply(data[,c("mean_est", "df", "total_var")], 1, function(x) x[1] + qt(0.05, df=x[2]) * sqrt(x[3]))
      upr_1t <- apply(data[,c("mean_est", "df", "total_var")], 1, function(x) x[1] + qt(0.95, df=x[2]) * sqrt(x[3]))     
      ot_sig_change <- as.numeric(lwr_1t > 0 | upr_1t < 0)
      rm(lwr_1t, upr_1t)
    }
    data <- data[do.call("order", data[idvars]), c(idvars, "mean_est", "total_var", "lwr", "upr")]
    names(data) <- c(idvars, c("direct", "var", "lb", "ub"))
    rownames(data) <- 1:nrow(data)
    if (ot_sig) {
      data$ot_sig_change <- ot_sig_change
      rm(ot_sig_change)
    }
    
    if (estvar == "prev") {
      data$lb <- pmax(0, data$lb)
      data$ub <- pmin(1, data$ub)
    }
    if (estvar == "rank") { 
      data$direct <- round(data$direct)
      data$lb <- pmax(1, floor(data$lb))
      data$ub <- pmin(length(unique(data$fips)), ceiling(data$ub))
    }
    data    
  } 
  
  est <- combine_imputations(all[[1]], c("fips", "sex", "year", "svy_sample"), "prev", F)
  diff <- combine_imputations(all[[2]], c("fips", "sex", "years", "svy_sample"), "diff", T)
  pct_diff <- combine_imputations(all[[3]], c("fips", "sex", "years", "svy_sample"), "pct_diff", T)
  arc <- combine_imputations(all[[4]], c("fips", "sex", "years", "svy_sample"), "arc", T)    
  rank_est <- combine_imputations(all[[5]], c("fips", "sex", "year", "svy_sample"), "rank", F)
  rank_diff <- combine_imputations(all[[6]], c("fips", "sex", "years", "svy_sample"), "rank", F)
  rank_pct_diff <- combine_imputations(all[[7]], c("fips", "sex", "years", "svy_sample"), "rank", F)
  rank_arc <- combine_imputations(all[[8]], c("fips", "sex", "years", "svy_sample"), "rank", F)        
  rm(all); gc() 
  
  # recalculate change point estimates for consistency with prevalence point estimates
  recalc <- lapply(diff_years, function(yy) { 
    y1 <- as.numeric(substr(yy, 1, 4))
    y2 <- as.numeric(substr(yy, 6, 9))
    temp <- reshape(est[est$year %in% c(y1, y2), c("fips", "sex", "year", "svy_sample", "direct")], direction="wide", idvar=c("fips", "sex", "svy_sample"), timevar="year")
    temp$diff <- temp[, paste0("direct.", y2)] - temp[, paste0("direct.", y1)]
    temp$pct_diff <- temp$diff/temp[, paste0("direct.", y1)]
    temp$arc <- log(temp[, paste0("direct.", y2)]/temp[, paste0("direct.", y1)])/(y2-y1)
    temp <- temp[, c("fips", "sex", "svy_sample", "diff", "pct_diff", "arc")]
    temp$years <- yy
    temp
  })
  recalc <- do.call("rbind", recalc)
  
  for (object in c("diff", "pct_diff", "arc")) { 
    orig <- get(object)  
    orig <- merge(orig, recalc)
    orig$direct <- orig[, object]
    orig <- orig[, names(get(object))]
    assign(object, orig)
    rm(orig)
  }
  
  # recalculate rank point estimates for consistency with prevalence/change point estimates 
  for (object in c("est", "diff", "pct_diff", "arc")) { 
    recalc <- get(object)
    recalc <- data.table(recalc[recalc$fips > 100,])
    recalc[, rank := rank(direct), by=c("sex", grep("year", names(recalc), value=T), "svy_sample")]
    recalc <- recalc[, c("fips", "sex", grep("year", names(recalc), value=T), "svy_sample", "rank"), with=F]
    orig <- get(paste0("rank_", object))
    orig <- merge(orig, recalc)
    orig$direct <- orig$rank
    orig <- orig[, names(get(paste0("rank_", object)))]
    orig$lb <- pmin(orig$lb, orig$direct)
    orig$ub <- pmax(orig$ub, orig$direct)
    orig <- orig[do.call("order", orig),]
    assign(paste0("rank_", object), orig)
    rm(orig)
  }
} 
  
save(est, diff, pct_diff, arc, rank_est, rank_diff, rank_pct_diff, rank_arc, file=paste("05_small_area_models/output/", outcome, "/predictions_model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, ".rdata", sep=""))

## Load validation results 
if (!is.na(validation_type)) {
  validation <- NULL
  for (samp in sort(unique(c(sampling, 9999)))) { 
    for (iter in if(samp == 9999) 0 else 1:iterations) { 
      load(paste("05_small_area_models/inter/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, "/preds_imp_0_samp_", samp, "_iter_", iter, ".rdata", sep=""))
      validation <- rbind(validation, data.frame(samp = samp, iter = iter, val_output))
    } 
  } 

  save(validation, file=paste("05_small_area_models/output/", outcome, "/validation_results_model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, ".rdata", sep=""))  
} 
