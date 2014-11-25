

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
## Description: Make plots of model parameters (first imputation only)
########################################################################################################################

library(lattice)
library(latticeExtra)
library(plyr)

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
  
## Get county and divisions codes
divisions <- read.csv("01_counties/census_regions_divisions.csv", stringsAsFactors=F)[, c("division2", "division2_fips", "state_fips")] 
divisions$abr <- sapply(divisions$division2, function(x) paste(strsplit(x, "")[[1]][gregexpr("[[:upper:]]", x)[[1]]], collapse=""))

load("01_counties/output/merged_counties.rdata")
codes <- unique(codes[, c("state_fips", "m.fips")])
names(codes)[2] <- "fips"
codes <- merge(codes, divisions[, c("state_fips", "abr")], all.x=T)
names(codes)[3] <- "division"
codes <- codes[order(codes$division, codes$fips),]
codes$cnty <- unlist(tapply(codes$fips, codes$division, function(x) 1:length(x)))

### Load full models and compile model components across sub-models ----------------------------------------------------

load(paste("05_small_area_models/inter/", outcome, "/model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, "/fit_imp_", min(imputations, 1), "_samp_9999_iter_0.rdata", sep=""))

## Fixed effects  
fixed <- do.call("rbind", lapply(names(fit), function(m) data.frame(sm = m, param = rownames(fit[[m]]$summary.fixed), mean = fit[[m]]$summary.fixed[, "mean"], lb = fit[[m]]$summary.fixed[, "0.025quant"], ub = fit[[m]]$summary.fixed[, "0.975quant"])))
fixed$param <- factor(fixed$param, levels=c("(Intercept)", paste("factor(age)", 1:(length(age_groups)-2), sep=""),
                                            unlist(lapply(c("race_black", "race_hisp", "race_native", "race_other"), function(x) grep(x, unique(fixed$param), value=T))),
                                            unlist(lapply(c("marriage_current", "marriage_never"), function(x) grep(x, unique(fixed$param), value=T))),
                                            unlist(lapply(c("edu_HS_grad", "edu_some_college", "edu_college_grad"), function(x) grep(x, unique(fixed$param), value=T))),
                                            grep("phone", unique(fixed$param), value=T),
                                            area_vars, "time_ln", "offset"))
fixed$division <- factor(as.numeric(substr(fixed$sm, 4, 4)), levels=unique(divisions$division2_fips), labels=unique(divisions$abr))
fixed$sex <- factor(substr(fixed$sm, 5, 5), levels=c("1", "2"), labels=c("M", "F"))
fixed$svy_sample <- factor(substr(fixed$sm, 6, 6), levels=c("l", "b"), labels=c("landline", "combined"))
fixed$sex_sample <- factor(paste(fixed$sex, " (", fixed$svy_sample, ")", sep=""))
rownames(fixed) <- 1:nrow(fixed)

## Time random effects
time <- lapply(names(fit), function(m) data.frame(sm = m, year = fit[[m]]$summary.random$time$ID, mean = fit[[m]]$summary.random$time[, "mean"], lb = fit[[m]]$summary.random$time[, "0.025quant"], ub = fit[[m]]$summary.random$time[, "0.975quant"])) 
time <- do.call("rbind", time)
time$division <- factor(as.numeric(substr(time$sm, 4, 4)), levels=unique(divisions$division2_fips), labels=unique(divisions$abr))
time$sex <- factor(substr(time$sm, 5, 5), levels=c("1", "2"), labels=c("M", "F"))
time$svy_sample <- factor(substr(time$sm, 6, 6), levels=c("l", "b"), labels=c("landline", "combined"))
time$year <- time$year - min(time$year)
if (correct_cells) time$year <- time$year + ifelse(time$svy_sample == "landline", year1, 2011) else time$year <- time$year + year1
rownames(time) <- 1:nrow(time)
  
## County random effects 
cnty_iid <- lapply(names(fit), function(m) data.frame(sm = m, cnty_iid = fit[[m]]$summary.random$cnty_iid$ID, mean = fit[[m]]$summary.random$cnty_iid[, "mean"], lb = fit[[m]]$summary.random$cnty_iid[, "0.025quant"], ub = fit[[m]]$summary.random$cnty_iid[, "0.975quant"])) 
cnty_iid <- do.call("rbind", cnty_iid)
cnty_iid$division <- factor(as.numeric(substr(cnty_iid$sm, 4, 4)), levels=unique(divisions$division2_fips), labels=unique(divisions$abr))
cnty_iid$sex <- factor(substr(cnty_iid$sm, 5, 5), levels=c("1", "2"), labels=c("M", "F"))
cnty_iid$svy_sample <- factor(substr(cnty_iid$sm, 6, 6), levels=c("l", "b"), labels=c("landline", "combined"))
rownames(cnty_iid) <- 1:nrow(cnty_iid)

if ("cnty_icar" %in% names(fit[[1]]$summary.random)) { 
  cnty_icar <- lapply(names(fit), function(m) data.frame(sm = m, cnty_icar = fit[[m]]$summary.random$cnty_icar$ID, mean = fit[[m]]$summary.random$cnty_icar[, "mean"], lb = fit[[m]]$summary.random$cnty_icar[, "0.025quant"], ub = fit[[m]]$summary.random$cnty_icar[, "0.975quant"])) 
  cnty_icar <- do.call("rbind", cnty_icar)
  cnty_icar$division <- factor(as.numeric(substr(cnty_icar$sm, 4, 4)), levels=unique(divisions$division2_fips), labels=unique(divisions$abr))
  cnty_icar$sex <- factor(substr(cnty_icar$sm, 5, 5), levels=c("1", "2"), labels=c("M", "F"))
  cnty_icar$svy_sample <- factor(substr(cnty_icar$sm, 6, 6), levels=c("l", "b"), labels=c("landline", "combined"))
  rownames(cnty_icar) <- 1:nrow(cnty_icar)
}

## Interactions 
int <- lapply(names(fit), function(m) data.frame(sm = m, int = fit[[m]]$summary.random$int$ID, mean = fit[[m]]$summary.random$int[, "mean"], lb = fit[[m]]$summary.random$int[, "0.025quant"], ub = fit[[m]]$summary.random$int[, "0.975quant"])) 
int <- do.call("rbind", int)
int <- ddply(int, "sm", function(x) {
  if (correct_cells & substr(x$sm[1], 6, 6) == "b") years <- 2011:year2 else years <- year1:year2
  n <- nrow(x)/length(years)
  cbind(x, year = rep(years, n), cnty = rep(1:n, each=length(years)))
})
int$division <- factor(as.numeric(substr(int$sm, 4, 4)), levels=unique(divisions$division2_fips), labels=unique(divisions$abr))
int$sex <- factor(substr(int$sm, 5, 5), levels=c("1", "2"), labels=c("M", "F"))
int$svy_sample <- factor(substr(int$sm, 6, 6), levels=c("l", "b"), labels=c("landline", "combined"))
rownames(int) <- 1:nrow(int)  

### Make plots ---------------------------------------------------------------------------------------------------------

pdf(paste("05_small_area_models/output/", outcome, "/model_parameters_model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, ".pdf", sep=""), width=14, height=8)
if (correct_cells == 1) {
  colors1 <- brewer.pal(6, "Paired")[c(6,5,2,1)] 
  colors2 <- brewer.pal(6, "Paired")[c(4,3)]
} else { 
  colors1 <- brewer.pal(6, "Paired")[c(6,2)]
  colors2 <- brewer.pal(6, "Paired")[4]
}

## Plot the fixed effects (except time) 
plotfe <- function(fixed, main, ci=TRUE) { 
  p <- xyplot(mean ~ division|param, data=fixed, groups=sex_sample, lwr=fixed$lb, upr=fixed$ub,
    scales=list(alternating=F, cex=0.8, y=list(rot=0), x=list(rot=45)), as.table=T, 
    xlab="", ylab="", main=paste(main, outcome, sep=" - "), 
    par.settings=list(superpose.symbol=list(col=colors1, pch=15), strip.background=list(col="gray")), 
    par.strip.text=list(cex=0.9),
    prepanel=function(x, y, lwr, upr, subscripts, ...) list(ylim = if (ci) range(upr[subscripts], lwr[subscripts]) else range(y), finite=T),
    panel=panel.superpose,
    panel.groups=function(x, y, upr, lwr, subscripts, col.symbol="black", group.number, ...) {
      panel.abline(h=0, lwd=2, col="gray")
      panel.abline(v=1:7+0.5, col="gray80", lty=3)
      plot_x <- x + (group.number-2.5)/8
      if (ci) panel.segments(plot_x, lwr[subscripts], plot_x, upr[subscripts], col=col.symbol)
      panel.xyplot(plot_x, y, col=col.symbol, pch=19, cex=1)
    },
    auto.key=list(x=0.99, y=0.01, corner=c(1,0), columns=2, background="white")) 
  print(p)
}  

try(plotfe(fixed[grepl("Intercept|offset|factor.age.[[:digit:]]$", fixed$param),], main="Fixed effects: Intercept and age dummies", F))  
try(plotfe(fixed[grepl("Intercept|offset|factor.age.[[:digit:]]$", fixed$param),], main="Fixed effects: Intercept and age dummies"))
try(plotfe(fixed[grepl("race", fixed$param),], main="Fixed effects: Race", F))
try(plotfe(fixed[grepl("race", fixed$param),], main="Fixed effects: Race"))
try(plotfe(fixed[grepl("marriage", fixed$param),], main="Fixed effects: Marriage", F))
try(plotfe(fixed[grepl("marriage", fixed$param),], main="Fixed effects: Marriage"))
try(plotfe(fixed[grepl("edu", fixed$param),], main="Fixed effects: Education", F))
try(plotfe(fixed[grepl("edu", fixed$param),], main="Fixed effects: Education"))
try(plotfe(fixed[grepl("phone", fixed$param),], main="Fixed effects: Phone", F))
try(plotfe(fixed[grepl("phone", fixed$param),], main="Fixed effects: Phone"))
try(plotfe(fixed[!grepl("Intercept|offset|factor.age.|race|marriage|edu|phone|time", fixed$param),], main="Fixed effects: Area-level covariates", F))
try(plotfe(fixed[!grepl("Intercept|offset|factor.age.|race|marriage|edu|phone|time", fixed$param),], main="Fixed effects: Area-level covariates"))

## Plot the time fixed + random effects
time <- merge(time, fixed[fixed$param == "time_ln", c("sm", "mean")], all.x=T, by="sm")
time$mean <- time$mean.x + time$mean.y*as.numeric(as.factor(time$year))
time <- time[order(time$sm, time$year),]

p <- xyplot(mean ~ year|division, groups=paste(sex, " (", svy_sample, ")", sep=""), data=time, 
  scales=list(alternating=F, cex=0.8, y=list(rot=0), x=list(rot=45)), as.table=T,
  xlab="Year", ylab="", main=paste("Random & fixed effects: Time", outcome, sep=" - "), 
  par.settings=list(superpose.line=list(col=colors1, lwd=2), strip.background=list(col="gray")),
  par.strip.text=list(cex=0.9), 
  panel = function(type="l", ...) { 
    panel.xyplot(type="g", ...)
    panel.abline(h=0, lwd=3, col="gray40")
    panel.xyplot(type="l", ...)
  },
  auto.key=list(points=F, lines=T, x=0.99, y=0.01, corner=c(1,0), columns=2, background="white"))
print(p) 

## Plot the space random effects 
p <- xyplot(mean ~ cnty_iid|division + sex, groups=svy_sample, data=cnty_iid, 
  scales=list(alternating=F, cex=0.8, axs="i", y=list(rot=0), x=list(relation="free", rot=45)),
  as.table=T, xlab="County", ylab="", main=paste("Random effects: Space (IID)", outcome, sep=" - "), 
  par.settings=list(superpose.symbol=list(col=colors2, pch=19, cex=0.25), strip.background=list(col="gray")),
  par.strip.text=list(cex=0.9),
  panel = function(type="l", ...) { 
    panel.xyplot(type="g", ...)
    panel.abline(h=0, lwd=3, col="gray40")
    panel.xyplot(type="p", ...)
  },
  auto.key= if(correct_cells) list(points=T, x=0.99, y=0.01, corner=c(1,0), background="white") else FALSE)
print(useOuterStrips(combineLimits(resizePanels(p, w=tapply(cnty_iid$cnty_iid, cnty_iid$division, max)), extend=F)))

if (!is.null(cnty_icar)) { 
  p <- xyplot(mean ~ cnty_icar|division + sex, groups=svy_sample, data=cnty_icar, 
    scales=list(alternating=F, cex=0.8, axs="i", y=list(rot=0), x=list(relation="free", rot=45)),
    as.table=T, xlab="County", ylab="", main=paste("Random effects: Space (ICAR)", outcome, sep=" - "), 
    par.settings=list(superpose.symbol=list(col=colors2, pch=19, cex=0.25), strip.background=list(col="gray")),
    par.strip.text=list(cex=0.9),
    panel = function(type="l", ...) { 
      panel.xyplot(type="g", ...)
      panel.abline(h=0, lwd=3, col="gray40")
      panel.xyplot(type="p", ...)
    },
    auto.key= if(correct_cells) list(points=T, x=0.99, y=0.01, corner=c(1,0), background="white") else FALSE)
  print(useOuterStrips(combineLimits(resizePanels(p, w=tapply(cnty_icar$cnty_icar, cnty_icar$division, max)), extend=F)))
}
    
## Plot the interaction random effects 
for (sex in c("M", "F")) { 
  p <- xyplot(mean ~ cnty|division + factor(year), groups=svy_sample, data=int[int$sex == sex,],
    scales=list(alternating=F, cex=0.8, axs="i", y=list(rot=0), x=list(relation="free", rot=45)),
    as.table=T, xlab="County", ylab="", main=paste("Random effects: Interaction (sex=", sex, ") - ", outcome, sep=""), 
    par.settings=list(superpose.symbol=list(col=colors2, pch=19, cex=0.25), strip.background=list(col="gray")),
    par.strip.text=list(cex=0.9),
    panel = function(type="l", ...) { 
      panel.xyplot(type="g", ...)
      panel.abline(h=0, lwd=3, col="gray40")
      panel.xyplot(type="p", ...)
    },
    auto.key= if(correct_cells) list(points=T, x=0.99, y=0.01, corner=c(1,0), background="white") else FALSE)  
  print(useOuterStrips(combineLimits(resizePanels(p, w=tapply(int$cnty, int$division, max)), extend=F)))
}

dev.off()
