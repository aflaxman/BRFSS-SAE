

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
## Description: Make maps and plots of model estimates
########################################################################################################################

library(maptools)
library(sp)
library(lattice)
library(latticeExtra) 

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

mapcolors <- c("#001C73", "#0070FF", "#00A884", "#A3FF73", "#FFFF00", "#FFBF00", "darkorange", "#FF0000", "#730000")     

## Load predictions
load(file=paste("05_small_area_models/output/", outcome, "/predictions_model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, ".rdata", sep=""))
rm(diff, pct_diff, arc, rank_est, rank_diff, rank_pct_diff, rank_arc)

## Load shape files
load("01_counties/output/county_shape_file.rdata")
counties@data <- counties@data[, c("m.fips", "id")]
names(counties@data)[1] <- "fips" 
ids <- counties@data

load("01_counties/output/state_shape_file.rdata")
states <- list("sp.lines", as(states, "SpatialLines"), lwd=0.4)

## Start pdf
pdf(paste("05_small_area_models/output/", outcome, "/graphs_and_maps_model_", model, "_race_", race, "_marriage_", marriage, "_edu_", edu, ".pdf", sep=""), width=14, height=8)
  
## Make levels maps
years <- c(year1, round(median(year1:year2)), year2)
data <- est[est$fips > 100 & est$year %in% years & est$sex < 3, c("fips", "sex", "year", "svy_sample", "direct")] 
data$direct <- 100*data$direct
ylim <- c(floor(min(data$direct)), ceiling(max(data$direct)))
data <- reshape(data, direction="wide", idvar=c("fips", "year", "svy_sample"), timevar="sex", v.names="direct")
data <- reshape(data, direction="wide", idvar=c("fips", "svy_sample"), timevar="year", v.names=c("direct.1", "direct.2"))
data <- merge(data, ids, by="fips", all=T)
data <- data[order(data$id, data$svy_sample),]

for (ss in sort(unique(data$svy_sample))) {     
  map <- counties
  map@data <- data[data$svy_sample == ss,]
  print(spplot(map, zcol=unique(c(paste("direct.1.", years, sep=""), paste("direct.2.", years, sep=""))),
               as.table=T, layout=c(length(unique(years)), 2), 
               at=seq(ylim[1], ylim[2], length.out=101), 
               col=NA, main=paste("Estimates by sex and year (", ss, "), model ", model, sep=""), 
               par.settings=list(regions=list(col=colorRampPalette(mapcolors)(100))),
               strip=strip.custom(bg="gray", fg="gray", factor.levels=unique(c(paste("Male", years), paste("Female", years)))),
               sp.layout=states)) 
  rm(map); gc()   
} 

## Plot estimates and confidence intervals  
range <- range(c(est$lb[est$fips > 100], est$ub[est$fips > 100]), na.rm=T)
for (year in sort(unique(est$year))) { 
  data <- est[est$fips > 100 & est$year == year & est$sex < 3, ]
  data$sex <- factor(data$sex, levels=1:2, labels=c("Males", "Females")) 
  data$svy_sample <- factor(factor(data$svy_sample, levels=c("b", "l"), c("combined sample", "landline only sample"))) 
  p <- xyplot(direct ~ factor(fips)|svy_sample+sex, data=data, 
              xlab="Fips (sorted by level)", ylab=outcome, main=paste("Estimates in ", year, " with CI, model", model), 
              scales=list(alternating=F, x=list(draw=F)), as.table=T, ylim=range, 
              par.settings=list(strip.background=list(col="gray")), 
              panel=function(x, y, subscripts, all=data) { 
                temp <- all[subscripts,]
                temp <- temp[order(temp$direct),]
                panel.grid()
                panel.segments(1:nrow(temp), temp$lb, 1:nrow(temp), temp$ub, col="lightpink")
                panel.xyplot(1:nrow(temp), temp$direct, subscripts=subscripts, col="red", pch=19)
              })
  print(useOuterStrips(p))
} 

## Plot cell phone correction
if (correct_cells) { 
  data <- reshape(est[est$fips > 100 & est$sex < 3, c("fips", "sex", "year", "svy_sample", "direct")], direction="wide", idvar=c("fips", "sex", "year"), timevar="svy_sample")
  bwplot(100*(direct.b - direct.l) ~ factor(year)|factor(sex), data=data, horizontal=F,
         scales=list(x=list(rot=90)), par.settings=list(strip.background=list(col="gray")), 
         xlab="", ylab="Corrected - Uncorrected Estimates (% points)", main="Effect of cell phone correction", 
         panel=function(...) {
           panel.abline(h=0, lwd=3, col="gray")
           panel.abline(h=c(-5,5), lwd=1, col="gray")
           panel.bwplot(...)
         })    
}   
dev.off()
