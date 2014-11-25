

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
## Description: Compile and plot the validation results 
########################################################################################################################

library(data.table)
library(reshape)
library(lattice)
library(latticeExtra)

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

load(file=paste("05_small_area_models/inter/", outcome, "/models.rdata", sep=""))  
  
## Load all validation results 
validation <- lapply(split(models, rownames(models)), function(mm) { 
  load(paste("05_small_area_models/output/", outcome, "/validation_results_model_", mm$family, "_race_", mm$race, "_marriage_", mm$marriage, "_edu_", mm$edu, ".rdata", sep=""))
  data.table(cbind(validation, model=mm$family, race=mm$race, marriage=mm$marriage, edu=mm$edu))
})
validation <- rbindlist(validation)
  
## Calculate means and quantiles
validation <- validation[, c(list(type = c("lwr", "med", "upr")), lapply(.SD, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))), by='model,race,marriage,edu,sex,samp', .SDcols=c("rmse", "mre", "cor", "coverage")]
validation <- melt(validation, id.vars=c("model", "race", "marriage", "edu", "sex", "samp", "type"))
validation <- reshape(validation, direction="wide", idvar=c("model", "race", "marriage", "edu", "sex", "samp", "variable"), timevar="type")
names(validation) <- gsub("value.", "", names(validation))

## Make plots of the results  
validation$model <- factor(validation$model, levels=1:4, labels=c("Naive", "Geospatial", "Covariate", "Full"))
validation$variant <- factor(paste("R=", validation$race, ", M=", validation$marriage, ", E=", validation$edu, sep=""))
validation$samp <- factor(validation$samp, levels=sort(sampling), labels=c(sort(sampling[sampling!=9999]), "In Sample"))
validation$sex <- factor(validation$sex, levels=1:2, labels=c("Male", "Female"))

pdf(paste("05_small_area_models/output/", outcome, "/validation_results.pdf", sep=""), width=14, height=8)
for (M in sort(unique(validation$variable))) {

# plot 1: variant in panel, sampling size on the x-axis, model family in color 
  colors <- brewer.pal(5, "Set1")
  temp <- validation[validation$variable == M,] 
	p <- xyplot(med ~ samp|variant + sex, groups=model,
              data=temp, lwr=temp$lwr, upr=temp$upr,
              xlab="Sampling Level", ylab=M, main=paste(outcome, ": ", M, sep=""), type=c("g","b"),
              panel = panel.superpose,
              panel.groups=function(x,y,lwr,upr,col.line="black",subscripts,group.number,...) {
                shift <- (group.number - 2.5)/20
                panel.arrows(x+shift, lwr[subscripts], x+shift, upr[subscripts], length=0.05, unit="native", angle=90, code=3, col=col.line)
                panel.xyplot(x+shift,y,col.line=col.line,...)
              },
              scales=list(alternating=F),
              par.strip.text=list(cex=0.7),
              par.settings=list(superpose.symbol=list(col=colors, pch=19),
                                superpose.line=list(col=colors),
                                strip.background=list(col="gray80")),
              auto.key=list(space="right"))
  print(useOuterStrips(p))

# plot 2: sampling size in panel, variant on the x-axis, and model-type in colors
  if (length(unique(temp$variant)) > 1) { 
    colors <- rep(brewer.pal(4, "Set1"), 3)
    symbols <- rep(c(17,15,18), each=4)
    lines <- rep(c(1,2,3), each=4)
  	p <- xyplot(med ~ model|samp + sex, groups=variant,
                data=temp, lwr=temp$lwr, upr=temp$upr,
                xlab="Model Family", ylab=M, main=paste(outcome, ": ", M, sep=""), type=c("g","p"),
                prepanel = function(x,y,upr,lwr,subscripts,...) list(ylim=range(c(lwr[subscripts], upr[subscripts]))),
                panel = panel.superpose,
                panel.groups=function(x,y,lwr,upr,col.line="black",subscripts,group.number,...) {
                  shift <- (group.number - 4.5)/10
                  if(sum(lwr[subscripts]) < sum(upr[subscripts])) panel.arrows(x+shift, lwr[subscripts], x+shift, upr[subscripts], length=0.05, unit="native", angle=90, code=3, col=col.line)
                  panel.xyplot(x+shift,y,col.line=col.line,...)
                },
                scales=list(alternating=F, y=list(relation="free"), x=list(rot=30)),
                par.settings=list(superpose.symbol=list(col=colors, pch=symbols),
                                  superpose.line=list(col=colors, lty=lines),
                                  strip.background=list(col="gray80")),
                auto.key=list(space="right"))
    print(useOuterStrips(combineLimits(p, margin.y=2, adjust.labels=F)))
  }
}      

dev.off()
