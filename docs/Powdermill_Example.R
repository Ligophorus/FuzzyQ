# Powdermill example
# Data from Merritt (2013) Long Term Mammal Data from Powdermill Biological
# Station 1979-1999. Environmental Data Initiative.
# https://doi.org/10.6073/pasta/83c888854e239a79597999895bb61cfe 
# 
# Data License Data and documentation is copyrighted by The Virginia Coast
# Reserve LTER project of the University of Virginia (UVA), and ownership
# remains with the UVA. The UVA grants a license to use the data and
# documentation for academic, and research purposes only, without a fee.
#
# This script uses derivative work which is a modified version and not the
# original data and documentation distributed by the UVA.
# The original data was obtained with the support of NSF Grants BSR-8702333-06,
# DEB-9211772, DEB-9411974, DEB-0080381 and DEB-0621014
#
# Package ID: knb-lter-vcr.67.17 Cataloging System:https://pasta.lternet.edu.
# Data set title: Long Term Mammal Data from Powdermill Biological Station 1979-1999.
# Data set creator:  Joseph Merritt -  
# Metadata Provider:    - Virginia Coast Reserve Long-Term Ecological Research Project 
# Contact:    - Information Manager LTER Network Office  - tech-support@lternet.edu
# Contact:  Joseph Merritt -    - jmerritt@illinois.edu
# Contact:    - Information manager - Virginia Coast Reserve Long-Term Ecological Research Project   - jporter@lternet.edu
# Metadata Link: https://portal.lternet.edu/nis/metadataviewer?packageid=knb-lter-vcr.67.17
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

# Download dataset
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-vcr/67/17/d97866b0e6804e6357cead3543554290" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

dt1 <-read.csv(infile1,header=F 
               ,skip=21
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "ID",     
                 "NEW",     
                 "SPECIES",     
                 "PERIOD",     
                 "TIME",     
                 "DATE",     
                 "QUADR",     
                 "SEX",     
                 "NO",     
                 "WEIGHT",     
                 "OT",     
                 "SC",     
                 "IN",     
                 "PR",     
                 "OP",     
                 "LG"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and
# nominal columns read as numeric or dates read as strings

if (class(dt1$ID)!="factor") dt1$ID<- as.factor(dt1$ID)
if (class(dt1$NEW)!="factor") dt1$NEW<- as.factor(dt1$NEW)
if (class(dt1$SPECIES)!="factor") dt1$SPECIES<- as.factor(dt1$SPECIES)
if (class(dt1$PERIOD)=="factor") dt1$PERIOD <-as.numeric(levels(dt1$PERIOD))[as.integer(dt1$PERIOD) ]               
if (class(dt1$PERIOD)=="character") dt1$PERIOD <-as.numeric(dt1$PERIOD)
if (class(dt1$TIME)!="factor") dt1$TIME<- as.factor(dt1$TIME)
if (class(dt1$DATE)=="factor") dt1$DATE <-as.numeric(levels(dt1$DATE))[as.integer(dt1$DATE) ]               
if (class(dt1$DATE)=="character") dt1$DATE <-as.numeric(dt1$DATE)
if (class(dt1$QUADR)!="factor") dt1$QUADR<- as.factor(dt1$QUADR)
if (class(dt1$SEX)!="factor") dt1$SEX<- as.factor(dt1$SEX)
if (class(dt1$NO)!="factor") dt1$NO<- as.factor(dt1$NO)
if (class(dt1$WEIGHT)=="factor") dt1$WEIGHT <-as.numeric(levels(dt1$WEIGHT))[as.integer(dt1$WEIGHT) ]               
if (class(dt1$WEIGHT)=="character") dt1$WEIGHT <-as.numeric(dt1$WEIGHT)
if (class(dt1$OT)!="factor") dt1$OT<- as.factor(dt1$OT)
if (class(dt1$SC)!="factor") dt1$SC<- as.factor(dt1$SC)
if (class(dt1$IN)!="factor") dt1$IN<- as.factor(dt1$IN)
if (class(dt1$PR)!="factor") dt1$PR<- as.factor(dt1$PR)
if (class(dt1$OP)!="factor") dt1$OP<- as.factor(dt1$OP)
if (class(dt1$LG)!="factor") dt1$LG<- as.factor(dt1$LG)

# Manipulate dataframe to end up with a site x species matrix
# together with year and quadrant
# Start by setting date as factor
dt1$DATE <- as.Date(as.character(dt1$DATE), format='%Y%m%d')
# The use date to create year as factor 
dt1$YR <- as.factor(format(as.Date(dt1$DATE, format="%d/%m/%Y"),"%Y"))
# delete plot names recorded erroneously
keep.quadr <- which(dt1$QUADR %in% as.vector(sapply(LETTERS[1:10],
                                          function(x) paste(x, 1:10, sep=""))))
dt1 <- dt1[keep.quadr, ]
# Add count column in order to aggregate occurrences with reshape
dt1$count <- rep(1, nrow(dt1))
# Reshape as site x species matrix with package reshape
library(reshape)
pdm.dat <-  cast(dt1, YR+QUADR ~ SPECIES, value="count", fun.aggregate=sum)
# Delete species data not included in metadata (see Metadata link above)
keep.spp <- which(colnames(pdm.dat) %in% c("BB", "CG", "DV", "GV", "MF", "MM",
                                           "NI", "NM", "PL", "PM", "SC", "SD",
                                           "SF", "SH","TS"))
pdm.dat <- pdm.dat[, c(1:2, keep.spp)]
# Ready to start analysis
library(FuzzyQ)
# Compute clustering estimates per spp per yr
yr.spp <- by(pdm.dat[, 3:16], pdm.dat[, "YR"], fuzzyq, rm.absent = FALSE,
              sorting = FALSE)
# Extract species info
yr.spp <- lapply(1:length(yr.spp), function(x) yr.spp[[x]]$spp)
# Extract commonness indices per species
yr.spp <- t(sapply(yr.spp, function(x) x[,'Common.I']))
colnames(yr.spp) <- colnames(pdm.dat[, 3:16])
rownames(yr.spp) <- levels(pdm.dat[, "YR"])
# Plot Commonness Indices per yr
par(mar=c(4, 4.1, 1, 5.5), xpd=TRUE)
spp.col <- c("firebrick2", "limegreen", "chocolate1", "royalblue", 
             "darkorchid", "turquoise",  "greenyellow",
             "turquoise4", "plum2", "tan4", "darkred",
             "darkseagreen1", "gold4", "navyblue")
matplot(yr.spp, type='l', xaxt = 'n', ylab= "Commonness index", col=spp.col,
        lty=1)
abline(h=0.5, lty=2, xpd=FALSE, col="orange3")
axis(1, at = 1:21, labels = rownames(yr.spp), las = 2)
legend(legend = colnames(yr.spp), x="topright", inset=c(-0.2,0), lty = 1,
       lwd = 2, col = spp.col, title='Species', cex=0.8)
#######################################################################
# GAM model of yr variation in Commonness Index of GV Glaucomys volans
library(gamlss)

# Set up a grid of m evenly-spaced points in 21 yrs of time series 
# on which to evaluate a spline
m.grid=100
eval.grid <- seq(from=1,to=21,length.out=m.grid) 
# fit model to data
dat <- as.data.frame(cbind(yr.spp[,'GV'],1:21))
colnames(dat) <- c("comm", "yr")
# family distribution chosen: beta distribution (BE)
fam = BE
mod.GV <- gamlss(comm ~ yr, data= dat, family =fam, trace =FALSE,n.cyc=100)
# Evaluate fit with diagnostic plots
      wp(mod.GV)
      plot(mod.GV)
# Generate a smooth spline to predict 5 yrs ahead
x  <- as.vector(mod.GV$mu.x[,2])
y  <- mod.GV$mu.fv 
sm <- smooth.spline(x, y) 
pred <- predict(sm, x=22:26) 

# Ad hoc function to compute confidence intervals for the model
BS.mod.GV <- function(commInx, eval.grid) {
  data.in <- as.data.frame(cbind(commInx,1:21))
  colnames(data.in) <- c("comm", "yr")
  bsmod <- gamlss(comm ~ yr, data= data.in, family =fam,  trace =FALSE, 
                  n.cyc=40)  
  x   <- as.vector(bsmod$mu.x[,2])
  y   <- bsmod$mu.fv               
  sm  <- smooth.spline(x, y)         # Fit spline to data
  fit <- predict(sm, x=eval.grid)    # Fitted values over grid values
  bspred <- predict(sm, x=22:26)     # Prediction 5yr ahead
  predict.out <- list(bsfit= fit$y, bspred=bspred$y)
}   

# Use fuzzyqBoot to compute bootstrap estimates of commonness
Nboot =1e3 # set no. bootstr replications
set.seed(2020)
BS.yr.spp <- by(pdm.dat[, 3:16], pdm.dat[, "YR"], fuzzyqBoot, N = Nboot,
                     level='spp', maxit=500)
BS.yr.spp <- lapply(BS.yr.spp, function(x) x$fq.rep)
# Extract commonness estimates of GV
BS.commGV <- sapply(1:length(BS.yr.spp), function(x) BS.yr.spp[[x]][ ,"GV"])
BSGV.pred <- apply(BS.commGV, 1, BS.mod.GV, eval.grid=eval.grid) # Adjust a model to each replicate
# Compute 95% percentiles of fitted values
BSGV.pred1 <- sapply(1:length(BSGV.pred), function(x) BSGV.pred[[x]][1])
BSGV.pred1 <- do.call(rbind, BSGV.pred1)
lo.GV.fit <- apply(BSGV.pred1, 2, quantile, probs=0.025)
hi.GV.fit <- apply(BSGV.pred1, 2, quantile, probs=0.975)
# compute 95% percentiles of 5yr predicted values
BSGV.pred2 <- sapply(1:length(BSGV.pred), function(x) BSGV.pred[[x]][2])
BSGV.pred2 <- do.call(rbind, BSGV.pred2)
lo.GV.pred <- apply(BSGV.pred2, 2, quantile, probs=0.025)
hi.GV.pred <- apply(BSGV.pred2, 2, quantile, probs=0.975)
# Plotting Commonness of GV and fitted model
plot(1:21, yr.spp[,'GV'], xaxt = 'n', ylab= "Commonness index",
     xlab="", xlim=c(1,26), cex.axis=0.8, cex.lab=1.4)
bg.col = "papayawhip"
polygon(c(rev(eval.grid),eval.grid), c(rev(hi.GV.fit), lo.GV.fit), border= NA, 
        col=bg.col)
lines(sm, col="darkred", lwd=2)
lines(eval.grid, lo.GV.fit, col="darkred", lty=1)
lines(eval.grid, hi.GV.fit, col="darkred", lty=1)
polygon(c(26:22, 22:26), c(rev(hi.GV.pred), lo.GV.pred), border= NA, 
        col=bg.col)
points(1:21, yr.spp[,'GV'],type="b", lty=3, pch=16, col="royalblue", cex=1.2)
lines(pred, col="darkred", lwd=2, lty=2)
lines(22:26, lo.GV.pred, col="darkred", lty=1)
lines(22:26, hi.GV.pred, col="darkred", lty=1)
abline(h=0.5, lty=2, xpd=FALSE, col="orange3")
axis(1, at=1:26, labels=1979:2004, las=2, cex=0.8)      
mtext("Year", side=1, line=3, cex=1.2)





