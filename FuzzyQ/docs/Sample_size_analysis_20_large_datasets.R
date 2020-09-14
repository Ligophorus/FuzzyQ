library(FuzzyQ)
# Two ad hoc functions:
# Neffect computes fuzzyq global metrics with n sites chosen randomly
# Parameters
#   x: number of sites taken at random to be dropped
#   M: site-by-species abundance matrix
#   Maxit: Maximum number of iterations (used for fuzzy clustering)
# Return
#   Community-level metrics computed by fuzzyq
Neffect <- function (x, M, Maxit) { 
  x <- sample(1:nrow(M), x)
  M <- M[-x, ] # remove x rows
  M <- fuzzyq(M, rm.absent = TRUE, maxit = Maxit, sorting = FALSE)
  return(M$global)
}

# PlotWrapup plots results of Neffect with increasing numbers of sites
# Parameters
#   M: site-by-species abundance matrix
#   maxI: Maximum number of iterations (used for fuzzy clustering)
#   drop.lim: Minimum number of sites to be analyzed
# Return
#   Line plot showing the variation of community level parameters with the
#   number of sites, increasing from drop.lim to the total number of sites
PlotWrapUp <- function (M, maxI = 500, drop.lim = 10) {
  nsites <- nrow(M)
  drop <- 1:(nsites - drop.lim)
  M <- t(sapply(drop, Neffect, M, Maxit = maxI))
  M <- M[nrow(M):1, ]  # reverse results for plotting
  matplot(M, type = "l", col = Colors, lty = 1, ylab= "", xaxt = "n")
  # Re-scale x axis to place drop.lim at 1
  if(nsites <= 200) axis(1, at = seq(1, nsites-drop.lim, 20), 
                         labels = seq(drop.lim, nsites, 20)) else
                    axis(1, at = seq(1, nsites-drop.lim, 50), 
                         labels = seq(drop.lim, nsites, 50)) 
  abline(v = 30 - drop.lim + 1, lty = 3) # draw lines at 30 & 50
  abline(v = 50 - drop.lim + 1, lty = 3)
}

# Colors to plot each community-level metric
Colors <- c("darkgreen", "darkred", "darkblue", "chartreuse3", "brown2",
            "steelblue")
# set.seed required if one wishes to make results reproducible.
set.seed(2020)

################################
#### Read datasets #############
################################
# 20 databases in which no. sites >=87
# of these, 16 from  Jeliazkov et al. (2020) 
# Scientific Data, 7, 6. doi: 10.1038/s41597-019-0344-7

load("CESTES.RData") # Downloaded from
# https://idata.idiv.de/ddm/Data/ShowData/286
# doi: 10.25829/idiv.286-21-2695

brindamour <- LSmatred$BrindAmour2011a$comm # relative abundance
pavoine <- LSmatred$Pavoine2011$comm
diaz <- LSmatred$Diaz2008$comm
lowe <- LSmatred$Lowe2018a$comm
bartonova <- LSmatred$Bartonova2016$comm # abundance class (1-8)
jeliazkov <- LSmatred$Jeliazkov2014$comm
charbonier <- LSmatred$Charbonnier2016a$comm
barbaro <- LSmatred$Barbaro2009a$comm
jeliazkov2 <- LSmatred$Jeliazkov2013$comm
barbaro2 <- LSmatred$Barbaro2009b$comm # relative abundance
charbonier2 <- LSmatred$Charbonnier2016b$comm
fried <- LSmatred$Fried2012$comm # densities
goncalves <- LSmatred$Goncalves2014a$comm
goncalves2 <- LSmatred$Goncalves2014b$comm
chmura <- LSmatred$Chmura2016$comm
ribera <- LSmatred$Ribera2001$comm

# Four additional datasets from Calatayud et al. (2020)
# Nature Ecology & Evolution, 1–6. doi:10.1038/s41559-019-1053-5
load("abundace_matrices.RData") # Downloaded from
# https://doi.org/10.6084/m9.figshare.9906092

antsA <- abun.mat$ants_data_Xavi_Darwin_A
antsB <- abun.mat$ants_data_Xavi_Darwin_B
antsC <- abun.mat$ants_data_Xavi_Darwin_C
antsD <- abun.mat$ants_data_Xavi_Darwin_D

# FIGURE 5 in publication. Plot 8 datasets
par(mfrow=c(4,2), mar=c(2,3,1,1), oma=c(4,1.5,0.5,0))
PlotWrapUp(brindamour, drop.lim = 20)
PlotWrapUp(pavoine)
PlotWrapUp(jeliazkov, drop.lim = 20)
PlotWrapUp(barbaro, maxI = 1000)
PlotWrapUp(chmura)
PlotWrapUp(ribera)
PlotWrapUp(goncalves, drop.lim = 20)
PlotWrapUp(antsA)
mtext("Number of sites", side=1, at=c(0.275,0.775), outer=TRUE, line=0.5, 
      cex = 0.9)
mtext("Community indices", side=2, at=c(0.14,0.38,0.63,0.88), outer=TRUE, 
      line=-0.5, cex = 0.9)
mtext(c("a)", "c)","e)","g)"), side=2, at=c(0.99,0.74,0.49,0.24), outer=TRUE,
      line=0.1, las=2)
mtext(c("b)", "d)","f)","h)"), side=2, at=c(0.99,0.74,0.49,0.24), outer=TRUE,
      line=-27.5, las=2)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# expressions to pass to legend:
SR <- expression(italic(bar(S)['R']))
SC <- expression(italic(bar(S)['C']))
S <- expression(italic(bar(S)))
CR <- expression(italic(bar(C)['R']))
CC <- expression(italic(bar(C)['C']))
D <- expression(italic("D'"))
legend(x = "bottom",inset = 0, legend = c(SR, SC, S, CR, CC, D), 
       col= Colors, cex=1, horiz = TRUE, lwd=2, xpd = TRUE, bty = "n")

# Figure S3 in Supporting Information. Plot 12 datasets
par(mfrow=c(6,2), mar=c(2,3,1,1), oma=c(4,1.5,0.5,0))
PlotWrapUp(diaz)
PlotWrapUp(lowe)
PlotWrapUp(bartonova)
PlotWrapUp(jeliazkov2)
PlotWrapUp(barbaro2)
PlotWrapUp(charbonier)
PlotWrapUp(charbonier2)
PlotWrapUp(fried)
PlotWrapUp(goncalves2, drop.lim = 20)
PlotWrapUp(antsB)
PlotWrapUp(antsC)
PlotWrapUp(antsD)
mtext("Number of sites", side=1, at=c(0.275,0.775), outer=TRUE, line=0.5, 
      cex = 0.6)
mtext("Community indices", side=2, at=c(0.09, 0.26, 0.425, 0.595, 0.755, 0.925),
      outer = TRUE, line= -0.5, cex = 0.6)
mtext(c("a)", "c)","e)","g)", "i)", "k)"), side=2,
      at=c(0.99, 0.82, 0.666, 0.5, 0.333, 0.166), outer = TRUE,
      line = 0.1, las = 2, cex = 0.6)
mtext(c("b)", "d)","f)","h)", "j)", "l)"), side=2,
      at=c(0.99, 0.82, 0.666, 0.5, 0.333, 0.166), outer=TRUE,
      line=-27.5, las=2, cex = 0.6)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# expressions to pass to legend:
SR <- expression(italic(bar(S)['R']))
SC <- expression(italic(bar(S)['C']))
S <- expression(italic(bar(S)))
CR <- expression(italic(bar(C)['R']))
CC <- expression(italic(bar(C)['C']))
D <- expression(italic("D'"))
legend(x = "bottom",inset = 0, legend = c(SR, SC, S, CR, CC, D), 
       col= Colors, cex=1, horiz = TRUE, lwd=2, xpd = TRUE, bty = "n")
