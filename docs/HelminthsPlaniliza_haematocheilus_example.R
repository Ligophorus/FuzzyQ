# Comparison of common-rare helminth species AORs of 
# Planiliza haematocheilus
# Japan Sea (native area) - Sea of Azov (introduced area)
library(FuzzyQ)
data(Azov) # helminth sp. abundance dataset 12 surveys Azov
data(Japan) # idem 7 surveys Japan

# Silhouette barplots
# Get silhouette widths per sp. per survey in both areas:
fuzzyq.azov <- by(Azov[, -1], Azov[, "sample"], fuzzyq, rm.absent = FALSE)
sppsilw.azov <- lapply(fuzzyq.azov, function(x) x$spp)

fuzzyq.japan <- by(Japan[, -1], Japan[, "sample"], fuzzyq, rm.absent = FALSE)
sppsilw.japan <- lapply(fuzzyq.japan, function(x) x$spp)

# Ad hoc function to plot silhouette widths
sil.barplot <- function(M, maint, col.bars=c("brown2", "turquoise3"),...) {
  barplot(M[,2], col=col.bars[M[,1]+1], main=maint, names.arg= rownames(M),...)
}
# Plot silhouette widths
# Azov
par(mfcol=c(4,3), mar = c(2.5, 1, 1.5, 1), oma=c(0,2.5,0,0), cex.axis=0.8,
    las=2, mgp=c(3,0.5,0), font.main=1)
mapply(sil.barplot, sppsilw.azov, levels(Azov$sample))
mtext("Silouhette widths", side =2, line= 1, outer=TRUE,
      at= c(0.13,0.381, 0.633,0.875), las=3, cex=0.6)
# Japan
par(mfcol=c(4,2), mar = c(2.5, 1, 1.5, 1), oma=c(0,2.5,0,0), cex.axis=0.8,
    las=2, mgp=c(3,0.5,0), font.main=1)
mapply(sil.barplot, sppsilw.japan, levels(Japan$sample))
mtext("Silouhette widths", side =2, line= 1, outer=TRUE,
      at= c(0.125,0.381, 0.63,0.882), las=3, cex=0.8)

# Collate global estimates per survey
Results.azov <- t(sapply(fuzzyq.azov, function(x) x$global))
Results.japan <- t(sapply(fuzzyq.japan, function(x) x$global))
## join 
Results.all <- as.data.frame(rbind(Results.azov,Results.japan))
Results.all$Area <- as.factor(c(rep("Azov", nrow(Results.azov)), 
                               rep("Japan", nrow(Results.japan))))

# Test differences in parameters between areas 
# Global average silhouette width:
    wilcox.test(silw.all ~ Area, data = Results.all)
# Common species average silhouette width:
    wilcox.test(silw.com ~ Area, data = Results.all)
# Rare species average silhouette width
    wilcox.test(silw.rar ~ Area, data = Results.all)
# Normalized Dunn's index:   
wilcox.test(N.Dunn ~ Area, data = Results.all)
# Common species average Commonness index:
    wilcox.test(commI.com ~ Area, data = Results.all)
# Rare species average Commonness index:
    wilcox.test(commI.rar ~ Area, data = Results.all)

# Compute 95% confidence intervals of global parameters by bootstrap
Nboot =1e3 # set no. bootstrap replicates
set.seed(2020) # randomization seed - optional

BSazov <- by(Azov[, -1], Azov[, "sample"], fuzzyqBoot, N = Nboot,
             level= "global")
# 
BSazov <- lapply(seq_len(length(BSazov)),
                 function(x) fuzzyqCI(BSazov[[x]],
                             fuzzyq.azov[[x]], method = "bca"))
                                            
#
BSjapan <- by(Japan[, -1], Japan[, "sample"], fuzzyqBoot,
              N = Nboot, level = "global")
           # increase number iterations (maxit) to reach convergence 
BSjapan <- lapply(seq_len(length(BSjapan)), function(x) fuzzyqCI(BSjapan[[x]],
                  fuzzyq.japan[[x]], method = "bca"))
# Collate results from both areas
BSall <- c(BSazov, BSjapan)

# To plot estimates +- 95% CIs:
# Prepare data for boxplots
silw.all <- cbind(sapply(BSall, function(x) x[, "silw.all"]))
silw.com <- cbind(sapply(BSall, function(x) x[, "silw.com"]))
silw.rar <- cbind(sapply(BSall, function(x) x[, "silw.rar"]))
N.Dunn <- cbind(sapply(BSall, function(x) x[, "N.Dunn"]))
commI.com <- cbind(sapply(BSall, function(x) x[, "commI.com"]))
commI.rar <- cbind(sapply(BSall, function(x) x[, "commI.rar"]))

# 2 colors by area to plot results
col.AJ <- c("darkorange2", "slategray")
col.tag <- col.AJ[Results.all$Area]
# Plot estimates +- 95% confidence intervals
par(mfcol=c(3,2), mar = c(0.2, 6, 0.5, 1), oma=c(3.5,0.25,0,0), cex.axis=0.8)
plot(Results.all$silw.all, col=col.tag, ylab="Average silhouette width",
    xlab="", xaxt='n', pch=16, ylim=c(min(silw.all),max(silw.all)))
  arrows(1:nrow(Results.all),silw.all[1,],1:nrow(Results.all),silw.all[2,],
       length=0, angle=90, code=3, col=col.tag)
plot(Results.all$silw.com, col=col.tag, ylab="Average silhouette width",
    xlab="", xaxt='n', pch=16, ylim=c(min(silw.com),max(silw.com)))
  arrows(1:nrow(Results.all),silw.com[1,],1:nrow(Results.all),silw.com[2,],
       length=0, angle=90, code=3, col=col.tag)
plot(Results.all$silw.rar, col=col.tag, ylab="Average silhouette width",
    xlab="", xaxt='n', pch=16, ylim=c(min(silw.rar), max(silw.rar)))
  arrows(1:nrow(Results.all),silw.rar[1,],1:nrow(Results.all),silw.rar[2,],
         length=0, angle=90, code=3, col=col.tag)  
  axis(side=1, at=1:nrow(Results.all),rownames(Results.all), las=2)
plot(Results.all$N.Dunn, col=col.tag, ylab="Normalized Dunn's coefficient",
    xlab="", xaxt='n', pch=16, ylim=c(min(N.Dunn),max(N.Dunn)))
  arrows(1:nrow(Results.all), N.Dunn[1,],1:nrow(Results.all),N.Dunn[2,],
         length=0, angle=90, code=3, col=col.tag)
legend(16, 0.9, c("Azov", "Japan"), col = col.AJ, pch = 16)
plot(Results.all$commI.com, col=col.tag, ylab="Average commonness index",
       xlab="", xaxt='n', pch=16, ylim=c(min(commI.com),max(commI.com)))
  arrows(1:nrow(Results.all), commI.com[1,],1:nrow(Results.all),commI.com[2,],
         length=0, angle=90, code=3, col=col.tag)
plot(Results.all$commI.rar, col=col.tag, ylab="Average commonness index",
       xlab="", xaxt='n', pch=16, ylim=c(min(commI.rar),max(commI.rar)))
  arrows(1:nrow(Results.all), commI.rar[1,],1:nrow(Results.all),commI.rar[2,],
         length=0, angle=90, code=3, col=col.tag)
  axis(side=1, at=1:nrow(Results.all),rownames(Results.all), las=2)
mtext("Survey", side = 1, line = 2, outer = TRUE, at = c(0.30, 0.80), las = 1)

