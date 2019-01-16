#compute bootstrap for each sample
# codtrace data adult helminths
library(cluster)
# read data 
X <- read.table("codtrace.txt", header=TRUE)
# Samp should be a factor -> perhaps not needed but just in case
X$Samp <- factor(X$Samp)
#new var combines Loc and Sampl
X$LxS <- factor(paste(X$Loc,X$Samp, sep="."))
# vector with new var levels
lxs <- levels(X$LxS)
# 
#### Function that boostraps data matrix 
boot_M <- function (M) {
  a <- sample(1:nrow(M), nrow(M), replace=TRUE)
  M <- M[a, ] # data bootsrapped by rows
} 
  
#### Function that computes
#### mean abundance and prevalence for clustering
cs_clus <- function (M) {  
  # mean abundance
  abund <- colMeans(M, na.rm=TRUE)
  # prevalence
  M[M > 0] <- 1 # convert to presence/ausence 
  prev  <- colMeans(M, na.rm=TRUE)
  # bind columns
  clus <- cbind(prev,abund)
  # remove absent species
  clus <- clus[-which(rowSums(clus)==0), ]
  return(clus)
}

#### Function that computes PAM 
#### It saves medioids and silouhette widths
#### inmput matrix (M - prevalence abundance in columns, species in rows) 
cs_pam <- function (clus) {
      # 1st controls that there are at least 3 spp.
  if (nrow(clus)<3) {
    medioid <- NULL
    sil.w <- NULL
    warning("Number of spp. < 3. NULL value produced")
  } else {
      # merge 
      # compute PAM with gower dist may include in future euclidean etc.
  gow <- cluster::daisy(clus, metric="gower")
  pamclus <- pam(gow, 2, stand=TRUE)
      # produce data
  sil.w <- pamclus$silinfo$widths[ ,c(1,3)]
  #sil.w <- sil.w[order(row.names(sil.w)),]
      # This part ensures that core and satellite are consistently related to 
      # clusters. (Rarest sp. is assumed to be satellite)
  sat.tag <- names(which(clus[,1]==min(clus[,1])))[1]
  sil.tag <- which(row.names(sil.w)==sat.tag)
  #0 = sat; 1 = core
  sil.w[which(sil.w[,1] == sil.w[sil.tag,1]), 1] <- 0
  sil.w[which(sil.w[,1] != sil.w[sil.tag,1]), 1] <- 1
  # order rownames alphabetically
  sil.w <- sil.w[ order(row.names(sil.w)), ]
  }
  return (sil.w)
}
#### Function computes FANNY 
#### Saves silouhette widths and memebrship coeff
#### inmput matrix (M - prevalence abundance in columns, species in rows) 
cs_fan <- function (clus, ...) {
  # controls that there are at least 6 spp. (k>= n/2-1)
  if (nrow(clus)<6) {
    sil.w <- NULL
    warning("Number of spp. < 6. NULL value produced")
  } else {
  # merge 
  # compute FANNY with gower dist may include in future euclidean etc.
  gow <- daisy(clus, metric="gower")
  # compute FAN 
  fanclus <- fanny(gow, 2, stand=TRUE, keep.diss=FALSE, keep.data=FALSE,
                   ...)
  # produce data
  sil.w <- fanclus$silinfo$widths[ ,c(1,3)]
  
  # This part ensures that core and satellite are consistently related to 
  # clusters. (Rarest sp. is assumed to be satellite)
  sat.tag <- names(which(clus[,1]==min(clus[,1])))[1]
  sil.tag <- which(row.names(sil.w)==sat.tag)
  # 0 = sat; 1 = core
  sil.w[which(sil.w[,1] == sil.w[sil.tag,1]), 1] <- 0
  sil.w[which(sil.w[,1] != sil.w[sil.tag,1]), 1] <- 1
  
  # Ensure that "membership" reflects "coreness"
  if (fanclus$membership[which(rownames(fanclus$membership)==sat.tag),1]<=
      fanclus$membership[which(rownames(fanclus$membership)==sat.tag),2])
    m <- 1 else m <- 2
  # sort prior to cbind with membership
  sil.w <- sil.w[row.names(fanclus$membership),]
  sil.w <- cbind(sil.w, fanclus$membership[,m])
  colnames(sil.w)[3] <- "Memb"
  # order rownames alphabetically
  sil.w <- sil.w[ order(row.names(sil.w)), ]
  }
  return (sil.w)
}

#  Functions to remove NULLs in list
#  posted by Josh O'Brien in StacK Overflow:
 # 1. A helper function that tests whether an object is either NULL _or_ 
 # a list of NULLs
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
 
 # 2. Recursively step down into list, removing all such objects 
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}
# set number of bootsrap replicates
Nboot=1e3

# Celtic.1
x <- X[which(X$LxS==lxs[4]), 7:30]
#compute with real data
abund <- colMeans(x, na.rm=TRUE)
# prevalence
x[x > 0] <- 1 # convert to presence/ausence 
prev  <- colMeans(x, na.rm=TRUE)
# bind columns
clus <- cbind(prev,abund)
# remove absent species
clus <- clus[-which(rowSums(clus)==0) ]
C1.gow <- daisy(clus, metric="gower")
# compute FAN 
C1.fanclus <- fanny(C1.gow, 2)
plot(C1.fanclus)


# abundance
#Celtic.1 <- cs_clus(x) 
#Celtic.1p <- cs_pam(Celtic.1) 
#Celtic.1f <- cs_fan(Celtic.1)

#compute bootstrap replicates
#set.seed(2)
Bclus <- replicate (Nboot, boot_M(x), simplify = FALSE)
Bclus <- lapply(Bclus, cs_clus)
Bpam <- lapply (Bclus, cs_pam)
Bfan <- lapply (Bclus, cs_fan)

#remove NULLs in list. Applicable if warnings indicate NULL values
#### CHECK THIS!!!!!!!
Bpam <- rmNullObs(Bpam)
Bpam[sapply(Bpam, is.null)] <- NULL

Bfan <- rmNullObs(Bfan)
Bfan[sapply(Bpam, is.null)] <- NULL

# matrix of cluster membership PAM
NA.matrix <- matrix(rep(NA, ncol(x)*length(Bpam)), length(Bpam), ncol(x))
colnames(NA.matrix) <- sort(colnames(x))
# 
# compute coreness pam
core.mp <- NA.matrix
for (i in 1:length(Bpam)) {
  core.mp[i, which(colnames(core.mp) %in% rownames(Bpam[[i]]))] <- Bpam[[i]][,1] 
}
#remove all NA columns
core.mp <- core.mp[ , ! apply( core.mp , 2 , function(x) all(is.na(x)) ) ]
#compute mean and confidence intervals
core.pam <- colMeans(core.mp, na.rm=TRUE)
sd.pam <- apply(core.mp, 2, sd, na.rm=TRUE)
ci.pam <- apply(core.mp, 2, quantile, probs=c(0.025,0.5, 0.975), na.rm=TRUE)
#boxplot(core.mp, las=2)

# compute coreness fanny
core.mf <- NA.matrix
for (i in 1:length(Bfan)) {
  core.mf[i, which(colnames(core.mf) %in% rownames(Bfan[[i]]))] <- Bfan[[i]][,1] 
} 
#remove all NA columns
core.mf <- core.mf[ , ! apply( core.mf , 2 , function(x) all(is.na(x)) ) ]
core.fan <- colMeans(core.mf, na.rm=TRUE) 

# compute member coreness fanny
core.Mf <- NA.matrix
for (i in 1:length(Bfan)) {
  core.Mf[i, which(colnames(core.mf) %in% rownames(Bfan[[i]]))] <- Bfan[[i]][,3] 
} 
#remove all NA columns
core.Mf <- core.Mf[ , ! apply( core.Mf , 2 , function(x) all(is.na(x)) ) ]
core.Mfan <- colMeans(core.Mf, na.rm=TRUE) 
#apply(core.Mf, 2, quantile, probs=c(0.025,0.5, 0.975), na.rm=TRUE)


core.Mf.list <- lapply(seq_len(ncol(core.Mf)), function(x) core.Mf[,x])
names(core.Mf.list) <- colnames(core.Mf)

stripchart(core.Mf.list, vertical=TRUE, las=2, method="jitter", pch="·",
           col="blue", xlab="Species", ylab="Coreness", cex.lab=1.1, cex=0.8,
           cex.axis=0.8, jitter=0.2)
points(core.Mfan, pch=20, col="red")


# Select survey - only abundace columns 
for(i in 4:length(lxs)) {
x <- X[which(X$LxS==lxs[i]), 7:30]
# abundance
abund <- colMeans(x, na.rm=TRUE)
# prevalence
x[x > 0] <- 1 # presence/aubsence 
prev  <- colMeans(x, na.rm=TRUE)
# bind data
clus <- cbind(prev,abund)
# remove 0 abundant spp.
clus <- clus[-which(rowSums(clus)==0), ]
#clus[,2] <- log(clus[,2]+0.01)
# run pam & fanny with gower distances
gow <- daisy(clus, metric="gower")
pamclus <- pam(gow, 2, stand=TRUE)
plot(pamclus, which.plots =2)
}
