# x es la matrix de abundancias
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