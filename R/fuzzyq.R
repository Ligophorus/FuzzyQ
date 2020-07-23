#' Fuzzy Quantification of Common and Rare Species in Ecological Communities
#'
#' Perform fuzzy clustering of each species based on its abundance and occupancy.
#' @param M A Matrix or dataframe of species abundaces (columns). Each row represents a sites or
#'     sampling unit.
#' @param nclus Integer. Indicate the number of clusters into which species should be allocated.
#'     Default is 2 (common and rare).
#' @param diss String. Specifiy the dissimilarity coefficient to be used. Default is "gower".
#'     The other options are "euclidean" and "manhattan".
#' @param rm.absent Logical. Indicate how to treat species absences in the community. If TRUE
#'     absences are considered as structural; absent species are removed from calculations. If
#'     FALSE (the default) absences are considered as random; absent species are included in the
#'     calculations.
#' @param sorting Logical. If TRUE (the default) species are sorted in the output by ascending
#'     silhouette widths within each cluster, else species are arranged in the same order as in the
#'     input matrix or dataframe.
#' @param keep.Diss Logical. Whether or not the species dissimilarity matrix shoudl be returned. The
#'     default is FALSE.
#' @param daisy.args Arguments to be passed to function daisy in package cluster.
#' @param ... Arguments to be passed to function fanny in package cluster.
#' @return A list of class fuzzyq with the abundance occupancy of each species, clustering metrics for each species and the whole community, and the species dissimilary matrix (if diss = TRUE).
#' \describe{
#'   \item{\code{A_O}}{Abundance-occupancy information for each species}
#'   \item{\code{Diss}}{Object of class dist with pairwise dissimilarities among species based on A_O}
#'   \item{\code{spp}}{Clustering metrics per species: Cluster membership (where 0 and 1 denote allocation
#'    to the rare and common category, respectively), Silhouette Widths and Commonness Indices)}
#'   \item{\code{global}}{Community level clustering metrics: Average silhouette widths per cluster and globally,
#'    Mean commonness indices per cluster and Normalized Dunn's coefficient}
#' }
#' @examples
#' data(antsA)
#' # Data set of 46 ant species colelcted in 100 sites.
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)

fuzzyq <- function(M, nclus = 2, diss = "gower", rm.absent = FALSE,
                   sorting = TRUE, keep.Diss = FALSE, daisy.args, ...) {
  if (length(dim(M)) != 2 || !(is.data.frame(M) || is.numeric(M)))
    stop("M is not a dataframe or a numeric matrix.")
  if (rm.absent == TRUE &&
      length(which(colSums(M) == 0)) != 0)  M <- M[, -which(colSums(M) == 0)]
  if (length(dim(M)) != 2) stop("Insufficent data: only one species")
  abund <- colMeans(M, na.rm = TRUE)
  M[M > 0] <- 1
  occ  <- colMeans(M, na.rm = TRUE)
  A_O <- cbind(occ, abund)
  colnames(A_O) <- c("frq.occ", "m.abund")
  if (missing(daisy.args)) D <- cluster::daisy(A_O, metric = diss) else
    D <- do.call(cluster::daisy, c(list(x = A_O, metric = diss), daisy.args))
 # check that there are at sufficient no. of spp for fanny (k >= n/2-1)
  n <- attr(D, "Size")
  if (nclus >= n / 2 - 1) {
    warning("Insufficient number of spp. for fuzzy clustering. NULLs produced")
    sil.w <- NULL
    global <- NULL
    cr.fan <- list(A_O = A_O, spp = sil.w, global = global)
    if (keep.Diss == TRUE) cr.fan <- append(cr.fan, list(Diss = D), 1)
  } else {
    fanclus <- cluster::fanny(D, k = nclus, keep.diss = FALSE,
                              keep.data = FALSE, ...)
    sil.w <- fanclus$silinfo$widths[, c(1, 3)]
    # This part ensures that common and rare are consistently related to
    # clusters. (Rarest sp. is assumed to belong to rare category)
    rar.tag <- names(which(A_O[, 1] == min(A_O[, 1])))[1]
    sil.tag <- which(row.names(sil.w) == rar.tag)
    # 0 = rare; 1 = common
    sil.w[which(sil.w[, 1] == sil.w[sil.tag, 1]), 1] <- 0
    sil.w[which(sil.w[, 1] != sil.w[sil.tag, 1]), 1] <- 1
    # Ensure that "membership" reflects "commonness"
    if (fanclus$membership[which(rownames(fanclus$membership) == rar.tag), 1] <=
        fanclus$membership[which(rownames(fanclus$membership) == rar.tag), 2])
    m <- 1 else m <- 2
    # sort prior to cbind by membership
    sil.w <- sil.w[row.names(fanclus$membership), ]
    sil.w <- cbind(sil.w, fanclus$membership[, m])
    colnames(sil.w)[3] <- "Common.I"
    # sorting options
    if (sorting == TRUE) {
      sil.w <- sil.w[order(sil.w[, 1], sil.w[, 2]), ]
      #order A_O rows as sil.w:
      A_O <- A_O[match(rownames(sil.w), rownames(A_O)), ]
    }
    # compute average sil widths and commoness coeffs
    sil.w <- as.data.frame(sil.w)
    silw <- c(mean(sil.w[sil.w$cluster == 0, 2]),
              mean(sil.w[sil.w$cluster == 1, 2]),
              mean(sil.w[, 2]))
    memb <- c(mean(sil.w[sil.w$cluster == 0, 3]),
              mean(sil.w[sil.w$cluster == 1, 3]))
    global <- c(silw, memb, fanclus$coeff[2])
    names(global) <- c("silw.rar", "silw.com", "silw.all",
                         "commI.rar", "commI.com", "N.Dunn")
    cr.fan <- list(A_O = A_O, spp = sil.w, global = global)
    if (keep.Diss == TRUE) cr.fan <- append(cr.fan, list(Diss = D), 1)
  }
  class(cr.fan) <- append(class(cr.fan), "fuzzyq")
  is.sorted <- sorting
  cr.fan$is.sorted <- is.sorted
  return(cr.fan)
}
