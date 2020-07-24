#' Apply Fuzzy Quantification of Common and Rare Species to Bootstrap Replicates
#'
#' Produce N replicates of the original site by species matrix or dataframe by taking bootstrap
#'     samples of sites (rows) and apply \code{fuzzyq} to each replicate.
#' @param M A matrix or dataframe of species abundaces (columns). Each row represents a site.
#' @param N Integer. Number of bootsrap replicates desired. Default is 1000.
#' @param level String. Specifiy the type of metrics to be computed for each bootstrap replicate.
#'     Either \code{"spp"} or \code{"global"}, corresponding to species or community-level metrics,
#'     respectively.
#' @param rm.absent Logical. Whether or not absent species are to be removed from the calculations.
#' @param sorting Logical. If \code{TRUE} (the default) species are sorted in the output by ascending
#'     silhouette widths within each cluster, else species are arranged in the same order as in the
#'     input matrix or dataframe.
#' @param daisy.args Arguments to be passed to function \code{daisy} in package \code{cluster}.
#' @param ... Arguments to be passed to function \code{fanny} in package \code{cluster}.
#' @return A list consisting of the following:
#' \describe{
#'   \item \code{bs.rep}{matrix of estimated metrics. Replicates are arranged in rows. If \code{level = "spp"},
#'     columns represent estimates of Commonness Indices per species. If \code{level = "global"}, columns
#'     represent estimates of Community level clustering metrics: Average silhouette widths per cluster
#'     and globally, Mean commonness indices per cluster and Normalized Dunn's coefficient.}
#'  \item \code{level}{Indicates whether the estimates are taken at species (\code{"spp"}) or
#'     community level (\code{"global"}).}
#' }
#' @examples
#' data(antsA)
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)
#' # Estimates of species Commonness Indices by species of 1,000 bootstrap replicates:
#' BS.FQAnts <- fuzzyqBoot (antsA, N = 1e3, level='spp')
#' # Estimates of global metrics of 1,000 boostrap replicates:
#' BS.global <- fuzzyqBoot (antsA, N = 1e3, level='global')
fuzzyqBoot <- function(M, N = 1e3, level="spp", rm.absent = FALSE,
                         daisy.args, ...) {
  if (length(dim(M)) != 2 || !(is.data.frame(M) || is.numeric(M)))
    stop("M is not a vector, a dataframe or a numeric matrix.")
  if (level %in% c("spp", "global") == FALSE)
    stop("Wrong community level specification.
         Valid choices are 'spp' or 'global'")
  if (rm.absent == TRUE &&
      length(which(colSums(M) == 0)) != 0)  M <- M[, -which(colSums(M) == 0)]
  if (length(dim(M)) != 2) stop("Insufficent data: only one species")
  # Bootstrap matrix by sites (rows)
  bootSxSP <- function(M) {
    boot_M <- function(M) {
      a <- sample(seq_len(nrow(M)), nrow(M), replace = TRUE)
      M <- M[a, ] # data bootstrapped by rows
      return(M)
    }
    M <- replicate(N, boot_M(M), simplify = FALSE)
    return(M)
  }
  spp.names <- colnames(M)
  M <- bootSxSP(M)
  if (missing(daisy.args))
    M <- suppressWarnings(lapply(M, fuzzyq, rm.absent = FALSE,
                                 keep.Diss = FALSE, sorting = FALSE, ...)) else
    M <- suppressWarnings(lapply(M, fuzzyq, rm.absent = FALSE, sorting = FALSE,
                                 daisy.args = dasiy.args, ...))
  null.control <- sapply(M, function(x) is.null(x$spp))
  sum.nulls <- sum(null.control, na.rm = TRUE)
  if (sum.nulls == N) stop("All replicates contained too few species
                           for fuzzy clustering.")
  if (sum.nulls > 0 && sum.nulls < N) {
    null.replicates <- which(null.control == TRUE)
    M <- M[-null.replicates]
    options(warning.length = 2000L)
    warning(paste(c("Replicates", null.replicates,
                  "contained too few species for fuzzy clustering and
                   have been removed."), collapse = " "))
  }
  if (level == "spp") {
     M <- t(sapply(M, function(x) cbind(x$spp[, "Common.I"])))
     colnames(M) <- spp.names
  } else {
     M <- t(sapply(M, function(x) x$global))
  }
 fq.bs <- list(fq.rep = M, level = level)
 return(fq.bs)
}
#Percentile bootstrap CIs function
.pctCI <- function(M, ci.level) {
  P1 <- ci.level / 2
  P2 <- 1 - P1
  bs.ci <- apply(M, 2, quantile, probs = c(P1, P2), na.rm = TRUE)
  return(bs.ci)
}
#bias-corrected CI function
.bcCI <- function(bt_x, fq, N = N, ci.level) {
  z <- qnorm(c(ci.level / 2, 1 - ci.level / 2)) # Std. norm. limits
  b <- qnorm((sum(bt_x > fq) + sum(bt_x == fq) / 2) / N)
  p <- pnorm(z - 2 * b) # bias-correct & convert to proportions
  bs.ci <- quantile(bt_x, probs = p)
  return(bs.ci)
}
#bias-corrected & accelerated CI function
.bcaCI <- function(bt_x, fq, N = N, ci.level) {
  # estimate acceleration constant
  n1 <- N - 1
  obsn <- fq * N
  pv <- sapply(1:N, function(x) sapply(bt_x[-x], mean))
  pv <- rowMeans(pv)
  pv <- obsn - n1 * pv
  je <- mean(pv) - pv
  a <- sum(je^3) / (6 * sum(je^2))^(3 / 2)
  b <- qnorm((sum(bt_x > fq) + sum(bt_x == fq) / 2) / N) #bias
  z <- qnorm(c(ci.level / 2, 1 - ci.level / 2)) # Std. norm. limits
  p <- pnorm((z - b) / (1 - a * (z - b)) - b) # correct & convert to proportions
  bs.ci <- quantile(bt_x, probs = p)
  return(bs.ci)
}
#' Compute Confidence Intervals of Clustering Metrics
#'
#' Computes confidence intervals of clustering metrics based on the bootstrap replicates produced by
#'     \code{fuzzyqBoot}.
#' @param fq.bs A list returned by \code{fuzzyqBoot}.
#' @param fq A list of class \code{fuzzyq} returned by \code{FuzzyQ::fuzzyq}. Required only if
#'     \code{method = "bc"} or \code{method = "bca"}.
#' @param method String. Specify the method to compute confidence intervals. Any of the following:
#'     "pct" (percentile), "bc" (bias corrected), "bca" (bias corrected and accelerated).
#' @param c.level Number within [0,1]. Specify the confidence interval level. Default is 0.95.
#' @return A matrix with upper and lower confidence interval limits of clustering metrics.
#' @examples
#' data(antsA)
#' FQAnts <- fuzzyq(antsA, sorting = TRUE)
#' # Estimates of species Commonness Indices by species of 1,000 bootstrap replicates:
#' BS.FQAnts <- fuzzyqBoot (antsA, N = 1e3, level='spp')
#' # Estimates of global metrics of 1,000 boostrap replicates:
#' BS.global <- fuzzyqBoot (antsA, N = 1e3, level='global')
#' ############# PENDING #################
fuzzyqCI <- function(fq.bs, fq = NULL, method = "pct", c.level = 0.95) {
  M <- fq.bs$fq.rep
  if (!(is.data.frame(M) || is.numeric(M)))
    stop("M is not a dataframe or a numeric matrix.")
  if (!missing(fq) && "fuzzyq" %in% class(fq) == FALSE)
    stop("fq is not a fuzzyq object.")
  if (method %in% c("pct", "bc", "bca") == FALSE)
    stop("Wrong confidence interval specification.
         Valid choices are 'pct','bc' or 'bca'")
  c.level <- 1 - c.level
  if (method == "pct") {
    BS.CI <- .pctCI(M, ci.level = c.level)
    colnames(BS.CI) <- colnames(M)
    rownames(BS.CI) <- c("Lower", "Upper")
    return(BS.CI)
  } else {
    if (missing(fq) && method != "pct")
      stop("fq required for BC or BCa confidence intervals")
    N <- nrow(M)
    if (fq.bs$level == "spp") {
      fq <- fq$spp
      fq <- fq[match(colnames(M), rownames(fq)), ]
      fq <- fq[, 3]} else fq <- fq$global
    if (length(fq) != ncol(M))
      stop("Number of observed values and bootstrap estimates mismatch")
    F2apply <- ifelse(method == "bc", .bcCI, .bcaCI)
    BS.CI <- sapply(seq_len(ncol(M)), function(x) F2apply(M[, x],
                      fq = fq[x], N = N, ci.level = c.level))
    colnames(BS.CI) <- colnames(M)
    rownames(BS.CI) <- c("Lower", "Upper")
    return(BS.CI)
  }
}
