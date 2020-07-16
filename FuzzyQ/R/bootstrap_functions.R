fuzzyqBoot <- function(M, N, level="spp", rm.absent = FALSE,
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
 return(M)
}
#Percentile bootstrap CIs function
.pctCI <- function(M, ci.level) {
  P1 <- ci.level / 2
  P2 <- 1 - P1
  bs.ci <- apply(M, 2, quantile, probs = c(P1, P2), na.rm = TRUE)
  return(bs.ci)
}
#bias-corrected CI function
.bcCI <- function(bt_x, fq.obj, N = N, ci.level) {
  z <- qnorm(c(ci.level / 2, 1 - ci.level / 2)) # Std. norm. limits
  b <- qnorm((sum(bt_x > fq.obj) + sum(bt_x == fq.obj) / 2) / N)
  p <- pnorm(z - 2 * b) # bias-correct & convert to proportions
  bs.ci <- quantile(bt_x, probs = p)
  return(bs.ci)
}
#bias-corrected & accelerated CI function
.bcaCI <- function(bt_x, fq.obj, N = N, ci.level) {
  # estimate acceleration constant
  n1 <- N - 1
  obsn <- fq.obj * N
  pv <- sapply(1:N, function(x) sapply(bt_x[-x], mean))
  pv <- rowMeans(pv)
  pv <- obsn - n1 * pv
  je <- mean(pv) - pv
  a <- sum(je^3) / (6 * sum(je^2))^(3 / 2)
  b <- qnorm((sum(bt_x > fq.obj) + sum(bt_x == fq.obj) / 2) / N) #bias
  z <- qnorm(c(ci.level / 2, 1 - ci.level / 2)) # Std. norm. limits
  p <- pnorm((z - b) / (1 - a * (z - b)) - b) # correct & convert to proportions
  bs.ci <- quantile(bt_x, probs = p)
  return(bs.ci)
}

fuzzyqCI <- function(M, fq.obj, level="spp", method = "pct", c.level = 0.95) {
  if (!(is.data.frame(M) || is.numeric(M)))
    stop("M is not a dataframe or a numeric matrix.")
  if ("fuzzyq" %in% class(fq.obj) == FALSE)
    stop("fq.obj is not a fuzzyq object.")
  if (level %in% c("spp", "global") == FALSE)
    stop("Wrong community level specification.
         Valid choices are 'spp' or 'global'")
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
    if (missing(fq.obj) && method != "pct")
      stop("fq.obj required for BC or BCa confidence intervals")
    N <- nrow(M)
    if (level == "spp") {
      fq.obj <- fq.obj$spp
      fq.obj <- fq.obj[match(colnames(M), rownames(fq.obj)), ]
      fq.obj <- fq.obj[, 3]} else fq.obj <- fq.obj$global
    if (length(fq.obj) != ncol(M))
      stop("Number of observed values and bootstrap estimates mismatch")
    F2apply <- ifelse(method == "bc", .bcCI, .bcaCI)
    BS.CI <- sapply(seq_len(ncol(M)), function(x) F2apply(M[, x],
                      fq.obj = fq.obj[x], N = N, ci.level = c.level))
    colnames(BS.CI) <- colnames(M)
    rownames(BS.CI) <- c("Lower", "Upper")
    return(BS.CI)
  }
}
