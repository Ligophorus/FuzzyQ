sortClus <- function(M, fq.obj) {
  if (length(dim(M)) != 2 || !(is.data.frame(M) || is.numeric(M)))
    stop("M is not a dataframe or a numeric matrix.")
  if ("fuzzyq" %in% class(fq.obj) == FALSE) stop("fq.obj is not a fuzzyq
                                                  object.")
  if (fq.obj$is.sorted == FALSE) stop("Common-rare species are not sorted in M.
                                      Run fuzzyq with sorting = TRUE")
  M <- M[, match(rownames(fq.obj$spp), colnames(M))]
  return(M)
}

AOplot <- function(fq, col.rc = c("red", "blue"), opacity = 0.1,
                  log.x = FALSE, log.y = FALSE,
                  xLab = "Fraction of sites occupied", yLab = "Mean abundance",
                    ...) {
  if ("fuzzyq" %in% class(fq) == FALSE) stop("fq is not a fuzzyq object.")
  if (fq$is.sorted == FALSE) stop("Common-rare species are not sorted in M.
                              Run fuzzyq with sorting= TRUE")
  if (length(col.rc) != 2) stop("Wrong col.rc format. Specify two colors")
  # check color format
   is.color <- function(x) {
    if (is.numeric(x)) return(x > 0 & (x %% 1 == 0))
    if (any(grepl("^[0-9]+$", x))) {
      x[grepl("^[0-9]+$", x)] <- palette()[(as.numeric(x[grepl("^[0-9]+$",
                                            x)]) - 1) %% 8 + 1]
    }
    y <- grepl("^\\#[a-fA-F0-9]{6}$", x) | grepl("^\\#[a-fA-F0-9]{8}$",
                                                 x) | (x %in% colors())
    return(y)
   }
  if (all(is.color(col.rc)) == FALSE) stop("Wrong col.rc specification.
                                           Use a valid color call.")
      x <- fq$A_O[, 1]
      y <- fq$A_O[, 2]
      if (log.x == TRUE) {
         x <- log10(x)
         xLab <- bquote(paste(log[10]~"("~.(xLab)~")"))
         }
      if (log.y == TRUE) {
         y <- log10(y)
         yLab <- bquote(paste(log[10]~"("~.(yLab)~")"))
         }
      plot(x, y, col = col.rc[fq$spp[, 1] + 1], xlab = xLab, ylab = yLab, ...)
      #draw convex hulls around common,rare spp
      xy <- cbind(x, y)
      Hrar <- chull(xy[which(fq$spp[, 1] == 0), ])
      Hcom <- chull(xy[which(fq$spp[, 1] == 1), ])
      polygon(xy[which(fq$spp[, 1] == 0), ][Hrar, ], border = col.rc[1],
              col = adjustcolor(col.rc[1], alpha.f = opacity))
      polygon(xy[which(fq$spp[, 1] == 1), ][Hcom, ], border = col.rc[2],
              col = adjustcolor(col.rc[2], alpha.f = opacity))
}
