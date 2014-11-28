bigtabulate <- function(x,
                        ccols, breaks=vector("list", length=length(ccols)),
                        table=TRUE, useNA="no",
                        summary.cols=NULL, summary.na.rm=FALSE,
                        splitcol=NULL, splitret="list") {

  if (!is.matrix(x) && !is.data.frame(x)) {
    if (class(x)!="big.matrix")
      stop("bigtabulate requires matrix, data.frame, or big.matrix objects.")
  }

  if (is.data.frame(x)) {
    for (i in 1:ncol(x)) {
      if (is.character(x[,i])) x[,i] <- factor(x[,i])
      if (is.factor(x[,i])) x[,i] <- as.integer(x[,i])
    }
    x <- as.matrix(x)
  }

  # Check and prepare ccols
  if (is.logical(ccols)) ccols <- which(ccols)
  if (length(ccols)!=length(breaks))
    stop("length(ccols) must equal length(breaks).")
  if (!is.numeric(ccols) & !is.character(ccols))
    stop("column indices must be numeric or character vectors.")
  if (is.character(ccols))
    if (is.null(colnames(x))) stop("column names do not exist.")
    else ccols <- bigmemory:::mmap(ccols, colnames(x))

  # Prepare breaks: could be a vector of length(ccols) of numbers of
  # breaks (or NA), assumed to span the ranges of the variables; or a list of
  # the same length containing a mixture of numbers of breaks (or NA) or triplets
  # of (min, max, breaks).  The result of the preparation is a matrix
  # with 3 rows and length(ccols) columns of (min, max, breaks) values (possibly NA).
  # NA indicates factor-like handling of the variable for the tabulations.
  breaks[sapply(breaks, is.null)] <- NA
  breakm <- matrix(NA, 3, length(breaks))
  if (is.numeric(breaks)) {
    if (is.matrix(x)) {
      breakm[1,!is.na(breaks)] <- apply(x[,ccols[!is.na(breaks)], drop=FALSE], 2, min, na.rm=TRUE)
      breakm[2,!is.na(breaks)] <- apply(x[,ccols[!is.na(breaks)], drop=FALSE], 2, max, na.rm=TRUE)
    } else {
      breakm[1,!is.na(breaks)] <- colmin(x, ccols[!is.na(breaks)], na.rm=TRUE)
      breakm[2,!is.na(breaks)] <- colmax(x, ccols[!is.na(breaks)], na.rm=TRUE)
    }
    breakm[3,] <- breaks
  }
  if (is.list(breaks)) {
    for (i in which(!sapply(breaks, is.na))) {
      if (length(breaks[[i]])==1) {
        if (is.matrix(x)) { 
          breakm[1,i] <- min(x[,ccols[i]], na.rm=TRUE)  
          breakm[2,i] <- max(x[,ccols[i]], na.rm=TRUE)
        } else {
          breakm[1,i] <- colmin(x, ccols[i], na.rm=TRUE)
          breakm[2,i] <- colmax(x, ccols[i], na.rm=TRUE)
        }
        breakm[3,i] <- breaks[[i]]
      } else {
        breakm[,i] <- breaks[[i]]
      }
    }
  }

  if (!is.logical(table)) stop("table must be logical.")
  table.useNA <- -1
  if (useNA=="no") table.useNA <- as.integer(0)
  if (useNA=="ifany") table.useNA <- as.integer(1)
  if (useNA=="always") table.useNA <- as.integer(2)
  if (table.useNA==-1) stop("invalid argument to useNA.")

  if (!is.logical(summary.na.rm)) stop("summary.na.rm must be logical.")
  if (is.logical(summary.cols)) summary.cols <- which(summary.cols)
  if (!is.numeric(summary.cols) && !is.character(summary.cols) &&
    !is.null(summary.cols)) {
    stop(paste("summary column indices must be numeric, logical,",
       "or character vectors."))
  }
  if (is.character(summary.cols))
    if (is.null(colnames(x))) stop("column names do not exist.")
    else summary.cols <- bigmemory:::mmap(summary.cols, colnames(x))
  if (!is.null(splitcol)) {
    if (!is.na(splitcol)) {
      if (is.logical(splitcol)) splitcol <- which(splitcol)
      if (!is.numeric(splitcol) & !is.character(splitcol))
        stop("splitcol must be numeric, logical, or character specifying one column, or NA or NULL.")
      if (is.character(splitcol))
        if (is.null(colnames(x))) stop("column names do not exist.")
        else splitcol <- bigmemory:::mmap(splitcol, colnames(x))
      if (length(splitcol)!=1) stop("splitcol must identify a single column or be NA or NULL.")
      splitcol <- as.numeric(splitcol)
    }
  }

  if (splitret!="vector" && splitret!="list" && splitret!="sparselist")
    stop("splitret must be 'vector' or 'list' or 'sparselist'")
  splitlist <- 0 # Was FALSE; this indicates vector return possibility
  if (splitret=="list") splitlist <- 1 # Was TRUE
  if (splitret=="sparselist") splitlist <- 2 # New option
  if (splitret=="vector" && is.numeric(splitcol))
    stop("splitting a specified column must return a list")
  splitlist <- as.integer(splitlist)

  # splitcol=NULL	Don't return any map type of anything.
  # splitcol=NA		Essentially split 1:nrow(x)
  # splitcol=a column   Split this single column.
  # splitlist=2: an option for the sparse list representation.
  # splitlist=1 by default: the return is a list of either split 1:nrow(x) or col entries
  # splitlist=0: only valid if splitcol==NA, in which case a vector of as.numeric(factor) entries
  summary <- is.numeric(summary.cols)
  if (is.numeric(splitcol) && splitlist==0)
    stop("vector split structure is not allowed on a column")

  if (!is.matrix(x)) {
    ans <- .Call("BigMatrixTAPPLY", x@address, as.numeric(ccols), as.numeric(breakm),
                 table, table.useNA,
                 summary, as.numeric(summary.cols), 
                 summary.na.rm, splitcol, splitlist)
  } else {
    if (is.integer(x)) {
      ans <- .Call("RIntTAPPLY", x, as.numeric(ccols), as.numeric(breakm),
                   table, table.useNA,
                   summary, as.numeric(summary.cols), 
                   summary.na.rm, splitcol, splitlist)
    } else {
      ans <- .Call("RNumericTAPPLY", x, as.numeric(ccols), as.numeric(breakm),
                   table, table.useNA,
                   summary, as.numeric(summary.cols), 
                   summary.na.rm, splitcol, splitlist)
    }
  }

  # The return will always contain
  # - ans$levels, a list of length(ccols) of factor levels possibly plus "NA"
  #
  # It will contain at least one of the following:
  # - ans$table:	vector of length prod(dim())
  # - ans$summary:	list of length prod(dim()) of cell summary matrices with 5 columns
  # - ans$split:	list of length prod(dim()) containing the split or map result;
  #                     nothing returned if is.null(splitcol), so don't return anything from C++
  #                     in that case.  Or a vector of the factor levels.

  dn <- lapply(ans$levels, function(x) { x[is.na(x)] <- "NA"; return(x) })
  ans$levels <- NULL
  if (table) ans$table <- array(ans$table, dim=sapply(dn, length), dimnames=dn)
  if (summary){
     ans$summary <- array(ans$summary, dim=sapply(dn, length), dimnames=dn)
  }
  #if (!is.null(splitcol)) {
  #  names(ans$split) <- unlist(dn) # This currently only works for 1-factor splits.
  #}

  if (length(ans)==1) return(ans[[1]])
  return(ans)

}

bigsplit <- function(x, ccols,
                     breaks=vector("list", length=length(ccols)), useNA="no", 
                     splitcol=NA, splitret="list") {

  return(bigtabulate(x, ccols=ccols, breaks=breaks,
                     table=FALSE, useNA=useNA,
                     splitcol=splitcol, splitret=splitret))
}

bigtable <- function(x, ccols,
                     breaks=vector("list", length=length(ccols)),
                     useNA="no") {

  return(bigtabulate(x, ccols=ccols, breaks=breaks,
                     table=TRUE, useNA=useNA,
                     splitcol=NULL))
}

bigtsummary <- function(x, ccols,
                        breaks=vector("list", length=length(ccols)), useNA="no",
                        cols, na.rm=FALSE) {

  return(bigtabulate(x, ccols=ccols, breaks=breaks,
                     table=FALSE, useNA=useNA,
                     summary.cols=cols, summary.na.rm=na.rm))

}

#
# April 25, 2010: we decided not to include bigaggregate() at this point.
# It just wasn't adding that much, and the performance is poor for all but
# the largest examples.
#
#bigaggregate <- function(x, stats, usesplit=NULL,
#                         ccols=NA, breaks=vector("list", length=length(ccols)), 
#                         useNA="no", distributed=FALSE, rettype="celllist", 
#                         simplify=TRUE) {
#  if (is.null(usesplit)) {
#    usesplit <- bigsplit(x, ccols=ccols, breaks=breaks, useNA=useNA, 
#      splitcol=NA, splitret="list")
#  }
#
#  # At this point I have usesplit, which is the map.  Everything else is much like I had
#  # previously in commented code, below.
#
#  if (is.data.frame(x)) {
#    for (i in 1:ncol(x)) {
#      if (is.character(x[,i])) x[,i] <- factor(x[,i])
#      if (is.factor(x[,i])) x[,i] <- as.integer(x[,i])
#    }
#    x <- as.matrix(x)
#  }
#
#  require(foreach)
#  if (is.null(getDoParName())) {
#    registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
#  }
#  if (!is.list(stats[[1]])) stats <- list(stats=stats)
#  if (is.null(names(stats))) stop("stats must be a named list")
#  if (length(unique(names(stats)))!=length(stats))
#    stop("names of stats list must be unique")
#
#  if (!is.null(getDoParName()) && getDoParName()!="doSEQ") {
#    require(bigmemory)
#    if (distributed) bf <- "" 
#    else bf <- NULL
#    if (is.matrix(x)) {
#      x <- as.big.matrix(x, backingfile=bf)
#      warning("Temporary shared big.matrix created for parallel calculations.")
#    }
#    if (!is.shared(x) && !is.filebacked(x)) {
#      x <- deepcopy(x, backingfile=bf)
#      warning("Temporary shared big.matrix created for parallel calculations.")
#    }
#  }
#
#  # Now prepare the arguments.
#  for (i in 1:length(stats)) {
#    thisname <- names(stats)[i]
#    args <- stats[[i]]
#    if (!is.list(args))
#      stop(paste("stats element", thisname, "needs to be list."))
#    if (!is.function(args[[1]]) && !is.character(args[[1]])) {
#      stop(paste("first argument of stats element", thisname, 
#        "needs to be a function."))
#    }
#    if (is.character(args[[2]])) {
#      if (is.null(colnames(x))) stop("column names do not exist.")
#      else args[[2]] <- bigmemory:::mmap(args[[2]], colnames(x))
#    }
#    if (!is.numeric(args[[2]])) args[[2]] <- as.numeric(args[[2]])
#    stats[[i]] <- args
#    names(stats)[i] <- thisname
#  }
#
#  # Here, process the chunks of data.
#  xdesc <- if (!is.matrix(x)) describe(x) else NULL
#  fans <- foreach(i=usesplit) %dopar% {
#    if (is.null(i)) {
#      temp <- as.list(rep(NA, length(stats)))
#      names(temp) <- names(stats)
#      return(temp)
#    }
#    if (!is.null(xdesc)) x <- attach.big.matrix(xdesc)
#    temp <- vector("list", length=0)
#    for (j in names(stats)) {
#      farg <- stats[[j]]
#      tempname <- names(formals(farg[[1]]))[1]
#      if (is.character(farg[[1]])) farg[[1]] <- as.symbol(farg[[1]])
#      farg[[2]] <- x[i,farg[[2]],drop=FALSE]
#      if (!is.null(tempname)) names(farg)[2] <- tempname
#      else names(farg)[2] <- ""
#      mode(farg) <- "call"
#      temp[[j]] <- eval(farg)
#    }
#    rm(farg)
#    gc()
#    
#    return(temp)
#  }
#
#  #temp <- array(temp, dim=sapply(dn, length), dimnames=dn)
#
#  if (rettype=="statlist") {
#    # Provide list of length(stats) of arrayed answers.
#    z <- NULL
#    for (j in names(stats)) {
#      res <- lapply(temp, function(x) return(x[[j]]))
#      if (all(sapply(res, length)==1) && simplify)
#        res <- array(unlist(res), dim=sapply(dn, length), dimnames=dn)
#      else {
#        # Here, there could be some empty cells with single NA values that need replication:
#
##        if (length(unique(sapply(res, length)))==2) {
##          nc <- max(unique(sapply(temp, length)), na.rm=TRUE)
##          usenames <- names(temp[[which(sapply(temp, length)==nc)[1]]])
##          for (k in which(sapply(temp, length)==1)) {
##            if (is.na(temp[[k]])) temp[[k]] <- as.numeric(rep(NA, nc))
##            names(temp[[k]]) <- usenames
##          }
##        }
#
#        res <- array(temp, dim=sapply(dn, length), dimnames=dn)
#      }
#      z[[j]] <- temp
#    }
#  }
#
#  z[is.null(z)] <- NULL
#
#  return(fans)
#
#}

