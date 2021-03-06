% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigtabulate.R
\name{bigtabulate}
\alias{bigtabulate}
\alias{bigsplit}
\alias{bigtable}
\alias{bigtsummary}
\title{Extended Tabular Operations for Both matrix and big.matrix Objects}
\usage{
bigtabulate(
  x,
  ccols,
  breaks = vector("list", length = length(ccols)),
  table = TRUE,
  useNA = "no",
  summary.cols = NULL,
  summary.na.rm = FALSE,
  splitcol = NULL,
  splitret = "list"
)

bigsplit(
  x,
  ccols,
  breaks = vector("list", length = length(ccols)),
  useNA = "no",
  splitcol = NA,
  splitret = "list"
)

bigtable(
  x,
  ccols,
  breaks = vector("list", length = length(ccols)),
  useNA = "no"
)

bigtsummary(
  x,
  ccols,
  breaks = vector("list", length = length(ccols)),
  useNA = "no",
  cols,
  na.rm = FALSE
)
}
\arguments{
\item{x}{a \code{\link[bigmemory]{big.matrix}} or a 
\code{\link{data.frame}} or a \code{\link{matrix}}.}

\item{ccols}{a vector of column indices or names specifying which 
columns should be used for conditioning (e.g. for building a contingency 
table or structure for tabulation).}

\item{breaks}{a vector or list of \code{length(ccols)}.  If a vector,
\code{NA} indicates that the associated column should be treated like a
factor (categorical variable), while an integer value indicates that the
range of the associated column should be broken into a specified number of
evenly-spaced bins (histogram-like).  If a list, \code{NA} triggers the
factor-like handling, a single number triggers bin-like behavior, while a
triplet (min,max,breaks) indicates that the bin-like behavior should be on
a restricted range rather than on the range of data for that column.  See
\code{\link[biganalytics]{binit}} for similar specification of this option.}

\item{table}{if \code{TRUE}, a list of table counts will be returned.}

\item{useNA}{whether to include extra '\code{NA}' levels in the table.}

\item{summary.cols}{column(s) for which table summaries will be calculated.}

\item{summary.na.rm}{if \code{TRUE}, \code{NA}s are removed from table
summary calculations.}

\item{splitcol}{if \code{NA}, the indices which correspond to
table-levels are returned.  If numeric, the corresponding column
values will be returned in a list corresponding to table-levels.  If
\code{NULL}, then there is no splitting at all.}

\item{splitret}{if \code{"list"}, the \code{splitcol} value is returned
as a list.  When \code{splitcol} is \code{NA}, \code{splitret} may
be \code{"vector"}.  Finally, \code{"sparselist"} may be a useful option
when the full-blown splitting structure has many unrepresented "cells";
this is like using the \code{drop=TRUE} option to \code{\link{split}}.}

\item{cols}{with \code{bigtsummary}, which column(s) should be conditionally
summarized?  This (or these) will be passed on as \code{summary.cols}.}

\item{na.rm}{an obvious option for summaries.}
}
\value{
array-like object(s), each similar to what is returned by
\code{\link{tapply}} and the associated \R functions.
}
\description{
This package extends the \pkg{bigmemory} package, but the 
functions may also be used with traditional \R \code{matrix} and 
\code{data.frame} objects. The function \code{\link{bigtabulate}} is 
exposed, but we expect most users will prefer the higher-level functions 
\code{bigtable}, \code{bigtsummary}, and \code{bigsplit}. Each of these
functions provides functionality based on a specified conditional 
structure.  In other words, for every cell of a (possibly multidimensional) 
contingency table, they provide (or tabulate) some useful conditional 
behavior (or statistic(s)) of interest.  At the most basic level, this 
provides an extremely fast and memory-efficient alternative to 
\code{\link{table}} for matrices and data frames.
}
\details{
This package concentrates on conditional stuctures and calculations,
much like \code{\link{table}}, \code{\link{tapply}}, and \code{\link{split}}.
The functions are juiced-up versions of the base \R functions;
they work on both regular \R matrices and data frames, but are specialized
for use with \pkg{bigmemory} and (for more advanced usage) \pkg{foreach}.
They are particularly fast and memory-efficient.  We have found that
\code{bigsplit} followed by \code{\link{lapply}} or \code{\link{sapply}}
can be particularly effective, when the subsets produced by the split
are of reasonable size.  For intensive calculations, subsequent use of
\code{foreach} can be helpful (think: parallel apply-like behavior).

When \code{x} is a \code{matrix} or a \code{data.frame}, some additional
work may be required.  For example, a character column of a \code{data.frame}
will be converted to a \code{\link{factor}} and then coerced to numeric
values (factor level numberings).

The conditional structure is specified via \code{ccols} and \code{breaks}.
This differs from the design of the base \R functions but is at the root
of the gains in speed and memory-efficiency.  The \code{breaks} may seem
distracting, as most users will simply condition on categorical-like columns.
However, it provides the flexibility to \dQuote{bin} \dQuote{continuous},
column(s) much like a histogram.  See \code{\link[biganalytics]{binit}} for
another example
of this type of option, which can be particularly valuable with massive 
data sets.

A word of caution: if a \dQuote{continuous} variable is not \dQuote{binned},
it will be treated like a factor and the resulting conditional structure will
be large (perhaps immensely so).
The function uses left-closed intervals [a,b) for the "binning" behavior,
when specified, except in the right-most bin, where the interval is entirely
closed.

Finally, \code{bigsplit} is somewhat more general than \code{split}.
The default behavior (\code{splitcol=NA})
returns a split of \code{1:nrow(x)} as a list
based on the specified conditional structure.  However, it may also
return a vector of cell (or category) numbers.  And of course it may
conduct a split of \code{x[,splitcol]}.
}
\examples{
data(iris)

# First, break up column 2 into 5 groups, and leave column 5 as a
# factor (which it is).  Note that iris is a data.frame, which is
# fine.  A matrix would also be fine.  A big.matrix would also be fine!
bigtable(iris, ccols=c(2, 5), breaks=list(5, NA))

iris[,2] <- round(iris[,2]) # So columns 2 and 5 will be factor-like
                            # for convenience in these examples, below:

ans1 <- bigtable(iris, c(2, 5))
ans1
# Same answer, but with nice factor labels from table(), because
# table() handles factors.  bigtable() uses the numeric factor
# levels only.
table(iris[,2], iris[,5])

# Here, our formulation is simpler than split's, and is faster and
# more memory-efficient:
ans2 <- bigsplit(iris, c(2, 5), splitcol=1)
ans2[1:3]
split(iris[,1], list(col2=factor(iris[,2]), col5=iris[,5]))[1:3]
}
