\name{exprset.log}
\alias{exprset.log}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Take log base 2 of Expression}
\description{
  Take log base 2 of the expression matrix in an \code{\link[Biobase]{exprSet}}
}
\usage{
exprset.log(exprset)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{exprset}{\code{\link[Biobase]{exprSet}}}
}
\details{
  This functions takes log base 2 of the \code{exprs} slot (expression
  matrix) in an  \code{\link[Biobase]{exprSet}}. Negatives are converted
  to the smallest non-negative entry.
}
\value{
 An \code{\link[Biobase]{exprSet}}
}
\author{Rafael A. Irizarry}
\keyword{manip}
