\name{exprset.log}
\alias{exprset.log}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Take log base 2 of Expression}
\description{
  Take log base 2 of the expression matrix in an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} 
}
\usage{
exprset.log(exprset)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{exprset}{\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}}
}
\details{
  This functions takes log base 2 of the expression matrix in an
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. Negatives
  are converted to the smallest non-negative entry.
}
\value{
 An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}
}
\author{Rafael A. Irizarry}
\keyword{manip}
