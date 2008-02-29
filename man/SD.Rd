\name{SD}
\alias{SD}
\alias{tableSD}
\alias{affycomp.figure7}
\alias{affycomp.compfig7}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{SD Assessment Functions}
\description{
  These functions create assessments, figures, and tables for expression
  standard errors}
\usage{
affycomp.figure7(l,main="Figure 7")
affycomp.compfig7(l,method.names=as.character(1:length(l)),
                  main="Figure 7")
tableSD(l,method.names=NULL)
}
\arguments{
  \item{l}{a list of lists with the necessary components to create the
    Figure. See details.}
  \item{method.names}{a character vector with the names of the
    expression measures methodologies being compared. }
  \item{main}{title of the Figure.}
}
\details{
  This uses the dilution data. The exprsets need to have standard error
  estimates in the \code{assayDataElement(exprset,"se.exprs")}. Read the
  vignette for more details. The functions work similarly to those assessing
  expression measures.

  All these files need the result of \code{\link{assessSD}}
}
\value{Depends on the call.}
\author{Rafael A. Irizarry}
\examples{
data(rma.sd.assessment) ##this was produced with affycomp.assess
data(lw.sd.assessment) ##this one too
affycomp.compfig7(list(rma.sd.assessment,lw.sd.assessment))
affycomp.figure7(rma.sd.assessment)
}
\keyword{manip}
