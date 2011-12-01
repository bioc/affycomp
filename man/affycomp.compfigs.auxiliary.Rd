\name{affycomp.compfigs.auxiliary}
\alias{affycomp.compfigs.auxiliary}
\alias{affycomp.compfig2}
\alias{affycomp.compfig2b}
\alias{affycomp.compfig3}
\alias{affycomp.compfig4a}
\alias{affycomp.compfig4b}
\alias{affycomp.compfig4c}
\alias{affycomp.compfig5a}
\alias{affycomp.compfig5b}
\alias{affycomp.compfig5cde}
\alias{affycomp.compfig5c}
\alias{affycomp.compfig5d}
\alias{affycomp.compfig5e}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Auxiliary functions to create comparitive Figures}
\description{
These functions are auxiliary function to
\code{\link{affycompPlot}}. These Figures are used to compare
expression measures. They take lists with components created by the
\code{\link{assessDilution}} and \code{\link{assessSpikeIn}}
functions.}
\usage{
affycomp.compfig2(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 2")

affycomp.compfig3(l, method.names = as.character(1:length(l)),
                  main = "Figure 3")

affycomp.compfig4a(l, method.names = as.character(1:length(l)),
                   add.legend = TRUE, main = "Figure 4a")

affycomp.compfig4b(l, method.names = as.character(1:length(l)),
                   add.legend = TRUE, main = "Figure 4b")

affycomp.compfig4c(l, method.names = as.character(1:length(l)),
                   add.legend = TRUE, rotate=TRUE, main = "Figure 4c")

affycomp.compfig5a(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 5a", maxfp=100)

affycomp.compfig5b(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 5b", maxfp=100)

affycomp.compfig5cde(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 5c", maxfp=100,
                  type=c("low","med","high"))

affycomp.compfig5c(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 5c", maxfp=100)

affycomp.compfig5d(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 5d", maxfp=100)

affycomp.compfig5e(l, method.names = as.character(1:length(l)),
                  add.legend = TRUE, main = "Figure 5e", maxfp=100)

}
\arguments{
  \item{l}{a list of lists with the necessary components to create the
    Figure. See details.}
  \item{method.names}{a character vector with the names of the
    expression measures methodologies being compared. }
  \item{add.legend}{logical. If TRUE a legend is added.}
  \item{main}{title of the Figure.}
  \item{rotate}{in the case of compfig4c one can eiher show the actual
    local slopes or the bias (local slope minus 1).}
  \item{maxfp}{range of the false positives in ROC will be from 0 to \code{maxfp}}
  \item{type}{compfig5cdef is the engine for 5c, 5d, and
5e. \code{type} tells is which of these 4 to run.}
}
\details{
  These are similar to the functions defined in
  \code{\link{affycomp.figures.auxiliary}}. Main difference is that
  here you send lists with the result of the assessment functions as
  components. 
}
\value{
Figures are produced.
}
\author{Rafael A. Irizarry}
\examples{
library(affycompData)
data(rma.assessment)
data(mas5.assessment)
affycomp.compfig2(list(rma.assessment$Dilution,mas5.assessment$Dilution))
affycomp.compfig3(list(rma.assessment$Dilution,mas5.assessment$Dilution))
affycomp.compfig4a(list(rma.assessment$Signal,mas5.assessment$Signal))
affycomp.compfig4b(list(rma.assessment$Dilution,mas5.assessment$Dilution))
affycomp.compfig5a(list(rma.assessment$FC,mas5.assessment$FC))
affycomp.compfig5b(list(rma.assessment$FC2,mas5.assessment$FC2))
}
\keyword{hplot}% at least one, from doc/KEYWORDS










