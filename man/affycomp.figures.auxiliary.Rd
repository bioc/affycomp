\name{affycomp.figures.auxiliary}
\alias{affycomp.figures.auxiliary}
\alias{affycomp.figure1}
\alias{affycomp.figure1b}
\alias{affycomp.figure2}
\alias{affycomp.figure2b}
\alias{affycomp.figure3}
\alias{affycomp.figure4a}
\alias{affycomp.figure4b}
\alias{affycomp.figure4c}
\alias{affycomp.figure5a}
\alias{affycomp.figure5b}
\alias{affycomp.figure5c}
\alias{affycomp.figure5d}
\alias{affycomp.figure5e}
\alias{affycomp.figure6a}
\alias{affycomp.figure6b}
\title{Auxiliary functions to create Figures}
\description{
These functions are auxiliary function to
\code{\link{affycompPlot}}. These Figures are used to assess an
expression measure. They take components created by the
\code{\link{assessDilution}} and \code{\link{assessSpikeIn}}
functions. }
\usage{
affycomp.figure1(l,main="Figure 1",xlim=NULL,ylim=NULL)
affycomp.figure1b(l,main="Figure 1b",xlim=NULL,ylim=NULL,cex=0.85,all=FALSE)
affycomp.figure2(l,main="Figure 2")
affycomp.figure2b(l,main="Figure 2b")
affycomp.figure3(l, main = "Figure 3")
affycomp.figure4a(l, main = "Figure 4a",equal.lims=FALSE)
affycomp.figure4b(l, main = "Figure 4b")
affycomp.figure4c(l, rotate=TRUE, main = "Figure 4c")
affycomp.figure5a(l, main = "Figure 5a",maxfp=100)
affycomp.figure5b(l, main = "Figure 5b",maxfp=100)
affycomp.figure5c(l, main = "Figure 5c",maxfp=100)
affycomp.figure5d(l, main = "Figure 5d",maxfp=100)
affycomp.figure5e(l, main = "Figure 5e",maxfp=100)
affycomp.figure6a(l, main = "Figure 6a",xlim = NULL, ylim = NULL)
affycomp.figure6b(l, main = "Figure 6b",xlim = NULL, ylim = NULL)
}
\arguments{
  \item{l}{A list with the necessary components to create the
    Figure. See details.}
  \item{main}{Title for the Figure.}
  \item{maxfp}{range of the false positives in ROC will be from 0 to \code{maxfp}}
  \item{xlim}{x-axis limits.}
  \item{ylim}{y-axis limits.}
  \item{cex}{size of numbers in figure 1b.}
  \item{all}{logical. If \code{TRUE} all spikeins are shown. Otherwise,
    only those resulting in smaller, realistic, fold changes are shown.}
  \item{equal.lims}{logical. If \code{TRUE} the limits of x-axis and
    y-axis will have same range.}
  \item{rotate}{in the case of compfig4c one can eiher show the actual
    local slopes or the bias (local slope minus 1).}
}
\details{
Read the vignette for more details on what each Figure is. You can read
\code{\link{assessSpikeIn}} and \code{\link{assessDilution}} to see
which assessments are needed.
}
\value{
Figures are produced.
}
\author{Rafael A. Irizarry}
\examples{
%%some examples
%%the functions not shown in example need another assessment function
data(rma.assessment)
affycomp.figure1(rma.assessment$MA)
affycomp.figure2(rma.assessment$Dilution)
affycomp.figure3(rma.assessment$Dilution)
affycomp.figure4a(rma.assessment$Signal)
affycomp.figure4b(rma.assessment$Dilution)
affycomp.figure5a(rma.assessment$FC)
affycomp.figure5b(rma.assessment$FC2)
affycomp.figure6a(rma.assessment$FC)
affycomp.figure6b(rma.assessment$FC)
}
\keyword{hplot}
