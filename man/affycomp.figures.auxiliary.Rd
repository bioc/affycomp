\name{affycomp.figures.auxiliary}
\alias{affycomp.figures.auxiliary}
\alias{affycomp.figure1}
\alias{affycomp.figure2}
\alias{affycomp.figure3}
\alias{affycomp.figure4a}
\alias{affycomp.figure4b}
\alias{affycomp.figure5a}
\alias{affycomp.figure5b}
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
affycomp.figure2(l,main="Figure 2")
affycomp.figure3(l, main = "Figure 3")
affycomp.figure4a(l, main = "Figure 4a")
affycomp.figure4b(l, main = "Figure 4a")
affycomp.figure5a(l, main = "Figure 5a",maxfp=100)
affycomp.figure5b(l, main = "Figure 5b",maxfp=100)
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
