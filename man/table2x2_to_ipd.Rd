\name{table2x2_to_ipd}
\alias{table2x2_to_ipd}
\title{Expand Diagnostic Accuracy Tables to Participant-Level Rows}
\usage{table2x2_to_ipd(data)}
\arguments{
  \item{data}{A data frame with \code{study_id}, \code{TP}, \code{FN},
  \code{FP}, and \code{TN}.}
}
\value{A data frame with one row per participant.}
\description{Expands raw study-level 2-by-2 tables into the participant-level
representation used by the package's bivariate GLMM.}
\examples{
data("scoff", package = "dtaIPDMA")
ipd <- table2x2_to_ipd(scoff)
head(ipd)
}

