\name{plot_dta_effects}
\alias{plot_dta_effects}
\title{Plot Empirical Bayes Effects}
\usage{plot_dta_effects(x, type = c("density", "qq"), weighted = TRUE, ...)}
\arguments{
  \item{x}{An object returned by \code{weighted_hz_test()}.}
  \item{type}{Either \code{"density"} or \code{"qq"}.}
  \item{weighted}{Plot weighted standardized effects when true.}
  \item{...}{Additional graphical arguments.}
}
\value{The input object, invisibly.}
\description{Produces marginal density or normal Q-Q plots for empirical Bayes
effects.}

