\name{frailty}
\alias{frailty}
\title{Extract group or spatial frailties}
\description{Summarize posterior frailties from a fitted exchangeable or spatial model.}
\usage{frailty(object, stage = c("survival", "split"), ...)}
\arguments{
  \item{object}{A fitted \code{stanSPsurv} object.}
  \item{stage}{The frailty process to summarize.}
  \item{...}{Reserved for extensions.}
}
\value{A data frame containing posterior summaries by group or area.}
