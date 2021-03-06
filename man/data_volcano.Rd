% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_simple_analysis.R
\name{data_volcano}
\alias{data_volcano}
\title{Draw volcano plot}
\usage{
data_volcano(
  prot_dat_1,
  prot_dat_2,
  meta,
  df,
  xleft = -5,
  xright = 5,
  ydown = 0,
  ytop = 5,
  confidence = 0.95,
  s0 = 0.1,
  smooth_curve = TRUE
)
}
\arguments{
\item{prot_dat_1}{A DataFrame of condition 1.}

\item{prot_dat_2}{A DataFrame of condition 2.}

\item{meta}{A DataFrame of meta file.}

\item{df}{A integer of freedom dgree.}

\item{xleft}{A numeric for x-axis limit on left.}

\item{xright}{A numeric for x-axis limit on right.}

\item{ydown}{A numeric for x-axis limit on bottom.}

\item{ytop}{A numeric for x-axis limit on top.}

\item{confidence}{A numeric of confidence interval.}

\item{s0}{A numeric of fudge factor s0.}

\item{smooth_curve}{A boolean of significant curve or not.}
}
\value{
A list of dataframe and ggplot object.
}
\description{
Draw volcano plot.
}
\details{
Input dataframe from previous processing, then return a ggplot object of volcano plot.
}
