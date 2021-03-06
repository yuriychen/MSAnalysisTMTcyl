% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_normalization.R
\name{data_norm_irs}
\alias{data_norm_irs}
\title{Normalize MS TMT data with IRS.}
\usage{
data_norm_irs(prot_dat, condition_num, repeat_num)
}
\arguments{
\item{prot_dat}{A DataFrame from read_maxquant_prot() with only numeric value.}
}
\value{
A DataFrame.
}
\description{
Normalize MS TMT data with Internal Reference Signal.
}
\details{
Input dataframe from read_maxquant_prot() function, then return a data set after normalization with IRS.
}
