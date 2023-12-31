% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_sps_data.R
\name{clean_for_sps}
\alias{clean_for_sps}
\title{Clean data for function sps()}
\usage{
clean_for_sps(X, scale = TRUE)
}
\arguments{
\item{X}{Site-level variables for the target population of sites. Row names should be names of sites. X cannot contain missing data.}

\item{scale}{TRUE or FALSE. Whether to standardize numeric variables to make each variable mean zero and standard deviation one.}
}
\value{
\code{clean_for_sps} returns an object of \code{matrix} class, which we supply to \code{sps()}.
 \itemize{
   \item \code{X}: Site-level variables for the target population of sites.
 }
}
\description{
Clean data for function sps()
}
\references{
Egami and Lee. (2023+). Designing Multi-Context Studies for External Validity: Site Selection via Synthetic Purposive Sampling. Available at \url{https://naokiegami.com/paper/sps.pdf}.
}
