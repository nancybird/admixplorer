# admixplorer

An R package for admixplorer, a new admixture dating approach to jointly analyse modern and ancient individuals sampled at different times.

## Installation

`if (!requireNamespace("devtools", quietly = TRUE)) {

install.packages("devtools")

}

devtools::install_github("nancybird/admixplorer")`

## Usage

admixplorer(

  infile,

  outfile,
 
  method = "GLOBETROTTER",
 
  ks = "1,2,3,4",
 
  mcmc_chains = 3,
 
  sample_age_est = TRUE,
 
  plot = TRUE
)

# Arguments

infile:	Path to input file (file with 5 columns: ind, sampling_age_lower, sampling_age_upper, inferred_admixture_date, inferred_admixture_sd). no header

outfile:	Output file prefix

method:	Method to use ("GLOBETROTTER" or "DATES", or other method (will use GLOBETROTTER thresholds ))

ks:	Comma-separated string of number of clusters to try (default "1,2,3,4", thresholds only go up to 4 clusters, after that is predicts 4+)

mcmc_chains:	Number of MCMC chains per k (default 3, use only 1 for faster running)

sample_age_est:	Whether to estimate sample ages for ancient individuals (default TRUE, if FALSE will use the midpoint of the range)

plot:	Whether to create plots (default TRUE, creates plots of coindence matrix and cluster means to inform on cluster numbers)
