


#### COVERAGE ANALYSIS ####
message('Analysing coverages ...')
source('coverage_virus.R')


### coverage summary plots
if(config$make.plots) {
  source('cov_stat_plots.R')
}
#### ####
##

#### write outputs ####
if (config$write.all) {
  
  fwrite(mapped.cov, paste0(outdir, '/mapped.cov.tsv'), sep = '\t')
  fwrite(merged_cov, paste0(outdir, '/merged_cov.tsv'), sep = '\t')
  
}

#### ####
##

