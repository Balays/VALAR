##### START PLOTTING FROM HERE
#cov.counts    <- fread(paste0(outdir, '/cov.counts.tsv'),    na.strings = '')

## read back in
prime5.counts <- fread(paste0(outdir, '/prime5.counts.tsv'), na.strings = '')
prime3.counts <- fread(paste0(outdir, '/prime3.counts.tsv'), na.strings = '')



## filter now
correct.only <- T
if(correct.only) {
  prime5.counts <- prime5.counts[correct_tss == T,]
  prime3.counts <- prime3.counts[correct_tes == T,]
}


### normalize for viral read counts?
## based on correct adaptered-only
norm <- F
if(norm) {
  colsby <- colnames(prime5.counts)[!colnames(prime5.counts) %in% c('pos', 'start', 'end', 'prime5', 'count')]
  prime5.counts[, sum_count := sum(count), by=colsby]
  prime5.counts[, count     := round((count/ sum_count) * 100, 3)]
  prime5.counts[, .(sum_count = sum(count)), by=colsby]
  prime5.counts[, sum_count := NULL]
  
  colsby <- colnames(prime3.counts)[!colnames(prime3.counts) %in% c('pos', 'start', 'end', 'prime3', 'count')]
  prime3.counts[, sum_count := sum(count), by=colsby]
  prime3.counts[, count     := round((count/ sum_count) * 100, 3)]
  prime3.counts[, .(sum_count = sum(count)), by=colsby]
  prime3.counts[, sum_count := NULL]
  
}

## check if all the samples have the same library size
all(round(prime5.counts[,.(sum_count = (sum(count))), by=.(sample)][,sum_count], 0) == 1)


##calculate the mean? otherwise, the plot will SUM the counts based on the "group" column
calc.mean <- F
##
if(calc.mean) {
  prime3.counts <- prime3.counts[,.(count = round(mean(count), 2)), by=.(seqnames, strand, pos, correct_tes, hpi, Time, cell_line, group, endtype)]
  prime5.counts <- prime5.counts[,.(count = round(mean(count), 2)), by=.(seqnames, strand, pos, correct_tss, hpi, Time, cell_line, group, endtype)]
  
  #prime5.counts <- unique(prime5.counts[,.(seqnames, strand, pos, correct_tss, count, hpi, Time, cell_line, group, endtype)])
  #prime3.counts <- unique(prime3.counts[,.(seqnames, strand, pos, correct_tes, count, hpi, Time, cell_line, group, endtype)])
  
}

##
#prime5.counts[,prime5 := pos][,prime5 := NULL]
#prime3.counts[,prime3 := pos][,prime3 := NULL]

##### plotting settings

##
plot.all.together <- F

## grouping
combine.groups <- 'group' ## 'sample' ## OR:

## adapter settings
adapter.setting <- 'v4'

##
tolower_gene_name <- F

##### PLOTTING

#### Normalized
file_name_suffix <- 'correct_norm'
##
ylim       <- NULL #c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')


#ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)
#ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW_ori.tif'), plot = FigY, height = 20, width = 36, limitsize = F, dpi=300)
#Fig1B <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')


######## DONE  FOR NOW





stop()





#Fig2B <- Fig2
#ggsave(paste0(fig.dir, '/', 'Figure 2B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)



#### ylim 0-500
file_name_suffix <- 'Mean_correct_ylim500'
##
ylim       <- c(0, 500) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'Figure 1A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig1A <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

ggsave(paste0(fig.dir, '/', 'Figure 2A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig2A <- Fig2




#### ylim 0-50
file_name_suffix <- 'Mean_correct_norm'
##
ylim       <- c(0, 50) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 1B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)
#ggsave(paste0(fig.dir, '/', 'Figure 1B_REVIEW_ori.tif'), plot = FigY, height = 20, width = 36, limitsize = F, dpi=300)
Fig1B <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

Fig2B <- Fig2

ggsave(paste0(fig.dir, '/', 'SuppFig 2B_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)


#### ylim 0-5000
file_name_suffix <- 'Mean_correct_ylim500'
##
ylim       <- c(0, 5000) # c(-1000, 1000)

## prime5
source('plot_prime5_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime5.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 1A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig1A <- Fig1


## prime3
source('plot_prime3_area_settings.R')
fig.dir <- 'Figures' # '../EHV-1 dynamic article/Review 2024 nov'
source('plot_prime3.R')

ggsave(paste0(fig.dir, '/', 'SuppFig 2A_REVIEW.jpg'), plot = FigY, height = 20, width = 30, limitsize = F, dpi=300)

Fig2A <- Fig2


