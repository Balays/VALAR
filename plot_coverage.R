

###### DO THE PLOTTING

fig.w.size    <- 5000
fig.h2w.ratio <- 0.30 #2/3

#### all samples on one plot
if (plot.all.together) {
  
  samples <- unique(metafilt[ , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- (fig.width * (fig.h2w.ratio * (length(unique(samples))))) / 4
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', 'ALL', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  Fig3 <- FigY
  
} else {

  #### hpi 1
  samples <- unique(metafilt[metafilt$Time == 1 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '1h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 2
  samples <- unique(metafilt[metafilt$Time == 2 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '2h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 4
  samples <- unique(metafilt[metafilt$Time == 4 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '4h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 6
  samples <- unique(metafilt[metafilt$Time == 6 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '6h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 8
  samples <- unique(metafilt[metafilt$Time == 8 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '8h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 12
  samples <- unique(metafilt[metafilt$Time == 12 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '12h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 24
  samples <- unique(metafilt[metafilt$Time == 24 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '24h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 48
  samples <- unique(metafilt[metafilt$Time == 48 , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '48h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)

}

### Ezeket kell allitgatni:
#visfrom; visto; fromto
#'bin_width.', bin_width
#ylim
#breakseq
#y.log10
#y.log2
#plot.title

#gene.sizes
