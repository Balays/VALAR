

###### DO THE PLOTTING

fig.w.size    <- 4500
fig.h2w.ratio <- 0.15 #2/3

#### all samples on one plot
if (plot.all.together) {
  
  samples <- unique(metafilt[ , combine.groups])
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- (fig.width * (fig.h2w.ratio * (length(unique(samples))))) / 4
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', 'ALL', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  Fig1 <- FigY
  
} else {

  #### hpi 0.5
  samples <- na.omit(unique(metafilt[metafilt$Time ==  0.5 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '0.5h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
   
  #### hpi 1
  samples <- na.omit(unique(metafilt[metafilt$Time == 1 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '1h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 1.5
  samples <- na.omit(unique(metafilt[metafilt$Time == 1.5 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '1.5h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
 
  #### hpi 2
  samples <- na.omit(unique(metafilt[metafilt$Time == 2 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '2h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 2.5
  samples <- na.omit(unique(metafilt[metafilt$Time == 2.5 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '2.5h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 4
  samples <- na.omit(unique(metafilt[metafilt$Time == 4 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '4h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 6
  samples <- na.omit(unique(metafilt[metafilt$Time == 6 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '6h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 8
  samples <- na.omit(unique(metafilt[metafilt$Time == 8 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '8h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 12
  samples <- na.omit(unique(metafilt[metafilt$Time == 12 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '12h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 24
  samples <- na.omit(unique(metafilt[metafilt$Time == 24 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  #ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '24h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)
  
  #### hpi 48
  samples <- na.omit(unique(metafilt[metafilt$Time == 48 , combine.groups]))
  FigY    <- plotfun.all(samples)
  
  fig.width  <- (visto - visfrom) / fig.w.size
  fig.height <- fig.width * (fig.h2w.ratio * (length(unique(samples))))  # max(length(unique(samples))*3 + 3, 12)
  
  #ggsave(paste0(fig.dir, '/', virus, '_', prime, '_', '48h', '_', file_name_suffix, '.jpg'), plot = FigY, height = fig.height, width = fig.width, limitsize = F)

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
