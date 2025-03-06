##### Import libraries and functions
library(cowplot)
library(grid)
library(ggplotify)
library(hrbrthemes, quietly = T)
library(gggenes)
library(tidygenomics)
library(ggh4x)

source('_WF.part0.R')




#### Settings
by <- c("seqnames", "start",  "end")


## colors and sizes
palette <-  c(pal_npg()(10), pal_aaas()(10))[c(1:10,14,18,15,13:11,16,17,20)]
alpha <- data.frame(cov_geom        = 0.6,
                    gene_geom       = 0.75,
                    unstranded_geom = 0.75)


#
visfrom <- 1
visto   <- l_genome
vis.add <- 0 # c(250,250)

sizes <- 14

gene.sizes <- data.frame(gene_arrowhead_height=8,
                         gene_arrow_body_height=6,
                         gene_label_height=4,   ## size of gene label text
                         gene_feature_height=5, ## lineheight of feature line
                         gene_feature_width=5,
                         gene_feature_label_height=1.5, ## position of gene label text
                         gene_feature_label_text=-2, ## inactive
                         unstranded_rect_height=0.5,
                         unstranded_feature_height=6,
                         unstranded_feature_width=6,
                         unstranded_feature_label_height=5,
                         unstranded_feature_label_text=6,
                         genome_feature_arrowhead_height=2.5,
                         genome_feature_arrow_body_height=2.5,
                         ## introns
                         linetype='dashed',
                         linewidth=0.01
)

angle   <- 90
## scale of annotation to coverage
genomplot.scale <- 5

breakseq <- 5000

## plot size
width <- 35; height <- 14


#### Gene annotation
genome.plotdata <- make.genome.plot.data(feature.df, feature.col = 8, feature.name = 'ID')
genome.plotdata$strand <- factor(genome.plotdata$strand, levels = c('+', '-', '*'))
setDT(genome.plotdata)
if (tolower_gene_name) { genome.plotdata[, gene_name := fifelse(grepl('ORF', gene), tolower(gene), gene)] }
genome.plotdata[, gene_name := gsub('_.*', '', gene_name)]
gene.plotdata <- as.data.frame(genome.plotdata)


##
make.fromto <- function(start=1, l=l_genome, n=4) {

  fromto    <- data.frame(from=round(seq(start, l, l/n),0))
  fromto$to <- c(fromto$from[-1], l)
  return(fromto)
}

fromto <- make.fromto(n=4)
##
add.unstranded  <- T ### Correct!!
##

##
cagefr.clust <- NA
add.cageTSS  <- F
##

#### coverage data
prime        <- 'prime3'
plot.data    <- data.table(data.frame(prime3.counts))
setnames(plot.data, old=c('pos', 'N'), new=c('prime3', 'count'), skip_absent=TRUE)
plot.data[,start := prime3][,end := prime3][,position := prime3][order(seqnames, strand)]

plot.data <- as.data.frame(plot.data)

if(prime == 'prime5') {correct_prime <- 'correct_tss'} else if (prime == 'prime3') {correct_prime <- 'correct_tes'} else {stop()}

plot.data$strand <- factor(plot.data$strand, levels = c('+', '-', '*'))


#### filter
#plot.data <- plot.data[!grepl('0hpi', plot.data$sample), ]


#### output directory
fig.dir <- paste0(outdir, '/', prime, '_plots'); dir.create(fig.dir)



#### Adapters

adapter.setting <- adapter.setting

### Dealing with TRUE and FALSE adapter-ed reads

## This combination will omit FALSE adaptered reads in those sites where there are no TRUE adaptered reads
if (adapter.setting == 'v1' ){
  keep.FALSE.at.TRUE   <- T
  crop.TRUE.from.FALSE <- F
  plot.all.ends        <- T
  crop.FALSE           <- T
  message('adapter.setting: ', adapter.setting )
} else

  ## This combination will keep FALSE adaptered reads  and plot them on another facet
  if (adapter.setting == 'v2' ){
    keep.FALSE.at.TRUE   <- F
    crop.TRUE.from.FALSE <- F
    plot.all.ends        <- F
    crop.FALSE           <- F
    message('adapter.setting: ', adapter.setting )
  } else

    ## This combination will keep FALSE adaptered reads but if there was a TRUE site it will plot them on another facet
    if (adapter.setting == 'v3' ){
      keep.FALSE.at.TRUE   <- F
      crop.TRUE.from.FALSE <- T
      plot.all.ends        <- F
      crop.FALSE           <- F
      message('adapter.setting: ', adapter.setting )
    } else

      ## This combination will keep FALSE adaptered reads
      if (adapter.setting == 'v4' ){
        keep.FALSE.at.TRUE   <- F
        crop.TRUE.from.FALSE <- F
        plot.all.ends        <- T
        crop.FALSE           <- T
        message('adapter.setting: ', adapter.setting )
      }

##
if (keep.FALSE.at.TRUE) {
  true_sites  <- unique(plot.data[plot.data[,correct_prime] == T, prime])
  plot.data   <- plot.data[is.element(plot.data[,prime], true_sites), ]
} else if (crop.TRUE.from.FALSE) {
  false_sites <- unique(plot.data[plot.data[,correct_prime] == F, prime])
  true_sites  <- unique(plot.data[plot.data[,correct_prime] == T, prime])
  false_sites <- setdiff(false_sites, true_sites)
  plot.data[is.element(plot.data[,prime], true_sites),  correct_prime] <- T
  plot.data[is.element(plot.data[,prime], false_sites), correct_prime] <- F
}
##
if (plot.all.ends) {
  plot.data$correct_tes <- T
  plot.data$correct_tss <- T
}
###



#### FigX Whole genome prime5 plots  ####
## Groups in the facets
combine.groups <- combine.groups # 'group'
if (combine.groups != F) {
  plot.data$hpi <- plot.data[,combine.groups]
} else {
  plot.data$hpi <- genome
}
##
#plot.data <-
#plot.data <- plot.data[,cnames]

## sum coverages!! IMPORTANT !!!
sum.counts <- T
if (sum.counts) {
  plot.sum <- plot.data %>% group_by(across(any_of(c('strand', 'hpi', correct_prime, prime)))) %>%
    summarise(count=sum(count)) %>%
    #summarise(count=n()) %>%
    spread(hpi, count, fill = 0) %>% gather(hpi, count, -c(1:3)) %>%
    spread(strand, count, fill=0) %>% gather(strand, count, -c(1:3))

  ### start, end, prime3,correct_tes,
  plot.sum[,'start']  <- plot.sum[,prime]
  plot.sum[,'end']    <- plot.sum[,prime]

} else {
  plot.sum <- plot.data
}

try({
  cagefr.clust <- cluster.prime3[cluster != 0,
                               .(seqnames, strand, start, end, gene=as.character(cluster), orientation = fifelse(strand == '+', 1, 0))]
})

visfrom  <- 1 # Ori.virus$visfrom[Ori.virus$ID == Ori]
visto    <- l_genome # 25000 # Ori.virus$visto[Ori.virus$ID   == Ori]

bin_width  <- 50
breakseq   <- 5000
#ylim       <- NULL #c(0, 50) # c(-1000, 1000)
ylims.gene <- c(-2, 3.5) # NULL # c(10, -20) ##
scales     <- 'free_y'

genomplot.scale <- 7

plotfun   <- function(samples=NA) {

  FigY.p    <- plot.genome.region(visfrom=visfrom,	visto=visto, samples = samples,
                                  facet_cropF = facet_nested(rows=vars(hpi), scales='fixed', drop=T),
                                  plot.data=plot.sum[plot.sum$strand == '+',],
                                  #[TIS.df.sum.all$tss == tss,],
                                  add.genome.plot=F, genome.only = F,
                                  gene.plotdata=gene.plotdata, genome=genome, prime = prime,
                                  sum.counts.in.window=T, bin_width = bin_width, add.all.pos=T, y.thresh=NA,
                                  gene.label=T,
                                  add.unstranded=T,
                                  add.feature=F,
                                  force.all.gene.down=T, angle = angle,
                                  add.cageTSS=F, cagefr.clust=NA,
                                  crop.FALSE = crop.FALSE, scales = 'free_y',
                                  geom    = geom_area(aes(x=prime3, y = count, fill=strand), position = 'identity', alpha=alpha$cov_geom, color = "black"),
                                  #geom    = geom_point(aes(x=prime5, y = count, color=strand), alpha=alpha$cov_geom, size=0.1),
                                  palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                  gene.sizes = gene.sizes, sizes = sizes,
                                  y.multip     = 2, y.log10 = F, ylim={if(is.null(ylim)) {ylim} else {c(0, ylim[2])}}, #ybreaks = -5:5,
                                  ylims.gene = ylims.gene,
                                  force.gene.y = T,
                                  gene.label.col='black',
                                  legend.position='top', plot.title = NA, #letters[2], #paste0(virus, '_', Ori),
                                  breakseq=breakseq,
                                  margins=unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                  gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
                                  #legend.position.prime='top',
                                  labels=labs(fill = "Strand"),
                                  vline=NULL,
                                  return.plot.data=F)

  FigY.ann    <- plot.genome.region(visfrom=visfrom,	visto=visto, samples = samples,
                                    facet_cropF = facet_nested(rows=vars(hpi), scales='fixed', drop=T),
                                    plot.data=plot.sum[plot.sum$strand == '+',],
                                    #[TIS.df.sum.all$tss == tss,],
                                    add.genome.plot=T, genome.only = T,
                                    gene.plotdata=gene.plotdata, genome=genome, prime = prime,
                                    sum.counts.in.window=T, bin_width = bin_width, add.all.pos=T, y.thresh=NA,
                                    gene.label=T,
                                    add.unstranded=T,
                                    add.feature=F,
                                    force.all.gene.down=T, angle = angle,
                                    add.cageTSS=add.cageTSS, cagefr.clust=cagefr.clust, cagefr.clust.dist=3,
                                    #cluster.aes = aes(xmin = region_start, xmax = region_end, y = y, forward = orientation, fill = strand, xsubmin = start, xsubmax = end, color=strand),
                                    crop.FALSE = '', scales = 'free_y',
                                    geom    = geom_area(aes(x=prime3, y = count, fill=strand), position = 'identity', alpha=alpha$cov_geom, color = "black"),
                                    #geom    = geom_point(aes(x=prime5, y = count, color=strand), alpha=alpha$cov_geom, size=0.1),
                                    palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                    gene.sizes = gene.sizes, sizes = sizes,
                                    y.multip     = 5, y.log10 = T, ylim=ylim, #ybreaks = -5:5,
                                    ylims.gene = ylims.gene, #
                                    force.gene.y = T,
                                    gene.label.col='black',
                                    legend.position='top', plot.title = NA, #letters[2], #paste0(virus, '_', Ori),
                                    breakseq=breakseq,
                                    margins=unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                    gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
                                    #legend.position.prime='top',
                                    labels=labs(fill = "Strand"),
                                    vline=NULL,
                                    return.plot.data=F)

  FigY.m    <- plot.genome.region(visfrom=visfrom,	visto=visto, samples = samples,
                                  facet_cropF = facet_nested(rows=vars(hpi), scales='fixed', drop=T),
                                  plot.data=plot.sum[plot.sum$strand == '-',],
                                  #[TIS.df.sum.all$tss == tss,],
                                  add.genome.plot=F, genome.only = F,
                                  gene.plotdata=gene.plotdata, genome=genome, prime = prime,
                                  sum.counts.in.window=T, bin_width = bin_width, add.all.pos=T, y.thresh=NA,
                                  gene.label=T,
                                  add.unstranded=T,
                                  add.feature=F,
                                  force.all.gene.down=T, angle = angle,
                                  add.cageTSS=F, cagefr.clust=NA,
                                  crop.FALSE = crop.FALSE, scales = 'free_y',
                                  geom    = geom_area(aes(x=prime3, y = count, fill=strand), position = 'identity', alpha=alpha$cov_geom, color = "black"),
                                  #geom    = geom_point(aes(x=prime5, y = count, color=strand), alpha=alpha$cov_geom, size=0.1),
                                  palette=palette[-1], alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                  gene.sizes = gene.sizes, sizes = sizes,
                                  y.multip     = 2, y.log10 = F, ylim={if(is.null(ylim)) {ylim} else {c(ylim[1], 0)}}, #ybreaks = -5:5,
                                  ylims.gene = ylims.gene,
                                  force.gene.y = T,
                                  gene.label.col='black',
                                  legend.position='top', plot.title = NA, #letters[2], #paste0(virus, '_', Ori),
                                  breakseq=breakseq,
                                  margins=unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                  gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
                                  #legend.position.prime='top',
                                  labels=labs(fill = "Strand"),
                                  vline=NULL,
                                  return.plot.data=F)

  ## overwrite x-limits
  xlim <- FigY.ann[["coordinates"]][["limits"]]$x


  FigY <- cowplot::plot_grid(FigY.p + xlim(xlim),
                             FigY.ann + xlim(xlim),
                             FigY.m + xlim(xlim),
                             rel_heights = c(3,1.25,3), align = 'v', axis=c('tlrb'), ncol = 1)
  #FigY <- FigY.ann
}


plotfun.all   <- function(samples=NA) {

  FigY    <- plot.genome.region(visfrom=visfrom,	visto=visto, samples = samples,
                                facet_cropF = facet_nested(rows=vars(hpi), scales=scales, drop=T),
                                plot.data=plot.sum, #[plot.sum$strand == '+',],
                                #[TIS.df.sum.all$tss == tss,],
                                add.genome.plot=T, genome.only = F,
                                gene.plotdata=gene.plotdata, genome=genome, prime = prime,
                                sum.counts.in.window=T, bin_width = bin_width, add.all.pos=T, y.thresh=NA,
                                gene.label=T, gene_name_col = 'gene_name',
                                add.unstranded=T,
                                add.feature=T,
                                force.all.gene.up.and.down = T,
                                force.all.gene.down=F, angle = angle,
                                add.cageTSS=add.cageTSS, cagefr.clust=F, cagefr.clust.dist=3,
                                crop.FALSE = crop.FALSE, scales = scales,
                                geom    = geom_area(aes(x=prime3, y = count, fill=strand, color = strand), position = 'identity', alpha=alpha$cov_geom),
                                #geom    = geom_point(aes(x=prime5, y = count, color=strand), alpha=alpha$cov_geom, size=0.1),
                                palette=palette, alpha=alpha, genomplot.scale=genomplot.scale,
                                gene.sizes = gene.sizes, sizes = sizes,
                                y.multip   = 3, y.log10 = F, ylim=ylim, #ybreaks = -5:5,
                                ylims.gene = NULL,
                                force.gene.y = T,
                                gene.label.col='black',
                                legend.position='top', plot.title = NA, #letters[2], #paste0(virus, '_', Ori),
                                breakseq=breakseq,
                                margins=unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                                gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene),
                                #legend.position.prime='top',
                                #labels=labs(fill = "Strand"),
                                vline=NULL,
                                xlab_name='Genomic Position', ylab_name='Read Count',
                                return.plot.data=F)


}

