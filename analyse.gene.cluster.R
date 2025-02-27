analyse.gene.clusters <- F
if(analyse.gene.clusters) {
  
  ## Import gene clusters
  
  gene.clusters.all <- fread('EHV-1.genes.TSS.TES.txt')
  #gene.clusters.all[,TES.canonic := fifelse(TES.canonic == cluster_TES, )]
  gene.clusters.all <- gene.clusters.all[order(CDS.start),]
  ORFs  <- unique(gene.clusters.all$gene)
  genes <- unique(gene.clusters.all$gene)
  
  ## add unknown to NA kinetic class
  gene.clusters.all[,Kinetic_class := fifelse(is.na(Kinetic_class), 'unknown', Kinetic_class)]
  
  ## add non-coding genes
  #
  add.nc.genes <- T
  if(add.nc.genes) {
    feature.nc <- as.data.frame(gene.clusters.all[type=='NC-gene',.(seqnames, start=CDS.start, end=CDS.end, strand, type, gene, part=1, ID=gene)])
    feature.nc$gene_name <- feature.nc$gene
    #feature.nc$gene_name[feature.nc$gene == 'NOIR']  <- 'NOIR1'
    #feature.nc$gene_name[feature.nc$gene == 'NOIR2'] <- 'NOIR1'
    #feature.nc$gene_name[feature.nc$gene == 'NOIR2-2'] <- 'NOIR2'
    #feature.nc$gene_name[feature.nc$gene == 'NOIR-2'] <- 'NOIR2'
    
    feature.df$gene_name <- feature.df$gene
    feature.df$gene_name <- gsub('_2', '', feature.df$gene_name)
    feature.df <- rbind(feature.df, feature.nc)
    
  } else {
    
    feature.df$gene_name <- feature.df$gene
    feature.df$gene_name <- gsub('_2', '', feature.df$gene_name)
  }
  
  feature.dt <- data.table(feature.df)
  
  
  ## Gene Regions (combination of genes and gene_clusters)
  # Function to append cluster name if different from last gene
  append_cluster <- function(genes, cluster) {
    if (tail(genes, 1) != cluster) {
      return(c(genes, cluster))
    } else {
      return(genes)
    }
  }
  
  # Apply the function by each gene_cluster
  genes_and_clusters <- gene.clusters.all[, .(CompleteList = append_cluster(gene, gene_cluster)), by = gene_cluster]
  genes_and_clusters[,gene_region := unlist(CompleteList)]
  
  genes_and_clusters <- merge(genes_and_clusters,
                              unique(gene.clusters.all[,.(gene_cluster, Kinetic_class)]),
                              by='gene_cluster')
  
  # Combine all results into one vector and ensure uniqueness
  gene_regions  <- unique(genes_and_clusters$gene_region)
  
}
