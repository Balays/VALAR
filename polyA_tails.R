
bamfiles <- list.files('C:/data/PRV_4cell/dorado_pass', '.bam', full.names = T)



params <- ScanBamParam(what = what, flag = flag, tag = 'pt')
bamfile <- bamfiles[1]

dt.from.bam2 <- function(bamfile, params) {
  bamname <- gsub('.*\\/', '', bamfile)
  bamname <- gsub(pattern, '', bamname)
  message('Start analyzing ', bamfile, '...')
  
  bam <- as.data.frame(scanBam(bamfile, param = params))
  setDT(bam)
  
  bam[,sample := ..bamname]

  return(bam)
}

#bam <- dt.from.bam2(bamfiles[1], params)

bam.list <- purrr::map(bamfiles, dt.from.bam2, params)

bam <- rbindlist(bam.list)

polya.dt <- bam[,.(sample, qname, pt)]

fwrite(polya.dt, 'dRNA_polyA_tail_lengths.tsv.gz')
