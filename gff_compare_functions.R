


# GFF-COMPARE FUNCTIONS


###
safeBPParam <- function(nworkers) {
  if (.Platform$OS.type=="windows") {
    BiocParallel::SerialParam(bpparam())
  } else {
    BiocParallel::MulticoreParam(workers=nworkers, tasks = 20, progressbar = TRUE)
  }
}


###
run_gff_compare <- function(ref_mRNA, 
                            viral.mrna,
                            query_gff = 'LoRTIA_SO.TR_reads.gff2', 
                            out_dir = 'gffcompare_iterative',
                            iterative=T,
                            gff_compare_path="/mnt/c/ubuntu/programs/gffcompare-0.12.6.Linux_x86_64/gffcompare") {
  
  dir.create(out_dir)
  
  ## subset the reference gff to one transcript
  if (iterative) {
    ref_mRNA.gff <- viral.mrna[transcript_id == ref_mRNA]
    
  } else {
    ref_mRNA <- 'ALL_TRs'
    ref_mRNA.gff <- viral.mrna
    
  }
  
  temp_gff_filename <- paste0(out_dir, '/', ref_mRNA, '.gff2')
  rtracklayer::export.gff2(ref_mRNA.gff, temp_gff_filename)
  
  output_file <- paste0(out_dir, '/', ref_mRNA, '.gffcompare')
  
  
  ##
  ## run gffcompare on this subset
  command = paste(gff_compare_path,
                  "-r", temp_gff_filename,
                  "-o", output_file,
                  "-V -e 25 --strict-match",
                  query_gff)
  
  system(command)
  
  ## remove the unnecessary "*.tracking" file
  file.remove(list.files(out_dir, "*.tracking", full.names=T))
  
}



###
## function
merge_gff_compare <- function(ref_mRNA, viral.mrna,
                              #temp_gff_filename = paste0(out_dir, '/', ref_mRNA, '.gff2'),
                              #query_gff = 'LoRTIA_SO.TR_reads.gff2', 
                              out_dir = 'gffcompare_iterative', 
                              output_file = paste0(out_dir, '/', ref_mRNA, '.gffcompare'),
                              result_file = paste0(output_file, '.annotated.gtf')
) {
  
  
  #result_file        <- paste0(output_file, '.annotated.gtf')
  message('processing: ', result_file, '...')
  try({
    if (!file.exists(result_file)) {
      message('Error! ', result_file, ' does not exist!')
      return(data.table(NULL))
    } else if (file.size(result_file) == 0) { 
      message('Error! ', result_file, ' is empty!') 
      return(data.table(NULL))
      
    } else {
      merged.result_gff.compare <- data.table(NULL)
      
      try({
        result_gff.compare <- data.table(as.data.frame(rtracklayer::import.gff(result_file)))
        result_gff.compare[,prime3 := ifelse(strand == '+', end,   start)]
        result_gff.compare[,prime5 := ifelse(strand == '+', start, end)]
        
        ###
        DT.tr <- result_gff.compare[type == 'transcript']
        setnames(DT.tr, old=c('start', 'end'), new=c('start.transcript', 'end.transcript'))
        annot_cnames <- colnames(DT.tr)[-c(1:9)]
        annot_cnames <- annot_cnames[!annot_cnames %in% c('gene_id', 'exon_number', 'prime3', 'prime5')]
        cnames <- c('seqnames', 'strand', 'start.transcript', 'end.transcript', annot_cnames)
        ## omit NA cmp_ref
        DTx   <- DT.tr[!is.na(cmp_ref), ..cnames]
        
        ###
        DT.tr <- viral.mrna[type == 'transcript']
        setnames(DT.tr, old=c('start', 'end'), new=c('start.transcript', 'end.transcript'))
        ##
        #DTy  <- DT.tr[, !c("type", "phase", 'source', ), with = FALSE]
        DTy  <- DT.tr[, .(seqnames, strand, start.transcript, end.transcript, transcript_id)]
        
        ### Merge reference transcripts with read.tr-s
        DTxy <- merge(DTx, DTy,
                      by.x = c('seqnames','cmp_ref'),
                      by.y = c('seqnames','transcript_id'),
                      suffixes=c('', '.ref')
                      #, all =T
        )
        
        stopifnot(length(unique(DTxy[,transcript_id])) == nrow(DTxy[,]))
        
        ### merge READ exons and transcripts
        DT.ex <- result_gff.compare[type == 'exon', .(seqnames, strand, start, end, transcript_id)]
        ## Add prime5 and prime3 columns
        DT.ex[, prime5 := ifelse(strand == '+', start, end)]
        DT.ex[, prime3 := ifelse(strand == '+', end, start)]
        ## order to calculate exon number
        DT.ex <- DT.ex[order(seqnames, strand, transcript_id, 
                             fifelse(strand == '+', prime5, -prime5))]
        ## Assign exon numbers
        DT.ex <- DT.ex[, exon_number := 1:.N, by = .(transcript_id)]
        ## exon count per transcript
        DT.ex[, exon_count := max(exon_number), by = .(transcript_id)]
        
        
        ## merge
        DTxy  <- merge(DT.ex, DTxy, 
                       by=c('seqnames', 'strand', 'transcript_id'),
                       suffixes=c('', '.transcript'))
        
        
        
        ### Merge ref exons and transcripts
        DT.ex <- viral.mrna[type == 'exon', .(seqnames, strand, start, end, transcript_id) ]
        ## Add prime5 and prime3 columns
        DT.ex[, prime5 := ifelse(strand == '+', start, end)]
        DT.ex[, prime3 := ifelse(strand == '+', end, start)]
        
        
        ## oredr to calculate exon number
        DT.ex <- DT.ex[order(seqnames, strand, transcript_id, 
                             fifelse(strand == '+', prime5, -prime5))]
        ## Assign exon numbers
        DT.ex <- DT.ex[, exon_number := 1:.N, by = .(transcript_id)]
        ## exon count per transcript
        DT.ex[, exon_count.ref := max(exon_number), by = .(transcript_id)]
        
        ## merge
        DTxy  <- merge(DTxy, DT.ex,
                       by.x=c('seqnames', 'strand.ref', 'cmp_ref', 'exon_number'),
                       by.y=c('seqnames', 'strand', 'transcript_id', 'exon_number'),
                       suffixes=c('', '.ref'),
                       allow.cartesian = T, all.x = T)
        
        
        ## order
        DTxy <- DTxy[order(seqnames, strand, transcript_id, 
                           fifelse(strand == '+', prime5, -prime5))]
        
        
        merged.result_gff.compare <- DTxy 
        
        ### DISTANCE CALCULATION
        merged.result_gff.compare <- merged.result_gff.compare[order(transcript_id, start)]
        
        ### calculate exon-level prime3 and prime5 distances
        merged.result_gff.compare[, distance.start  := start.ref - start]
        merged.result_gff.compare[, distance.end    := end - end.ref]
        
        merged.result_gff.compare[, distance_prime5 := ifelse(strand == '+', distance.start, distance.end)]
        merged.result_gff.compare[, distance_prime3 := ifelse(strand == '+', distance.end,   distance.start)]
        
        ## regardless of strand, distance is negative if the read is shorter than the reference and positive if it is higher
        
        merged.result_gff.compare[, abs_distance    := abs(distance_prime5) + abs(distance_prime3)]
        
        ### Calculate intron start and end and sizes
        ## for query
        merged.result_gff.compare[, intron_start := c(NA, head(end, -1) + 1), by=transcript_id]
        merged.result_gff.compare[, intron_end   := start - 1]
        merged.result_gff.compare[, intron_width := intron_end - intron_start + 1]
        merged.result_gff.compare[is.na(intron_width), intron_width := 0]
        
        ## for ref
        merged.result_gff.compare[, intron_start.ref := c(NA, head(end.ref, -1) + 1), by=transcript_id]
        merged.result_gff.compare[, intron_end.ref   := start.ref - 1]
        merged.result_gff.compare[, intron_width.ref := intron_end.ref - intron_start.ref + 1]
        merged.result_gff.compare[is.na(intron_width.ref), intron_width.ref := 0]
        
        
        ### Calculate transcript-level distances
        merged.result_gff.compare[, distance.start.tr  := start.transcript.ref - start.transcript]
        merged.result_gff.compare[, distance.end.tr    := end.transcript - end.transcript.ref]
        
        merged.result_gff.compare[, distance_prime5.tr := ifelse(strand == '+', distance.start.tr, distance.end.tr)]
        merged.result_gff.compare[, distance_prime3.tr := ifelse(strand == '+', distance.end.tr,   distance.start.tr)]
        
        merged.result_gff.compare[, abs_distance.tr    := abs(distance_prime5.tr) + abs(distance_prime3.tr)]
        
        
        ### Sum up intron size abs differences for each TR (junction distance)
        merged.result_gff.compare[, abs_junc_distance.tr := sum(abs(intron_width - intron_width.ref), na.rm = T), 
                                  by=.(seqnames, strand, transcript_id, cmp_ref, class_code)]
        ###
        
        
        ### write comparison output
        fwrite(merged.result_gff.compare, paste0(output_file, '.result.tsv'), sep='\t')
        
      })
      return(merged.result_gff.compare)
    }
  })
}



