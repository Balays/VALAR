




## read in GFF-compare results
all.gff.compare.outfile <- "all.merged.result_gff.compare.tsv.gz"
all.merged.result_gff.compare       <- fread(paste0(outdir, "/", all.gff.compare.outfile), na.strings = '')



all.merged.result_gff.compare.dRNA  <- fread(paste0('PRV-MDBIO-4cell_dRNA', "/", all.gff.compare.outfile), na.strings = '')



all.merged.result_gff.compare <- rbind(all.merged.result_gff.compare,
                                       all.merged.result_gff.compare.dRNA)



best.merged.result_gff.compare <- unique(all.merged.result_gff.compare[,
                                    .(seqnames, strand, start, end, exon_number, start.transcript, end.transcript, transcript_id)])

best.merged.result_gff.compare[,position := paste(start, end, sep = '-'), by=.I]

best.merged.result_gff.compare.sp <- dcast.data.table(best.merged.result_gff.compare, 
                                                      seqnames + strand + start.transcript + end.transcript + transcript_id ~ exon_number, 
                                                      value.var = c('position'))

nrow(best.merged.result_gff.compare.sp)




####

TR.gff <- data.table(as.data.frame(rtracklayer::import.gff2(config$TR.reads.gfffile)))
#TR.EX  <- fread(paste0(outdir, '/TR.EX.tsv'))
TR.gff.ori <- TR.gff
TR.gff <- data.table(as.data.frame(rtracklayer::import.gff2("PRV-MDBIO-4cell_dRNA/TR.reads.gff2")))

### compare transfrags with previous collection

#TR.gff.ori <- data.table(as.data.frame(rtracklayer::import.gff2("PRV-MDBIO-4cell/TR.reads.gff2")))
TR.gff.ori[,ori := T]
TR.gff.nov <- merge(TR.gff.ori[type == 'exon', c('seqnames', 'strand', 'start', 'end', 'exon_number', 'ori')], 
                    TR.gff[type == 'exon', c('seqnames', 'strand', 'start', 'end', 'exon_number', 'transcript_id', 'source', 'type', 'score', 'phase')], 
                    by=c('seqnames', 'strand', 'start', 'end', 'exon_number'), 
                    all.y=T)

TR.gff.nov <- TR.gff.nov[is.na(ori)]

length(unique(TR.gff.nov$transcript_id))
length(setdiff(unique(TR.gff.nov$transcript_id), unique(TR.gff.ori$transcript_id)))

length(unique(TR.gff.ori$transcript_id))
length(setdiff(unique(TR.gff.ori$transcript_id), unique(TR.gff.nov$transcript_id)))

length(unique(all.merged.result_gff.compare$transcript_id))


TR.gff <- rbind(TR.gff.nov, 
                TR.gff.ori[type == 'exon'], 
                fill = T)

TR.gff[, start.transcript := min(start), by=transcript_id]
TR.gff[, end.transcript := max(end), by=transcript_id]
TR.gff[, exon_number    := as.integer(exon_number)]

TR.gff[,position := paste(start, end, sep = '-'), by=.I]

TR.gff.sp <- dcast.data.table(TR.gff,
                              seqnames + strand + start.transcript + end.transcript + transcript_id ~ exon_number, 
                              value.var = c('position'))


####

merged.dt <-merge(best.merged.result_gff.compare, 
                  TR.gff, 
                  by=c('seqnames', 'strand', 'start', 'end', 'exon_number', 'start.transcript', 'end.transcript'),
                  all=T)
