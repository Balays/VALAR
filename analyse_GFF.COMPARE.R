

## import
#all.merged.result_gff.compare <- fread(paste0(outdir, "/all.merged.result_gff.compare.tsv"), sep='\t')
###

## which TRs have an "equal" hit in the refernece?
all.merged.result_gff.compare[, equal := ifelse(class_code == '=', T, F)]

equal.ref.TR <- unique(all.merged.result_gff.compare[class_code == '=', .(seqnames, strand, transcript_id, class_code, cmp_ref)])


#### Reclassify class codes

### Reclassify '='
## Set new columns based on 5-prime and 3-prime distance thresholds
thresh.eq.prime5 <- thresh.eq.prime5 ## 10
thresh.eq.prime3 <- thresh.eq.prime3 ## 10 


## Set new columns based on 5-prime and 3-prime distance thresholds
all.merged.result_gff.compare[, equal.5_prime.tr := ifelse(abs(distance_prime5.tr) <= thresh.eq.prime5, T, F)]
all.merged.result_gff.compare[, equal.3_prime.tr := ifelse(abs(distance_prime3.tr) <= thresh.eq.prime3, T, F)]

## Overwrite '~' and 'c' to '=' where it meets the threshold
# because by default, only the EXACT matches are '='
all.merged.result_gff.compare[class_code == '~' & equal.5_prime.tr == T & equal.3_prime.tr == T, class_code := '=']  ## 
all.merged.result_gff.compare[class_code == 'c' & equal.5_prime.tr == T & equal.3_prime.tr == T, class_code := '=']

#sus <- all.merged.result_gff.compare[equal.5_prime.tr == T & equal.3_prime.tr == T & strand == strand.ref & exon_count == exon_count.ref & !class_code %in% c('=', '~', 'c'), ]
#stopifnot(nrow(sus) == 0)
## OK!

### introduce 5`- and 3`- short and long categories (number of exons and strand have to match!)
all.merged.result_gff.compare[distance_prime5.tr < -thresh.eq.prime5 & equal.3_prime.tr == T & strand == strand.ref & exon_count == exon_count.ref, class_code := "5`-short" ]
all.merged.result_gff.compare[distance_prime5.tr >  thresh.eq.prime5 & equal.3_prime.tr == T & strand == strand.ref & exon_count == exon_count.ref, class_code := "5`-long"  ]

all.merged.result_gff.compare[distance_prime3.tr < -thresh.eq.prime3 & equal.5_prime.tr == T & strand == strand.ref & exon_count == exon_count.ref, class_code := "3`-short" ]
all.merged.result_gff.compare[distance_prime3.tr >  thresh.eq.prime3 & equal.5_prime.tr == T & strand == strand.ref & exon_count == exon_count.ref, class_code := "3`-long"  ]
##

### introduce splice-variant category for those that share their 3-primes and 5-primes but have different junctions
all.merged.result_gff.compare[equal.5_prime.tr == T & equal.3_prime.tr == T & strand == strand.ref & abs_junc_distance.tr > thresh.eq.junc, class_code := "Splice-variant" ]


#### statistics for all hits
best.ref.TR <- unique(all.merged.result_gff.compare[
  best_prime3 == T & best_prime5 == T & best_junc == T, 
  .(seqnames, transcript_id, cmp_ref, min_dist_prime5, min_dist_prime3, min_dist_junc)])


## all class code freq per ref transcript
gff.compare.ref.TR.class.freq <- all.merged.result_gff.compare[
  , .N, by=.(seqnames, cmp_ref, class_code)
]

gff.compare.ref.TR.class.freq.sp <- dcast(gff.compare.ref.TR.class.freq,
                                          seqnames + cmp_ref ~ class_code, value.var = 'N'
)


## all class code freq per transcript
gff.compare.TR.class.freq <- all.merged.result_gff.compare[
  , .N, by=.(seqnames, transcript_id, class_code)
]

gff.compare.TR.class.freq.sp <- dcast(gff.compare.TR.class.freq,
                                      seqnames + transcript_id ~ class_code, value.var = 'N'
)

## find the lowest and highest distance for each cmp_ref per class_code
class_code.ref.dist.stat <- all.merged.result_gff.compare[,
                                                          .(min_dist_class = min(abs_distance.tr), max_dist_class = max(abs_distance.tr)), # Selecting the lowest distance TR-read
                                                          by = .(cmp_ref, class_code, seqnames)
]



#### keep only the "best" ref transcript per read.TR -->> USE THIS !!
best.merged.result_gff.compare <- all.merged.result_gff.compare[
  order(abs_junc_distance.tr, abs_distance.tr, abs(distance_prime5.tr), abs(distance_prime3.tr)),
  .(cmp_ref, class_code, seqnames, strand, strand.ref, start, end, abs_distance, abs_junc_distance.tr, distance_prime5.tr, distance_prime3.tr, equal.5_prime.tr, equal.3_prime.tr), # Selecting the columns of interest
  by = transcript_id
][, .SD[1], by = transcript_id]
stopifnot(nrow(best.merged.result_gff.compare) == length(unique(best.merged.result_gff.compare$transcript_id)))
best.merged.result_gff.compare[,best_result:=T]

### check if the equal results are the 'best' results in all cases -->> DOESENT MAKE SENSE !!!
equal.ref.TR <- merge(equal.ref.TR, 
                      best.merged.result_gff.compare[,.(transcript_id, cmp_ref, best_result)], 
                      by=c('transcript_id', 'cmp_ref'), all.x=T)
tr_to_check  <- equal.ref.TR[class_code == '=' & is.na(best_result), transcript_id]

tr_to_check <- all.merged.result_gff.compare[transcript_id %in% tr_to_check, ]
#stopifnot(nrow(tr_to_check) == 0)



#### make exon-level data from best transcript results
best.merged.result_gff.compare.EX <- unique(all.merged.result_gff.compare
                                            [,.(seqnames, start, end,  strand, transcript_id, prime3, prime5# width, source, type, gene_id, gene_name, # score, phase,
                                                # cmp_ref, class_code, tss_id, exon_number
                                            )]) # contained_in, cmp_ref_gene, ref_gene_id, transcript_start, transcript_end,
best.merged.result_gff.compare.EX <- merge(best.merged.result_gff.compare.EX, 
                                           unique(best.merged.result_gff.compare[, !c("start", "end"), with = FALSE]),
                                           by=c('seqnames', 'transcript_id'))





stopifnot(nrow(best.merged.result_gff.compare) == length(unique(best.merged.result_gff.compare$transcript_id)))
## OK

#### Write out
fwrite(best.merged.result_gff.compare,     paste0(outdir, "/best.merged.result_gff.compare.tsv"),     sep = '\t')

