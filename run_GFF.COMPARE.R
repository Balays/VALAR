

#TR.gff.compare.refmap <- fread('TORMAETAL.corr.v2_vs_LoRTIA.SO.TR.READS_GFF.COMPARE.LoRTIA_SO.TR_reads.gff2.refmap')

#TR.gff.compare.tmap <- fread('TORMAETAL.corr.v2_vs_LoRTIA.SO.TR.READS_GFF.COMPARE.LoRTIA_SO.TR_reads.gff2.tmap')

#TR.gff.compare.tracking <- fread('TORMAETAL.corr.v2_vs_LoRTIA.SO.TR.READS_GFF.COMPARE.tracking')




## NEM TUDOM HASZNÁLNI A GFF COMPARE-T SEM, ÖSSZEVISSZA ANNOTAÁLJA A NESTED TRANSZKRIPTEKET CSAK A TELJSE EGYEZÉST TALÁLJA MEG.


## IRNOM KELL EGY KÓDOT, AMI ÖSSZEVETI FOVERLAPS2-VEL A TRANSZKRIPTEKET ES A LEGNAGYOBB OVERLAP-T KIVALSZTJA

## VAGY CSINÁLOK EGY FOR CIKLUS A GFFCOMPARE-HEZ AHOL MEGADOM MINDEGYIK TRANSZKRIPETT KULON KÜLÖN ÉS AZTÁN UTAN VALSZTOM KI 
## AZ OUTPUOKBOL  A LEGKISEBB TAVOLSAGUT TRANSZKRIPTET

# NOIR-2-2
# UL15-splice-L
# TRS1-S-2
# nc-CTO-S


require(rtracklayer)
require(BiocParallel)
source('gff_compare_functions.R')

if (.Platform$OS.type=="windows") { NULL } else {safeBPParam(nproc)}


#

### Round 1.)
## run the function on all transcripts simultaneously

run_gff_compare(ref_mRNA = 'ALL_TRs', viral.mrna, query_gff = config$TR.reads.gfffile,
                out_dir = paste0(outdir, '/gffcompare_iterative'), 
                iterative = F
)
## this is ti check the annotation


### Round 2.)
## if everything OK,
## run the function on each transcript individually
ref_mRNAs <- unique(viral.mrna[,transcript_id])

all.merged.result_gff.compare.list <- bplapply(ref_mRNAs[],
                                               run_gff_compare,
                                               viral.mrna, query_gff = TR.reads.gfffile,
                                               out_dir = paste0(outdir, '/gffcompare_iterative')
)


all.merged.result_gff.compare <- rbindlist(all.merged.result_gff.compare.list)


### Round 3.)
ref_mRNAs <- unique(viral.mrna[,transcript_id])
## check the output files

out_dir <-paste0(outdir, '/gffcompare_iterative')
files <- data.table(filename=list.files(out_dir, '*.annotated.gtf'))
files[,transcript_id := gsub('.gffcompare.annotated.gtf', '', filename)]
files$size <- unlist(purrr::map(files$filename, function(x) as.integer(file.size(paste0(out_dir, '/', x)) / 1024)))

files <- merge(unique(viral.mrna[,.(transcript_id)]), files, by='transcript_id', all.y=T)
stopifnot(sum(is.na(files))==0)

## outputs that have no ID in the reference (maybe from previous run?)
files[!files[,transcript_id] %in% ref_mRNAs]

## reference transcripts that have no output
ref_mRNAs[!ref_mRNAs %in% files[,transcript_id]]

## run the analysis again on them
ref_mRNAs <- ref_mRNAs[!ref_mRNAs %in% files[,transcript_id]]
## also on those that have an empty output
ref_mRNAs <- c(ref_mRNAs, files[size == 0,transcript_id])

##
if (length(ref_mRNAs) > 0) {
  all.merged.result_gff.compare.list <- bplapply(ref_mRNAs[],
                                                 run_gff_compare,
                                                 viral.mrna, query_gff = config$TR.reads.gfffile,
                                                 out_dir = paste0(outdir, '/gffcompare_iterative')
  )
}

## check the outputs again
files <- data.table(filename=list.files(out_dir, '*.annotated.gtf'))
files[,transcript_id := gsub('.gffcompare.annotated.gtf', '', filename)]
files$size <- unlist(purrr::map(files$filename, function(x) as.integer(file.size(paste0(out_dir, '/', x)) / 1024)))

files <- merge(unique(viral.mrna[,.(transcript_id)]), files, by='transcript_id', all.y=T)
stopifnot(sum(is.na(files))==0)

### Finally:
## subset the reference transcripts to those that have an output
ref_mRNAs <- unique(viral.mrna[,transcript_id])
ref_mRNAs <- files[transcript_id %in% ref_mRNAs, transcript_id]

tr.miss   <-  setdiff(viral.mrna[,transcript_id], ref_mRNAs)

message('NR of Reference Transcripts: ',  length(unique(viral.mrna[,transcript_id])))
message('NR of Reference Transcripts with    an output: ', length(ref_mRNAs))
message('NR of Reference Transcripts without an output: ', length(tr.miss))


## check if all ref mRNAs have an output
if(length(tr.miss) > 0) {
  
  message('Reference Transcripts without an output: \n', paste(tr.miss, collapse = '\n'))
  
  ## Start again on tr.miss !!!
  
}



####  Part 2.) Import gff-compare results
## Testing:
#all.merged.result_gff.compare.list <- merge_gff_compare(ref_mRNA = ref_mRNAs[1],
#                                                        viral.mrna,
#                                                        out_dir = paste0(outdir, '/gffcompare_iterative')
#)


## Run the function in parallel
all.merged.result_gff.compare.list <- bplapply(ref_mRNAs[],
                                               merge_gff_compare,
                                               viral.mrna,
                                               out_dir = paste0(outdir, '/gffcompare_iterative')
)

gc()

saveRDS(all.merged.result_gff.compare.list, paste0(outdir, '/all.merged.result_gff.compare.list.rds'))

## check for the outputs- was there any errors?
res.class <- unlist(lapply(lapply(all.merged.result_gff.compare.list, class), paste, collapse='::'))
nrows     <- unlist(lapply(all.merged.result_gff.compare.list, nrow))
ncols     <- unlist(lapply(all.merged.result_gff.compare.list, ncol))

tr.miss   <- ref_mRNAs[!grepl('data.frame', res.class) | nrows == 0 | ncols != 40 ]
match(tr.miss, ref_mRNAs)

## run both functions again on erroneous outputs
if(length(tr.miss) != 0) {
  ##
  all.merged.result_gff.compare.list2 <- bplapply(tr.miss,
                                                 run_gff_compare,
                                                 viral.mrna, query_gff = config$TR.reads.gfffile,
                                                 out_dir = paste0(outdir, '/gffcompare_iterative')
  )
  
  all.merged.result_gff.compare.list2 <- bplapply(tr.miss[],
                                                 merge_gff_compare,
                                                 viral.mrna,
                                                 out_dir = paste0(outdir, '/gffcompare_iterative')
  )
  
}

## filter out the corrected results from all the results
all.merged.result_gff.compare.list <- all.merged.result_gff.compare.list[ -match(tr.miss, ref_mRNAs) ]

## put back the corrected
all.merged.result_gff.compare.list <- c(all.merged.result_gff.compare.list, all.merged.result_gff.compare.list2)

# Check length
stopifnot(length(ref_mRNAs) == length(all.merged.result_gff.compare.list))

# Reorder
tr.names       <- sapply(all.merged.result_gff.compare.list, function(x) unique(x[,'cmp_ref']))
tr.names.index <- match(ref_mRNAs, tr.names)
all.merged.result_gff.compare.list3 <- all.merged.result_gff.compare.list[tr.names.index]

# Check names
tr.names        <- unlist(sapply(all.merged.result_gff.compare.list3, function(x) unique(x[,'cmp_ref'])))
names(tr.names) <- NULL
stopifnot(all.equal(tr.names, ref_mRNAs))
##

## make DT
all.merged.result_gff.compare <- rbindlist(all.merged.result_gff.compare.list3)

## done




#### Part 3.)

## find the ref.TR with lowest 5-prime and 3-prime distance for each read.TR
all.merged.result_gff.compare[, min_dist_prime5 := min(abs(distance_prime5.tr)), by=transcript_id]
all.merged.result_gff.compare[, min_dist_prime3 := min(abs(distance_prime3.tr)), by=transcript_id]

## 
all.merged.result_gff.compare[min_dist_prime3 == abs(distance_prime3.tr), best_prime3 := T]
all.merged.result_gff.compare[min_dist_prime5 == abs(distance_prime5.tr), best_prime5 := T]

##
all.merged.result_gff.compare[, min_dist_junc := min(abs(abs_junc_distance.tr)), by=transcript_id]
all.merged.result_gff.compare[min_dist_junc == abs(abs_junc_distance.tr), best_junc := T]


gc()

#### write out
fwrite(all.merged.result_gff.compare, paste0(outdir, "/all.merged.result_gff.compare.tsv"), sep='\t')



## Move all files into gff-compre results directory
file.rename(list.files(outdir, '*.refmap', full.names = T), 
            gsub(paste0(outdir, '/'), paste0(outdir, '/gffcompare_iterative/'), list.files(outdir, '*.refmap', full.names = T)) ) 

file.rename(list.files(outdir, '*.tmap',  full.names = T),
            gsub(paste0(outdir, '/'), paste0(outdir, '/gffcompare_iterative/'), list.files(outdir, '*.tmap', full.names = T)) )



# Remove unnessary files?
file.remove( list.files(out_dir, '*.refmap', full.names = T) )
file.remove( list.files(out_dir, '*.tmap', full.names = T) )
file.remove( list.files(out_dir, '*.loci', full.names = T) )

# Remove all temporary files?
file.remove( list.files(out_dir, '*.refmap', full.names = T) )
file.remove( list.files(out_dir, '*.tmap', full.names = T) )
file.remove( list.files(out_dir, '*.loci', full.names = T) )









