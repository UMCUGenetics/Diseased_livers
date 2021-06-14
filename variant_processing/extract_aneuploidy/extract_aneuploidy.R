options(stringsAsFactors=F)

## Path prefixes ================================
base_dir <- list(
   hpc='/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/',
   mnt='/Users/lnguyen/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/',
   local='/Users/lnguyen/Documents/Luan_projects/Liver_footprint/'
)



for(i in base_dir){
   if(dir.exists(i)){ 
      base_dir <- i 
      break
   }
}


## Functions ================================
CENTRO_POS_HG19 <- c(
   chr1=123035434,chr2=93826171,chr3=92004854,chr4=51160117,chr5=47905641,
   chr6=60330166,chr7=59554331,chr8=45338887,chr9=48867679,chr10=40754935,
   chr11=53144205,chr12=36356694,chr13=17500000,chr14=17500000,chr15=18500000,
   chr16=36835801,chr17=23763006,chr18=16960898,chr19=26181782,chr20=27869569,
   chr21=12788129,chr22=14500000,chrX=60132012,chrY=11604553
)

getChromArm <- function(
   df=NULL, centro.pos=CENTRO_POS_HG19, one.armed.chroms=c(13,14,15,21,22),
   arm.only=F, seq.levels.style='NCBI', centro.intervals.rough.fix=F, show.warnings=F
){
   
   colnames(df)[1:3] <- c('chrom','start','end')
   
   df$chrom <- as.character(df$chrom)
   GenomeInfoDb::seqlevelsStyle(df$chrom)<- GenomeInfoDb::seqlevelsStyle(names(centro.pos))[1]
   
   one.armed.chroms <- as.character(one.armed.chroms)
   GenomeInfoDb::seqlevelsStyle(one.armed.chroms)<- GenomeInfoDb::seqlevelsStyle(names(centro.pos))[1]
   
   df$centro_pos <- centro.pos[ df$chrom ]
   
   df$is_pq <- df$start < df$centro_pos & df$end > df$centro_pos
   df$is_p <- df$start < df$centro_pos & df$end < df$centro_pos
   df$is_q <- df$start >= df$centro_pos & df$end >= df$centro_pos
   
   if(any(df$is_pq) & show.warnings){
      warning(sum(df$is_pq)," intervals overlap the centromere. Assigning as 'pq'")
   }
   df$arm <- c('pq','p','q')[
      max.col(df[,c('is_pq','is_p','is_q')])
   ]
   
   if(centro.intervals.rough.fix){
      df$arm[df$arm=='pq'] <- 'q'
   }
   
   df$arm[df$chrom %in% one.armed.chroms] <- 'q'
   
   if(!arm.only){
      GenomeInfoDb::seqlevelsStyle(df$chrom)<- seq.levels.style
      paste0(df$chrom,df$arm)
   } else {
      df$arm
   }
}

selectRequiredCols <- function(df, required.cols, sel.cols=NULL){
   #required.cols=c('ResolvedType','ClusterId','PosStart','PosEnd')
   
   if(!is.data.frame(df) & !is.matrix(df)){
      stop("`df` must be a dataframe or matrix")
   }
   
   if(is.null(sel.cols) || is.na(sel.cols)){
      sel.cols <- structure(required.cols, names=required.cols)
   }
   if(length(names(sel.cols))==0){
      stop("`sel.cols` must be a named character vector")
   }
   
   missing_names <- required.cols[ !(required.cols %in% names(sel.cols)) ]
   if(length(missing_names)>0){
      stop("`sel.cols` is missing the names: ", paste(missing_names, collapse=', '))
   }
   
   missing_cols <- sel.cols[!(sel.cols %in% colnames(df))]
   if(length(missing_cols)>0){
      stop("Some columns are missing in the input table: ",paste(missing_cols, collapse=', '))
   }
   
   df <- df[,sel.cols,drop=F]
   colnames(df) <- names(sel.cols)
   return(df)
}

calcChromArmPloidies <- function(
   ## I/O
   cnv.file=NULL, cnv=NULL, out.file=NULL, sel.cols=NULL,
   
   ## Params
   mode=c('total_cn','minor_cn'),
   min.rel.cum.segment.size=if(mode[1L]=='minor_cn'){ 0.9 } else { 0.5 },
   max.rel.cum.segment.size.diff=0.1,
   
   ## Misc
   keep.chroms=paste0('chr',c(1:22,'X')),
   one.armed.chroms=c(13,14,15,21,22),
   verbose=F
){
   
   ## Load data --------------------------------
   if(!is.null(cnv.file)){
      cnv <- read.delim(cnv.file, check.names=F, stringsAsFactors=F)
   }
   
   if(is.null(sel.cols)){
      sel.cols <- c(
         chrom='chromosome',start='start',end='end',
         total_cn='copyNumber'#,major_cn='majorAllelePloidy',minor_cn='minorAllelePloidy'
      )
   }
   
   cnv <- selectRequiredCols(
      df=cnv,
      required.cols=c(
         'chrom','start','end','total_cn'#,
         #'major_cn','minor_cn'
      ),
      sel.cols=sel.cols
   )
   
   GenomeInfoDb::seqlevelsStyle(cnv$chrom)<- 'NCBI'
   if(!is.null(keep.chroms)){
      GenomeInfoDb::seqlevelsStyle(keep.chroms)<- 'NCBI'
      cnv <- cnv[cnv$chrom %in% keep.chroms,]
   }
   
   ## Pre-calculations --------------------------------
   cnv$segment_size <- (cnv$end - cnv$start) + 1
   
   cnv$total_cn[cnv$total_cn < 0] <- 0 ## Ceiling negative values
   cnv$total_cn_int <- round(cnv$total_cn)
   
   if(mode[1L]=='minor_cn'){
      cnv$minor_cn[cnv$minor_cn < 0] <- 0
      cnv$minor_cn_int <- round(cnv$minor_cn)
   }

   
   ##----------------------------------------------------------------
   if(verbose){ message('Splitting CN segments by chrom arm...') }
   cnv$chrom_arm <- getChromArm(
      data.frame(chrom=cnv$chrom, start=cnv$start, end=cnv$end),
      one.armed.chroms=one.armed.chroms,
      centro.intervals.rough.fix=T
   )
   cnv$chrom_arm <- gsub('chr','',cnv$chrom_arm)
   cnv$chrom_arm <- factor(cnv$chrom_arm, unique(cnv$chrom_arm))
   cnv_split <- split(cnv, cnv$chrom_arm)
   
   # if(mode[1L]=='minor_cn'){
   #    minor_cn_segment_support <- lapply(cnv_split, function(i){
   #       #i=cnv_split$`Yp`
   #       df <- aggregate(i$segment_size, by=list(i$minor_cn_int), FUN=sum)
   #       colnames(df) <- c('minor_cn_int','cum_segment_size')
   #       df <- as.data.frame(lapply(df, as.integer)) ## aggregate can return lists instead of vectors as output
   #       
   #       df$cum_segment_size_rel <- df$cum_segment_size / sum(df$cum_segment_size)
   #       
   #       #df[which.max(df$cum_segment_size_rel),]
   #       if(nrow(df)!=0L){
   #          df <- df[order(df$cum_segment_size_rel, decreasing=T),]
   #          return(df)
   #       } else {
   #          df[1L,] <- NA
   #          return(df)
   #       }
   #    })
   # }
   
   
   ##----------------------------------------------------------------
   if(verbose){ message('Calculating copy number segment support...') }
   calcSegmentSupport <- function(colname){
      #colname='total_cn_int'
      lapply(cnv_split, function(i){
         #i=cnv_split$`1q`
         df <- aggregate(i$segment_size, by=list(i[,colname]), FUN=sum)
         colnames(df) <- c('cn','cum_segment_size')
         #df <- as.data.frame(lapply(df, as.integer)) ## aggregate can return lists instead of vectors as output
         
         df$cum_segment_size_rel <- df$cum_segment_size / sum(df$cum_segment_size)
         
         #df[which.max(df$cum_segment_size_rel),]
         if(nrow(df)!=0L){
            df <- df[order(df$cum_segment_size_rel, decreasing=T),]
            return(df)
         } else {
            df[1L,] <- NA
            return(df)
         }
      })
   }
   
   if(mode[1L]=='total_cn'){
      cn_segment_support <- calcSegmentSupport('total_cn_int')
   } else if(mode[1L]=='minor_cn'){
      cn_segment_support <- calcSegmentSupport('minor_cn_int')
   } else {
      stop("`mode` must be 'total_cn' or 'minor_cn'")
   }
   
   ## CN with most frequent total segment support is preliminary CN
   arm_cn_prelim <- unlist(lapply(cn_segment_support, function(i){ i[1L,'cn'] }))
   
   if( sum(!is.na(arm_cn_prelim)) / length(arm_cn_prelim) < 0.6 ){
      if(verbose){ message('Warning: CN data contains too many NAs. Returning a vector of NAs') }
      genome_cn <- NA ## Genome ploidy cannot be calculated if there is too little CN data
   } else {
      genome_cn <- as.numeric(names(
         sort(table(arm_cn_prelim),decreasing=T)
      )[1])
   }
   
   if(is.na(genome_cn)){
      arm_cn <- structure(rep(NA, length(cn_segment_support)), names=names(cn_segment_support))
   } else {
      if(verbose){ message('Calculating final arm ploidies...') }
      arm_cn <- unlist(lapply(cn_segment_support, function(i){
         #i=cn_segment_support[['12p']]
         
         if(is.na(i[1L,1L])){ return(NA) }
         
         ## E.g. >=50% of arm has CN of 2, then this is the CN
         if( i[1L,'cum_segment_size_rel'] >= min.rel.cum.segment.size ){
            return(i[1L,'cn'])
         }
         
         ## When multiple CNs have similar segment support as the one with the highest, if one
         ## of these have the same CN as the genome CN, return the genome CN. Otherwise, simply
         ## return the one with the highest segment support (as is done above)
         i$diffs <- i[1L,'cum_segment_size_rel'] - i[,'cum_segment_size_rel']
         cn_doubt <- i[i$diffs < max.rel.cum.segment.size.diff,'cn']
         
         if(any(cn_doubt==genome_cn)){
            return(genome_cn)
         } else {
            return(i[1L,'cn'])
         }
      }))
   }
   
   ploidies <- c(arm_cn, genome=genome_cn)
   
   ## Ensure consistent output (e.g. for when CN data is missing)  --------------------------------
   if(verbose){ message('Preparing output...') }
   chrom_arm_names <- paste0(rep(keep.chroms,each=2L),c('p','q'))
   chrom_arm_names <- chrom_arm_names[!(chrom_arm_names %in% paste0(one.armed.chroms,'p'))] ## Keep only q arm for one arm chromosomes
   chrom_arm_names <- c(chrom_arm_names,'genome')
   
   out <- structure(rep(NA,length(chrom_arm_names)), names=chrom_arm_names)
   out[names(ploidies)] <- ploidies
   
   if(is.null(out.file)){
      return(out)
   } else {
      out <- data.frame(chrom=names(out),ploidy=out,row.names=NULL)
      if(verbose){ message('Writing output...') }
      write.tsv(out, out.file)
   }
   
}

## Main ================================
manifest <- read.delim(paste0(base_dir,'/metadata/manifest.txt.gz'))

cnv_files <- manifest$cnv
if(dir.exists('/Users/lnguyen/')){
   cnv_files <- paste0('/Users/lnguyen/',cnv_files)
}
names(cnv_files) <- manifest$sample_name

counter <- 0
m_ploidy <- lapply(cnv_files, function(i){
   #i=cnv_files[1]
   #i=cnv_files[['ALC3_CLONE32']]
   
   counter <<- counter + 1
   message('[',counter,']: ', basename(i))
   
   cnv <- read.delim(i, check.names=F)
   colnames(cnv)[1] <- 'chromosome'
   out <- calcChromArmPloidies(cnv=cnv)
   #out <- t(out)
   return(out)
})

m_ploidy <- do.call(rbind, m_ploidy)
write.table(
   m_ploidy, gzfile(paste0(base_dir,'/scripts/extract_aneuploidy/m_ploidy.txt.gz')),
   sep='\t',quote=F
)
