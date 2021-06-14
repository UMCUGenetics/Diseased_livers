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

wd <- paste0(base_dir,'/analysis3/qc_checks/')

library(openxlsx)
library(naturalsort)

## Load data ================================
snp_genotyping <- read.xlsx(paste0(base_dir,'/scripts/snp_genotyping/genotypes.xlsx'), sheet=1)

metadata_footprint <- read.xlsx(paste0(base_dir,'/metadata/vcf_paths_2020.xlsx'))

## Gather QC ================================
paths <- metadata_footprint[,c('qc_metrics','sample_name')]
paths <- na.exclude(paths)
paths <- structure(paste0('/Users/lnguyen/',paths$qc_metrics), names=paths$sample_name)

qc_metrics <- do.call(rbind, lapply(names(paths), function(i){
   #message(i)
   df <- tryCatch({ read.delim(paths[[i]]) }, error=function(e){ return(NULL) })
   
   if(is.data.frame(df)){
      df <- cbind(sample_name=i, df)
   }
   
   df[,1:7]
}))

qc_metrics$is_ref <- ifelse(
   grepl('REF|BLOOD|(LIVER$)|(SPLEEN$)|(Adjacent$)|(t0)',qc_metrics$sample),
   'ref','case'
)

qc_metrics <- cbind(
   exp_group=metadata_footprint[match(qc_metrics$sample_name, metadata_footprint$sample_name),'exp_group'],
   qc_metrics
)

qc_metrics <- qc_metrics[
   naturalorder(paste(
      qc_metrics$exp_group,
      qc_metrics$sample_name,
      qc_metrics$is_ref,
      sep='_'
   ))
,]

qc_metrics$exp_group <- factor(qc_metrics$exp_group, unique(qc_metrics$exp_group))
qc_metrics$sample_name <- factor(qc_metrics$sample_name, unique(qc_metrics$sample_name))

qc_metrics$qc_pass <- TRUE
qc_metrics$qc_pass[qc_metrics$MEAN_COVERAGE<15] <- FALSE

write.table(
   qc_metrics,
   paste0(wd, '/wgs_qc_metrics.txt'),
   sep='\t', quote=F, row.names=F
)

## Blacklist samples ================================
snp_genotyping <- read.xlsx(paste0(base_dir,'/scripts/snp_genotyping/genotypes.xlsx'), sheet=1)

sample_blacklist <- c(
   as.character(qc_metrics$sample_name[!qc_metrics$qc_pass]),
   snp_genotyping$sample[snp_genotyping$sample_swap]
)

write.table(
   data.frame(sample_blacklist),
   paste0(wd, '/sample_blacklist.txt'),
   sep='\t', quote=F, row.names=F, col.names=F
)
