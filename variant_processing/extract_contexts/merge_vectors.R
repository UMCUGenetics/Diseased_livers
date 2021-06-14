options(stringsAsFactors=F)

#========= Path prefixes =========#
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

devtools::load_all('/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/')

#========= Extract contexts from vcfs =========#
# vcf_paths <- read.xlsx(paste0(base_dir,'/metadata/vcf_paths_2020.xlsx'),sheet=1)
# 
matrices_dir <- paste0(base_dir,'/matrices/')
mut_types <- c('snv','indel','indel_pcawg','dbs','sv')


dirs <- paste0(matrices_dir,'/',mut_types,'/')
context_files <- lapply(dirs, list.files, full.names=T)
names(context_files) <- mut_types

contexts <- lapply(context_files, function(i){
   t(do.call(cbind,
      lapply(i, function(j){ read.delim(j, check.names=F) })
   ))
})

saveRDS(contexts, paste0(matrices_dir,'/merged_contexts.rds'))
write.table(
   do.call(cbind,contexts), 
   paste0(matrices_dir,'/merged_contexts.txt'),
   sep='\t',quote=F
)
