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

devtools::load_all('/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor')
#library(openxlsx)

#========= Extract contexts from vcfs =========#
matrices_dir <- paste0(base_dir,'/scripts/extract_sigs/matrices/')
dir.create(matrices_dir, showWarnings=F)

mut_types <- c('snv','indel','indel_pcawg','dbs','sv')
for(i in mut_types){
   dir.create(paste0(matrices_dir,'/',i), showWarnings=F)
}

write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

extractContextsAll <- function(vcf.som, vcf.sv, sample.name, out.dir=matrices_dir){
   
   #sample.name=vcf_paths$sample[1]
   #vcf.som=vcf_paths$som[1]
   #vcf.sv=vcf_paths$sv[1]
   
   if(dir.exists('/Users/lnguyen/')){
      vcf.som <- paste0('/Users/lnguyen/',vcf.som)
      vcf.sv <- paste0('/Users/lnguyen/',vcf.sv)
   }
   
   out_paths <- paste0(out.dir,'/',mut_types,'/',sample.name,'_',mut_types,'.txt')
   names(out_paths) <- mut_types
   
   contexts <- list()
   
   if(!file.exists(out_paths['snv'])){
      contexts$snv <- extractSigsSnv(
         vcf.som, output='contexts', vcf.filter='PASS', sample.name=sample.name, 
         keep.chroms=c(1:22,'X'), verbose=T
      )
      write.tsv(contexts$snv, out_paths['snv'])
   }
   
   if(!file.exists(out_paths['indel'])){
      contexts$indel <- extractSigsIndel(
         vcf.som, vcf.filter='PASS', sample.name=sample.name, 
         keep.chroms=c(1:22,'X'), verbose=T, method='CHORD'
      )
      write.tsv(contexts$indel, out_paths['indel'])
   }
   
   if(!file.exists(out_paths['indel_pcawg'])){
      contexts$indel_pcawg <- extractSigsIndel(
         vcf.som, vcf.filter='PASS', sample.name=sample.name, 
         keep.chroms=c(1:22,'X'), verbose=T, method='PCAWG'
      )
      write.tsv(contexts$indel_pcawg, out_paths['indel_pcawg'])
   }
   
   if(!file.exists(out_paths['dbs'])){
      contexts$dbs <- extractSigsDbs(
         vcf.som, vcf.filter='PASS', sample.name=sample.name, 
         keep.chroms=c(1:22,'X'), verbose=T
      )
      write.tsv(contexts$dbs, out_paths['dbs'])
   }
   
   if(!file.exists(out_paths['sv'])){
      contexts$sv <- extractSigsSv(
         vcf.sv, output='contexts', vcf.filter='PASS', sample.name=sample.name, 
         keep.chroms=c(1:22,'X'), verbose=T
      )
      write.tsv(contexts$sv, out_paths['sv'])
   }
   
}

args <- commandArgs(trailingOnly=T)
extractContextsAll(vcf.som=args[1], vcf.sv=args[2], sample.name=args[3])


# #========= Merge contexts =========#
# context_files <- list.files(matrices_dir, recursive=T, full.names=T)
# context_files <- lapply(mut_types, function(i){
#    grep(i, context_files, value=T)
# })
# names(context_files) <- mut_types
# 
# contexts <- lapply(context_files, function(i){
#    t(do.call(cbind,
#       lapply(i, function(j){ read.delim(j) })
#    ))
# })
# 
# saveRDS(contexts, paste0(matrices_dir,'/merged_contexts.rds'))
# write.tsv(do.call(cbind,contexts), paste0(matrices_dir,'/merged_contexts.txt'))




