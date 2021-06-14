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

devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor')

## Main ================================
manifest <- read.delim(paste0(base_dir,'/metadata/manifest.txt.gz'))
manifest <- subset(manifest, exp_group=='HCC_multibiopsy.raw')

paths <- manifest$som_vcf
if(dir.exists('/Users/lnguyen/')){
   paths <- paste0('/Users/lnguyen/',paths)
}
names(paths) <- manifest$sample_name

l_vcf <- lapply(paths, function(i){
   #i=paths[1]
   message('Reading: ',i)
   df <- readVcfFields(i)
   df[df$FILTER=='PASS',]
})

trunk_muts <- (function(){
   v <- lapply(l_vcf, function(i){
      paste(i$CHROM,i$POS,i$REF,i$ALT,sep='_')
   })
   v <- unlist(v, use.names=F)
   
   tab <- sort(table(v), decreasing=T)
   names(tab)[tab>=length(l_vcf)]
})()

l_vcf <- lapply(l_vcf, function(i){
   i$is_trunk_mut <- paste(i$CHROM,i$POS,i$REF,i$ALT,sep='_') %in% trunk_muts
   return(i)
})

## Subset ================================
l_vcf_2 <- lapply(l_vcf, function(i){
   df <- i[!i$is_trunk_mut,]
   df$is_trunk_mut <- NULL
   return(df)
})
names(l_vcf_2) <- paste0(names(l_vcf_2),'.non_trunk')

l_vcf_2$HCC_multibiopsy.trunk <- do.call(rbind, lapply(l_vcf, function(i){
   df <- i[i$is_trunk_mut,]
   df$is_trunk_mut <- NULL
   colnames(df)[ncol(df)] <- 'HCC_Trunk'
   return(df)
}))
l_vcf_2$HCC_multibiopsy.trunk <- l_vcf_2$HCC_multibiopsy.trunk[
   with(l_vcf_2$HCC_multibiopsy.trunk,{
      !duplicated(paste(CHROM,POS,REF,ALT,sep='_'))
   })
,]

header_vcf_paths <- c(paths, paths[[1]])
names(header_vcf_paths) <- names(l_vcf_2)

## Export ================================
writeVcf <- function(header.vcf, df, out.vcf){
   out.vcf.tmp <- paste0(out.vcf,'.tmp')
   
   get_header.sh <- sprintf("gzcat %s | grep ^## > %s", header.vcf, out.vcf.tmp)
   system(get_header.sh)
   
   #df <- l_vcf[[path.name]]
   #df <- df[df$is_trunk_mut==sel.trunk.muts,]
   
   colnames(df)[1] <- paste0('#',colnames(df)[1])
   suppressWarnings({
      write.table(df, out.vcf.tmp, sep='\t', quote=F, row.names=F, append=T)
   })
   
   compress.sh <- sprintf('gzip -c %s > %s && rm %s', out.vcf.tmp, out.vcf, out.vcf.tmp)
   system(compress.sh)
}

vcf_out_dir <- paste0(base_dir,'/scripts/filter_trunk_muts/vcf/')
dir.create(vcf_out_dir, showWarnings=F)

for(sample.name in names(l_vcf_2)){
   
   message(sample.name)
   
   writeVcf(
      header.vcf=header_vcf_paths[[sample.name]],
      df=l_vcf_2[[sample.name]],
      out.vcf=paste0(vcf_out_dir,'/',sample.name,'.vcf.gz')
   )
}
