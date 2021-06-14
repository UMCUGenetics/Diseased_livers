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

wd <- paste0(base_dir,'/analysis3/compare_mut_profiles/')

pcawg_contribs_raw <- list(
   snv=read.csv(paste0(base_dir,'/analysis3/select_sigs/PCAWG_sigProfiler_SBS_signatures_in_samples.csv')),
   indel=read.csv(paste0(base_dir,'/analysis3/select_sigs/PCAWG_SigProfiler_ID_signatures_in_samples.csv')),
   dbs=read.csv(paste0(base_dir,'/analysis3/select_sigs/PCAWG_sigProfiler_DBS_signatures_in_samples.csv'))
)

calcSigFreq <- function(cancer.type='Liver-HCC'){
   pcawg_contribs <- lapply(pcawg_contribs_raw, function(i){
      #i[i$Cancer.Types %in% c('Liver-HCC','Biliary-AdenoCA'),]
      i[i$Cancer.Types==cancer.type,]
   })
   
   sig_freq <- lapply(names(pcawg_contribs), function(i){
      #i='snv'
      m <- pcawg_contribs[[i]][,-(1:3)]
      
      # apply(m,2,function(i){
      #    i2 <- i[i>0]
      #    out <- median(i2)
      # })
      
      df <- data.frame(
         n_with_contrib=colSums(m!=0),
         n_total=nrow(m)
      )
      df$frac_with_contrib <- df$n_with_contrib / df$n_total
      df <- cbind(sig_name=rownames(df), sig_type=i, df)
      rownames(df) <- NULL
      
      return(df)
   })
   sig_freq <- do.call(rbind, sig_freq)
   sig_freq$cancer_type <- cancer.type
   return(sig_freq)
}

sig_freqs <- rbind(
   calcSigFreq('Liver-HCC'),
   calcSigFreq('Biliary-AdenoCA')
)

write.table(
   sig_freqs, paste0(base_dir,'/analysis3/select_sigs/sig_freqs.txt'),
   sep='\t',
)
