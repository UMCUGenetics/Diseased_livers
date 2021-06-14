options(stringsAsFactors=F)

## Path prefixes ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/')

wd <- paste0(base_dir,'/analysis3/compare_mut_profiles/')

WRITE_OUTPUT <- FALSE

## Main ================================
metadata <- read.delim(paste0(base_dir,'/analysis3/compare_mut_profiles/tables/S1_sample_metadata.txt'))
manifest <- read.delim(paste0(base_dir,'/metadata/manifest.txt.gz'))

manifest_ss <- manifest[match(metadata$sample_name, manifest$sample_name),]
write.table(
   manifest_ss,
   gzfile(paste0(base_dir,'/analysis3/variant_calling_output/manifest.ss.txt.gz')),
   sep='\t',quote=F,row.names=F
)
