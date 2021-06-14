#!/usr/bin/env Rscript
## Run on hpc

## Inputs ================================
args <- commandArgs(trailingOnly=T)
vcf.file <- args[1]
out.dir <- args[2]
sample.name <- args[3]

devtools::load_all('//hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/')

if(F){
   devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CUPs_classifier/processed/cuplr/geneDriverAnnotator/')
   vcf.file='/Users/lnguyen/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint_data/vcf/batch01//c110116R_ALC3CLONE32_ALC3BLOOD//somaticVariants/c110116RALC3BLOOD_c110116RALC3CLONE32/c110116RALC3BLOOD_c110116RALC3CLONE32_post_processed.vcf.gz'
   sample.name='ALC3_CLONE32'
   out.dir='/Users/lnguyen/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/scripts/annotate_som_variants/test/'
}

genes.bed.file=GENES_BED_FILE
exons.bed.file=EXONS_BED_FILE
java.path=JAVA_PATH
snpsift.path=SNPSIFT_PATH
snpeff.path=SNPEFF_PATH
verbose=T

dir.create(out.dir, showWarnings=F)

#bed_file <- read.delim(genes.bed.file, stringsAsFactors=F, check.names=F)

## Filter/annotate vcf ================================
vcf_ss_path <- paste0(out.dir,'/',sample.name,'.vcf.gz')
vcf_ann_path <- paste0(out.dir,'/',sample.name,'.ann.vcf.gz')
filt_ann_done <- paste0(out.dir,'/filt_ann.done')


if(!file.exists(filt_ann_done)){
   if(verbose){ message('> Filtering vcf on variants in genes in bed file and PASS variants...') }
   filterVcf(
      vcf.file = vcf.file,
      out.file = vcf_ss_path,
      mode = 'som',
      bed.file = NULL,
      java.path = java.path,
      snpsift.path = snpsift.path
   )
   
   if(verbose){ message('> Annotating variants...') }
   annotateVariantType(
      vcf.file = vcf_ss_path, 
      out.file = vcf_ann_path,
      genome = 'GRCh37.75',
      java.path = java.path, 
      snpeff.path = snpeff.path
   )
} else {
   if(verbose){ message('> Skipping filtering and annotating variants...') }
}


## Make txt ================================
all_done <- paste0(out.dir,'/all.done')

if(!file.exists(all_done)){
   txt_path <- paste0(out.dir,'/',sample.name,'.txt.gz')
   
   if(verbose){ message('> Extracting revelant fields in vcf to txt...') }
   extractVcfFields(
      vcf.file = vcf_ann_path,
      out.file = txt_path,
      java.path = java.path,
      snpsift.path = snpsift.path
   )
   
   txt <- read.delim(txt_path, stringsAsFactors=F)
   if(nrow(txt)!=0){
      txt$chrom <- as.character(txt$chrom)
      GenomeInfoDb::seqlevelsStyle(txt$chrom)<- 'NCBI'
   }
   
   # if(verbose){ message('> Subsetting for ENSG ids present in bed file...') }
   # txt <- txt[txt$ensembl_gene_id %in% bed_file$ensembl_gene_id,]
   # 
   # if(verbose){ message('> Adding HGNC gene ids...') }
   # #txt$hgnc_symbol <- ensgToHgncSymbol(txt$ensembl_gene_id)
   # txt$hgnc_symbol <- bed_file[match(txt$ensembl_gene_id, bed_file$ensembl_gene_id),'hgnc_symbol']
   
   if(nrow(txt)!=0){
      if(verbose){ message('> Subsetting for gene regions...') }
      txt <- txt[nchar(txt$snpeff_gene)!=0,]
      
      if(verbose){ message('> Adding ClinVar annotations...') }
      txt$clinvar_sig <- getClinSig(txt, CLINVAR_PATH) ## seqminer::tabix.read returns error if dataframe is empty
      txt$is_hotspot_mut <- detIsHotspotMut(txt, HOTSPOTS_PATH)
      
      #txt$clinvar_sig[ txt$is_hotspot_mut ] <- 'Pathogenic'
      
   } else {
      if(verbose){ message('> No variants remain after filtering. Skipping adding ClinVar annotations...') }
      txt <- cbind(txt, data.frame(clinvar_sig=character(), is_hotspot_mut=logical()))
   }
   
   if(verbose){ message('> Making mut_profike...') }
   txt <- mkMutProfileSnvIndel(
      txt,
      scoring=SCORING_MUT,
      keep.only.first.eff=T,
      include.hotspots=T,
      filter.no.impact.variants=F,
      verbose=verbose
   )
   
   write.table(txt, gzfile(txt_path), sep='\t', row.names=F, quote=F)
} else {
   if(verbose){ message('> Skipping making mut profile...') }
}

