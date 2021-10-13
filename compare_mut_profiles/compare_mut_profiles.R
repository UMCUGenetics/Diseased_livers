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

## Dependencies ================================
devtools::load_all('/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor')
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(nlme)
library(ggrepel)
library(cowplot)
library(matrixTests)
library(naturalsort)

## Helper functions
cacheAndReadData <- function(
   remote.path, local.path=NULL,
   cache.dir=path.expand('~/Documents/R_cache/'),
   overwrite=F
){
   
   #remote.path=paste0(base_dir,'/misc/processed/Chromatin_modifiers/scripts/annotate_svs/vis_sv_data_compact.txt.gz')
   
   ## Init --------------------------------
   ext <- c('.rds','.txt','.txt.gz','.csv','.csv.gz')
   regex <- gsub('[.]','[.]',ext)
   regex <- paste0(regex,'$')
   
   is_valid_ext <- sapply(regex, function(i){ grepl(i, remote.path) })
   if(all(!is_valid_ext)){
      stop('File extension must be one of the following:\n  ', paste(ext, collapse=', '))
   }
   
   ext <- ext[ which(is_valid_ext) ]
   regex <- regex[ which(is_valid_ext) ]
   
   set.seed(nchar(remote.path))
   
   ## Copy --------------------------------
   if(dir.exists('/hpc/')){
      local.path <- remote.path
   } else {
      if(is.null(local.path)){
         local.path <- paste0(
            cache.dir,
            sub(regex,'',basename(remote.path)),
            '.',paste(sample(letters, 8), collapse=''),ext
         )
      }
      
      local.path <- gsub('/+', '/', local.path)
      if(!file.exists(local.path) | overwrite){
         if(!file.exists(local.path)){
            message('Making local copy: ', local.path)
         } else {
            message('Updating local copy: ', local.path)
         }
         
         system(sprintf(
            'rsync -a %s %s',
            remote.path,
            local.path
         ))
      }
   }
   
   ## Read --------------------------------
   message('Reading local copy: ', local.path)
   if(ext=='.rds'){
      readRDS(local.path)
   } else if(grepl('txt',ext)){
      read.delim(local.path, check.names=F)
   } else {
      read.csv(local.path, check.names=F)
   }
}


## Metadata ================================
metadata_raw <- read.delim(paste0(base_dir,'/metadata/manifest.txt.gz'))

## Select and order samples --------------------------------
contexts_raw <- readRDS(paste0(base_dir,'/scripts/extract_sigs/matrices/merged_contexts.rds'))
contexts_raw$indel <- contexts_raw$indel_pcawg
contexts_raw$indel_pcawg <- NULL

sample_qc_blacklist <- read.delim(paste0(base_dir,'/analysis3/qc_checks/sample_blacklist.txt'))[,1]

SAMPLE_BLACKLIST <- unlist(c(
   sample_qc_blacklist, ## Sample swaps and low coverage
   
   ## Hypermutators
   names(which(rowSums(contexts_raw$snv)>=50000)), 
   names(which(rowSums(contexts_raw$indel)>=5000)),
   
   ## Hypomutators
   names(which(rowSums(contexts_raw$snv)<500)), 
   names(which(rowSums(contexts_raw$indel)<100))
))

DISEASED_LIVER_GROUPS <- c("Healthy","ALC","NASH","PSC","PCAWG_HCC","PCAWG_CCA")

metadata <- subset(
   metadata_raw, 
   !(sample_name %in% SAMPLE_BLACKLIST) & exp_group %in% DISEASED_LIVER_GROUPS, 
   select=-c(germ_vcf,som_vcf,sv_vcf,cnv_tsv)
)

#metadata <- metadata[naturalorder(metadata$sample_name),]

metadata$exp_group <- factor(metadata$exp_group, DISEASED_LIVER_GROUPS)
2 <- metadata[order(metadata$exp_group),]

metadata$sample_group <- factor(metadata$sample_group, unique(metadata$sample_group))
metadata$sample_name <- factor(metadata$sample_name, unique(metadata$sample_name))
metadata$sample_name_2 <- factor(metadata$sample_name_2, unique(metadata$sample_name_2))

## Add short sample name
metadata$sample_name_short <- (function(){
   v <- as.character(metadata$sample_name_2)
   for(i in unique(metadata$exp_group)){
      v <- gsub(i,'',v)
   }
   v <- gsub('^\\.','',v)
   #v <- gsub('CLONE','',v)
   return(v)
})()

## Add mut load
metadata <- (function(){
   m_mut_load <- do.call(cbind, lapply(contexts_raw, rowSums))
   colnames(m_mut_load) <- paste0('mut_load.',colnames(m_mut_load))
   cbind(
      metadata,
      m_mut_load[as.character(metadata$sample_name),]
   )
})()

## Write supp table
# write.table(
#    metadata, paste0(base_dir,'/analysis3/compare_mut_profiles/tables/S1_sample_metadata.txt'),
#    sep='\t', quote=F, row.names=F
# )

DISEASED_LIVER_COLORS <- c(
   Healthy='darkgrey',
   
   ## Set3
   ALC="#E41A1C",
   NASH="#377EB8",
   PSC="#4DAF4A",
   #HCC_multibiopsy="#984EA3",
   PCAWG_HCC="#FF7F00",
   PCAWG_CCA="#F781BF"
)
#brewer.pal(n=8, name='Set1')

## Metadata helper functions --------------------------------
getMetadata <- function(samples, keys=NULL, df=metadata, strings.as.factors=T){
   #samples=c('Healthy-CLONE8.1-PO','DRUP01010113T')
   #keys=c('sample_group','age')
   
   if(is.null(keys)){
      keys <- colnames(df)[colnames(df)!='sample_name']
   }
   
   main <- function(key){
      v <- df[match(samples, df$sample_name),key]
      if(strings.as.factors & (is.character(v) | is.factor(v))){
         exist_values <- unique(df[,key])
         exist_values <- exist_values[exist_values %in% v]
         v <- factor(v, exist_values)
      }
      return(v)
   }
   
   if(length(keys)==1){ 
      main(keys)
   } else {
      out <- as.data.frame(lapply(keys, main), stringsAsFactors=strings.as.factors)
      colnames(out) <- keys
      return(out)
   }
}

## Mutation contexts ###############################################################################
## Format data ================================
contexts <- lapply(contexts_raw, function(i){ i[as.character(metadata$sample_name),] })

MUT_SUBTYPES <- (function(){
   l <- lapply(contexts_raw, colnames)
   
   snv_types <- structure( gsub('(^\\w\\[)|(\\]\\w$)','',l$snv), names=l$snv )
   indel_types <- structure( gsub('[.]\\d+\\+*$','',l$indel), names=l$indel )
   sv_types <- structure( gsub('_.+$','',l$sv), names=l$sv )
   
   dbs_types <- mutSigExtractor::DBS_TYPES[ match(l$dbs, mutSigExtractor::DBS_TYPES$context),'ref' ]
   dbs_types <- paste0(dbs_types,'>NN')
   names(dbs_types) <- l$dbs
   
   c(snv_types, dbs_types, indel_types, sv_types)
})()

meltContexts <- function(contexts){
   df <- melt(contexts)
   colnames(df) <- c('sample','feature','count','mut_type')
   
   df$mut_subtype <- MUT_SUBTYPES[as.character(df$feature)]
   
   df$count_rel <- melt(lapply(contexts, function(i){ i/rowSums(i) }))$value
   df$count_rel[is.na(df$count_rel)] <- 0
   
   df <- as.data.frame(lapply(df, function(i){
      if(is.character(i)){ factor(i, unique(i)) }
      else if(is.factor(i)){ 
         i <- as.character(i) 
         factor(i, unique(i)) 
      }
      else { i }
   }))
   df <- cbind(
      df,
      getMetadata(df$sample, keys=c('sample_group','exp_group','cohort'))
   )
   
   return(df)
}

contexts_melt <- meltContexts(contexts)

## Calc no. of mutations found --------------------------------
# df <- subset(contexts_melt, cohort=='Footprint' & mut_type!='sv')
# sum(df$count)
# aggregate(df$count, list(df$exp_group), sum)

## Plot mut load ================================
## Prep data --------------------------------
calcMutLoad <- function(contexts_melt, age.na.fill=25){
   df <- aggregate(
      contexts_melt$count, 
      by=list(contexts_melt$sample, contexts_melt$mut_type),
      sum
   )
   colnames(df) <- c('sample','mut_type','count')
   
   #metadata_ss <- getMetadata(df$sample, keys=c('sample_group','exp_group','age','cohort'))
   #df <- cbind(df,metadata_ss)
   
   df <- cbind( df, getMetadata(df$sample) )
   
   df$age[is.na(df$age)] <- age.na.fill
   
   df <- do.call(rbind, lapply(split(df, df$exp_group), function(i){
      i$color_group <- as.character(as.integer(
         factor(i$sample_group, naturalsort(unique(i$sample_group)))
      ))
      return(i)
   }))
   
   row_order <- naturalorder(paste0(
      as.integer(factor(df$mut_type,levels=names(contexts))),
      as.character(df$sample)
   ))
   df <- df[row_order,]
   df <- do.call(rbind, split(df, df$exp_group))
   rownames(df) <- NULL
   df$sample <- factor(as.character(df$sample), unique(as.character(df$sample)))
   
   return(df)
}

mut_load <- calcMutLoad(contexts_melt)

# (function(){
#    df <- subset(mut_load, cohort=='Footprint')
#    aggregate(df$count, list(df$exp_group), sum)
#    sum(df$count)
# })()

## Total mut type load --------------------------------
p_mut_type_load <- (function(){
   pd <- subset(mut_load, exp_group %in% c('Healthy','ALC','NASH','PSC'))
   levels(pd$mut_type)[levels(pd$mut_type)=='snv'] <- 'sbs'
   levels(pd$mut_type) <- toupper(levels(pd$mut_type))
   
   ggplot(pd, aes(x=exp_group, y=count)) +
      facet_grid(mut_type~., scales='free_y') +
      scale_y_continuous(name='No. of mutations', limits=c(0,NA)) +
      scale_x_discrete(name='') +
      geom_boxplot(outlier.shape=NA) +
      geom_jitter(height=0, width=0.1) +
      theme_bw()
})()

if(WRITE_OUTPUT){
   pdf(paste0(wd,'/plots/mut_type_load.pdf'), 8, 6)
   plot(p_mut_type_load)
   dev.off()
}


## Mut rate per year --------------------------------
calcMutsPerYr <- function(mut_load){
   ## For one mut type:
   ##   Fit linear model
   ##   Compare mut accumulation in disease vs healthy liver
   
   main <- function(df){
      #df=subset(mut_load, mut_type=='snv')
      df_split <- split(df, df$exp_group)
      df_split <- df_split[sapply(df_split,nrow)!=0]
      lm_summ <- lapply(df_split, function(i){
         #i=df_split[[1]]
         #model_lm <- lm(count ~ 0 + age, i)
         model <- lme(count~0+age, data=i,random=~+1|sample_group)
         
         summ <- summary(model)
         ci <- intervals(model, 0.95, which='fixed')$fixed
         
         out <- cbind(
            summ$tTable,
            ci[,'lower'],
            ci[,'upper'],
            as.numeric(VarCorr(model)['(Intercept)','Variance'])
         )
         colnames(out) <- c('slope','se','df','tvalue','pvalue','ci_lower','ci_upper','variance')
         return(out)
      })
      lm_summ <- as.data.frame(do.call(rbind, lm_summ))
      lm_summ$n <- lm_summ$df+1
      
      rownames(lm_summ) <- unique(df$exp_group)
      lm_summ$pvalue <- NULL ## This is the pvalue for the lines themselves, not the z-test between healthy vs disease
      
      ## Calc z statistic and pvalues, Disease vs. Healthy
      ## Adapted from:
      ## https://stats.stackexchange.com/questions/435644/is-there-a-method-to-look-for-significant-difference-between-two-linear-regressi
      ## 2-sided z-test
      z_test <- function(b1, se1, b2, se2){
         2*pnorm(-abs(
            (b1-b2)/sqrt(se1^2+se2^2)
         ))
      }

      lm_summ$z_pvalue <- z_test(
         b1= lm_summ[,'slope'],
         se1=lm_summ[,'se'],
         b2= lm_summ['Healthy','slope'],
         se2=lm_summ['Healthy','se']
      )
      
      ## Two-tailed F-test, Disease vs. Healthy
      ## Adapted from:
      ## https://www.statisticshowto.com/probability-and-statistics/hypothesis-testing/f-test/
      ## https://www.statology.org/f-critical-value-r/
      f_test <- function(variance1, df1, variance2, df2, alternative='two.sided'){
         # variance1=lm_summ[,'variance']
         # df1=lm_summ[,'df']
         # 
         # variance2=lm_summ['Healthy','variance']
         # df2=lm_summ['Healthy','df']
         
         # ## Make sure that higher variance divides by lower variances
         # variances <- cbind(variance1, variance2)
         # variances <- t(apply(variances,1,sort,decreasing=T))
         # f_statistic <- variances[,1] / variances[,2]
         
         ## two tailed F-test, see code for `stats:::var.test.default()`
         f_statistic <- variance1 / variance2
         f_pvalue <- pf(f_statistic, df1=df1, df2=df2)
         
         if(alternative=='two.sided'){
            f_pvalue <- 2*pmin(f_pvalue, 1-f_pvalue) 
         } else if(alternative=='greater'){
            f_pvalue <- 1-f_pvalue
         }
         cbind(f_statistic, f_pvalue)
      }
      
      # ## Test if F-test code is correct
      # x <- rnorm(50, mean = 0, sd = 2)
      # y <- rnorm(30, mean = 1, sd = 1)
      # var.test(x,y, alternative='greater')
      # 
      # var_x <- var(x)
      # var_y <- var(y)
      # df_x <- length(x)-1
      # df_y <- length(y)-1
      # f_test(var_x, df_x, var_y, df_y, alternative='greater')
      
      lm_summ <- cbind(
         lm_summ,
         f_test(
            variance1=lm_summ[,'variance'], 
            df1=lm_summ[,'df'], 
            variance2=lm_summ['Healthy','variance'], 
            df2=lm_summ['Healthy','df'],
            alternative='greater'
         )
      )
      
      ## Z-test power analysis
      ## Cohen's d = abs(mean1-mean2)/pooled_sd;
      # cohens_d <- function(n1, b1, var1, n2, b2, var2){
      #    ## Pooled variance: https://www.statology.org/pooled-variance-in-r/
      #    pooled_var <- ((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2)
      #    abs(b1-b2) / sqrt(pooled_var)
      # }
      # 
      # lm_summ$cohens_d <- cohens_d(
      #    n1=  lm_summ[,'n'],
      #    b1=  lm_summ[,'slope'],
      #    var1=lm_summ[,'variance'],
      #    n2=  lm_summ['Healthy','n'],
      #    b2=  lm_summ['Healthy','slope'],
      #    var2=lm_summ['Healthy','variance']
      # )
      # 
      # pwr::cohen.ES(test='t', size='large')
      # 
      # pwr::pwr.t2n.test(n1=5, n2=8, sig.level=0.05, power=0.8, alternative='greater')
      # pwr::pwr.t2n.test(n1=3, n2=8, sig.level=0.05, power=0.8, alternative='greater')
      # 
      # 
      # pwr::pwr.t.test(
      #    n=NULL,
      #    d=0.5,
      #    sig.level=0.05, power=0.8, alternative='greater'
      # )

      
      return(lm_summ)
   }
   
   l_lm_summ <- split(mut_load, mut_load$mut_type)
   l_lm_summ <- l_lm_summ[sapply(l_lm_summ,nrow)!=0]
   l_lm_summ <- lapply(l_lm_summ, main)
   
   l_lm_summ <- lapply(names(l_lm_summ), function(i){
      lm_summ <- l_lm_summ[[i]]
      cbind(mut_type=i, exp_group=rownames(lm_summ), lm_summ)
   })
   
   out <- do.call(rbind, l_lm_summ)
   out$mut_type <- factor(out$mut_type, levels(mut_load$mut_type))
   out$exp_group <- factor(out$exp_group, levels(mut_load$exp_group))
   
   rownames(out) <- NULL
   
   return(out)
}
if(F){
   param_grid <- expand.grid(
      c('Healthy','ALC','NASH','PSC'),
      c('snv','dbs','indel','sv')
   )
   colnames(param_grid) <- c('exp_group','mut_type')
   do.call(rbind, lapply(1:nrow(param_grid), function(i){
      #i=1
      i_exp_group <- as.character(param_grid$exp_group[i])
      i_mut_type <- as.character(param_grid$mut_type[i])
      
      df <- subset(mut_load, exp_group==i_exp_group & mut_type==i_mut_type)
      group_means <- aggregate(df$count, list(df$sample_group), mean)$x
      test_out <- shapiro.test(group_means)
      
      # test_out <- shapiro.test(
      #    subset(mut_load, exp_group==i_exp_group & mut_type==i_mut_type, count, drop=T)
      # )
      data.frame(
         exp_group=i_exp_group,
         mut_type=i_mut_type,
         pvalue=test_out$p.value
      )
   }))
   

}

plotMutsPerYr <- function(df_raw, show.sample.numbers=F, exp.group.colors=NULL, drop.legend.levels=F){
   #exp.group.colors=DISEASED_LIVER_COLORS
   #df_raw <- subset(mut_load, mut_type %in% c('snv','indel') & (cohort=='Footprint' | exp_group=='Healthy') & exp_group!='HCC_multibiopsy')
   
   ## Main --------------------------------
   levels(df_raw$mut_type)[levels(df_raw$mut_type)=='snv'] <- 'sbs'
   levels(df_raw$mut_type) <- toupper(levels(df_raw$mut_type))
   
   df <- df_raw
   df$exp_group <- getMetadata(df$sample, 'exp_group')
   
   df_split <- split(df, df$exp_group)
   df_split <- df_split[sapply(df_split, nrow)!=0]
   
   ## Pair data for healthy group with all disease groups
   df <- do.call(rbind, lapply( df_split[names(df_split)!='Healthy'], function(i){
      out <- rbind(i, df_split[['Healthy']])
      out$plot_group <- out$exp_group[1]
      return(out)
   }))
   
   ## Convert sample groups to color number code
   df_split <- split(df, df$exp_group)
   df <- do.call(rbind, lapply(df_split, function(i){
      #i=df_split[[1]]
      i <- i[order(i$age),]
      i$color_group <- as.character(as.integer(
         factor(i$sample_group, unique(i$sample_group))
      ))
      return(i)
   }))
   rownames(df) <- NULL
   
   df$exp_group <- factor(df$exp_group, names(exp.group.colors))
   df$color_label <- df$color_group
   df$color_label[ df$cohort=='PCAWG' ] <- ''
   
   ## Sumary stats --------------------------------
   # fit <- calcMutsPerYr(df_raw)
   # fit <- fit[fit$exp_group %in% df$exp_group,]
   # fit$exp_group <- droplevels(fit$exp_group)
   # 
   # fit$ci_error <- fit$slope - fit$ci_lower
   # 
   # max_age <- max(df$age)
   # df_line_tmp <- fit[,c('mut_type','exp_group','slope','ci_error')]
   # df_line <- rbind(
   #    within(df_line_tmp,{ 
   #       age <- 0
   #       count <- 0
   #       age_lower <- -ci_error
   #       age_upper <- ci_error
   #    }),
   #    
   #    within(df_line_tmp,{ 
   #       age <- max_age
   #       count <- max_age * slope 
   #       age_lower <- age-ci_error
   #       age_upper <- age+ci_error
   #    })
   # )
   # 
   # df_line_split <- split(df_line, df_line$exp_group)
   # df_line_split <- df_line_split[names(df_line_split)!='Healthy']
   # df_line <- lapply(names(df_line_split), function(i){
   #    out <- rbind(
   #       subset(df_line, exp_group=='Healthy'),
   #       df_line_split[[i]]
   #    )
   #    out$plot_group <- i
   #    return(out)
   # })
   # df_line <- do.call(rbind, df_line)
   # rm(df_line_split)
   # 
   # ggplot(df, aes(x=age, y=count, group=exp_group)) +
   #    facet_grid(mut_type~plot_group, scales='free_y') +
   #    
   #    geom_line(data=df_line, aes(color=exp_group)) +
   #    geom_ribbon(data=df_line, aes(fill=exp_group, xmin=age_lower, xmax=age_upper), alpha=0.5) +
   #    
   #    geom_point(aes(fill=exp_group), shape=21, size=1.2) +
   #    
   #    theme_bw()
   

   
   ######
   fit <- calcMutsPerYr(df_raw)
   fit <- fit[fit$exp_group %in% df$exp_group,]
   fit$exp_group <- droplevels(fit$exp_group)
   
   fit$ci_error <- fit$slope - fit$ci_lower
   
   ## Rounding
   fit$slope_string <- sub('0+$','',signif(fit$slope,3))
   fit$n_decimal_places <- nchar(sub('[.]','',stringr::str_extract(fit$slope_string, '[.]\\d+$')))
   fit$n_decimal_places[is.na(fit$n_decimal_places)] <- 0
   
   fit$ci_error_string <- round(fit$ci_error,fit$n_decimal_places)
   fit$ci_error_string <- sub('0+$','',fit$ci_error_string )
   
   fit$label <-  (function(){
      l <- split(fit, droplevels(fit$mut_type))
      unlist(lapply(l, function(i){
         #i=l[[1]]
         paste0(
            'Muts. / year:\n',
            ' ',i$exp_group,': ', i$slope_string,' ±',i$ci_error_string,'\n',
            ' Healthy: ', i$slope_string[1],' ±',i$ci_error_string[1],'\n',
            '\n',
            'Z-test: p = ',signif(i$z_pvalue,2),'\n',
            'F-test: p = ',signif(i$f_pvalue,2)
         )
         
      }))
   })()
   #fit <- fit[fit$exp_group!='Healthy',]
   #levels(fit$exp_group) <- levels(df$exp_group)
   fit$plot_group <- fit$exp_group
   
   non_healthy_levels <- levels(fit$plot_group)[levels(fit$plot_group)!='Healthy']
   fit_non_healthy <- do.call(rbind, lapply(non_healthy_levels, function(i){
      fit_ss <- fit[fit$plot_group=='Healthy',]
      fit_ss$label <- ''
      fit_ss$plot_group <- i
      return(fit_ss)
   }))
   
   fit <- rbind(
      fit[fit$plot_group!='Healthy',],
      fit_non_healthy
   )
   label_ypos <- aggregate(df$count, list(df$mut_type), max)
   label_ypos <- structure(label_ypos$x, names=as.character(label_ypos$Group.1))
   fit$label_ypos <- label_ypos[as.character(fit$mut_type)]
   
   ## Plot --------------------------------
   p <- ggplot(df, aes(x=age, y=count, group=exp_group)) +
      facet_grid(mut_type~plot_group, scales='free_y') +
      
      geom_abline(data=fit, aes(intercept=0, slope=slope, color=exp_group), show.legend=F) +
      geom_point(aes(fill=exp_group), shape=21, size=1.2) +
      geom_text(data=fit, aes(x=0, y=label_ypos, label=label), hjust=0, vjust=1, size=2.5) +
      
      #geom_ribbon(aes(ymin=slope*age-ci_error, ymax=slope*age+ci_error, fill=exp_group), alpha=0.5) +
      
      #xlim(0, NA) + 
      #xlim(0, max(AGES, na.rm=T)) + 
      xlim(0, max(df$age, na.rm=T)+5) + 
      ylim(0,NA) +
      ylab('# mutations') + xlab('Age (years)') +
      
      theme_bw()
   
   if(!is.null(exp.group.colors)){
      p <- p + 
         scale_fill_manual(values=exp.group.colors, name='Group', drop=drop.legend.levels) +
         scale_color_manual(values=exp.group.colors, name='Group', drop=drop.legend.levels)
   }
   
   if(show.sample.numbers){
      p <- p + 
         geom_text_repel(
            aes(label=sample_name_short, color=exp_group), 
            max.overlaps=20, size=2, segment.size=0.3, show.legend=F
         )
   }
   
   return(p)
}

if(WRITE_OUTPUT){
   ## Footprint only
   pdf(paste0(wd,'/plots/muts_per_yr.footprint.pdf'), 8, 4.5)
   plotMutsPerYr(
      subset(mut_load, mut_type %in% c('snv','indel') & (cohort=='Footprint' | exp_group=='Healthy') & exp_group!='HCC_multibiopsy'),
      exp.group.colors=DISEASED_LIVER_COLORS,
      show.sample.numbers=T, drop.legend.levels=T
   )
   plotMutsPerYr(
      subset(mut_load, mut_type %in% c('dbs','sv') & (cohort=='Footprint' | exp_group=='Healthy') & exp_group!='HCC_multibiopsy'),
      exp.group.colors=DISEASED_LIVER_COLORS,
      show.sample.numbers=T, drop.legend.levels=T
   )
   dev.off()
}



## Plot contexts ================================
MUT_SUBTYPE_COLORS <- c(
   "C>A"="#06BAEB",
   "C>G"="#737373",
   "C>T"="#E12825",
   "T>A"="grey",
   "T>C"="#A1CD63",
   "T>G"="#EDC4C5",
   
   "AC>NN"="#03BCEC",
   "AT>NN"="#0165CB",
   "CC>NN"="#9FCC62",
   "CG>NN"="#006401",
   "CT>NN"="#FE9796",
   "GC>NN"="#E22823",
   "TA>NN"="#FCB066",
   "TC>NN"="#FC8100",
   "TG>NN"="#CC98FD",
   "TT>NN"="#4C0199",
   
   "del.1.C"="#F3BF7B",
   "del.1.T"="#EE8633",
   "ins.1.C"="#B9DA94",
   "ins.1.T"="#569D40",
   "del.2.rep"="#F5CBB7",
   "del.3.rep"="#EE8E72",
   "del.4.rep"="#DE523E",
   "del.5+.rep"="#AD2D25",
   "ins.2.rep"="#D2E0EF",
   "ins.3.rep"="#9DC3DC",
   "ins.4.rep"="#5D97C4",
   "ins.5+.rep"="#2E65A5",
   "del.2.mh"="#E1E1ED",
   "del.3.mh"="#B5B6D5",
   "del.4.mh"="#8585B8",
   "del.5+.mh"="#5D4494",
   "del.mh"="#5D4494",
   
   "DEL"="#9DD1C7",
   "DUP"="#FFFCBB",
   "INV"="#BDBAD7",
   "TRA"="#EB8677"
)

plotContexts <- function(
   df, value.col='count_rel', ctrl.group='Healthy',exp.group.order=DISEASED_LIVER_GROUPS, rel.widths=c(9,1)
){
   #df=subset(contexts_melt_ss, mut_type=='indel')
   df$feature <- factor(df$feature, unique(df$feature))
   
   ## Format data for bar plots
   df_agg <- aggregate(
      df[[value.col]], 
      list(feature=df$feature,exp_group=df$exp_group), 
      function(x){ 
         c(median=median(x), min=min(x), max=max(x), q1=unname(quantile(x,0.25)), q3=unname(quantile(x,0.75)))
      }
   )
   
   df_agg <- cbind(
      df_agg[,1:2],
      as.data.frame(df_agg$x)
   )
   
   df_agg$mut_subtype <- df[match(df_agg$feature, df$feature),'mut_subtype']
   df_agg$feature <- factor(df_agg$feature, unique(df$feature))
   
   ## Wilcox tests vs Healthy
   matrices <- lapply(unique(df$exp_group), function(i){
      df_ss <- df[df$exp_group==i,]
      m <- dcast(df_ss[,c('sample','feature',value.col)], sample~feature, value.var=value.col)
      rownames(m) <- m[,1]
      m[,1] <- NULL
      m <- as.matrix(m)
      return(m)
   })
   names(matrices) <- unique(df$exp_group)
   
   case_groups <- names(matrices)[names(matrices)!=ctrl.group]
   suppressWarnings({
      wilcox_tests <- do.call(rbind, lapply(case_groups, function(i){
         #i=case_groups[[1]]
         out <- col_wilcoxon_twosample(matrices[[i]], matrices[[ctrl.group]], alternative='greater')
         data.frame(
            feature=rownames(out), exp_group=i, 
            pvalue=out$pvalue, qvalue=p.adjust(out$pvalue, 'bonferroni')
         )
      }))
   })
   
   df_agg$qvalue <- wilcox_tests[
      match(
         paste0(df_agg$exp_group,':',df_agg$feature),
         paste0(wilcox_tests$exp_group,':',wilcox_tests$feature)
      )
      ,'qvalue']
   df_agg$qvalue[is.na(df_agg$qvalue)] <- 1
   df_agg$signif <- ''
   df_agg$signif[df_agg$qvalue<0.01] <- '*'
   
   ##
   yvalue_offset <- max(df_agg$median) / 200
   df_agg$median[df_agg$median==0] <- yvalue_offset
   
   if(!is.null(exp.group.order)){
      df_agg$exp_group <- factor(
         df_agg$exp_group, 
         exp.group.order[exp.group.order %in% df_agg$exp_group]
      )
   } else {
      df_agg$exp_group <- factor(df_agg$exp_group, unique(df_agg$exp_group))
   }
   
   ##
   df_agg$mut_subtype_2 <- df_agg$mut_subtype
   levels(df_agg$mut_subtype_2)[grep('del.*mh',levels(df_agg$mut_subtype_2))] <- 'del.mh'
   
   ## Plot contexts
   p1 <- ggplot(df_agg, aes(x=feature, y=median)) +
      facet_grid(exp_group~mut_subtype_2, scales='free_x', space='free_x') +
      
      geom_bar(aes(fill=mut_subtype),stat='identity') +
      scale_fill_manual(values=MUT_SUBTYPE_COLORS, name='Mut. subtype', guide=F) +
      geom_linerange(aes(ymin=q1, ymax=q3)) +
      geom_text(aes(y=0, label=signif), vjust=0.25, size=5) +
      
      ylab('Relative counts\n(median, q1, q3)') +
      xlab('Mutation context') +
      theme_bw() +
      theme(
         panel.grid.minor.y=element_blank(),
         axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5),
         panel.spacing.x=unit(0, "lines"),
         legend.position='left'
      )
   
   ## Plot cos sim
   cosSim <- function(x, y) { x %*% y / (sqrt(x %*% x) * sqrt(y %*% y)) }
   
   medians <- split(df_agg$median, df_agg$exp_group)
   cos_sims <- sapply(medians[names(medians)!=ctrl.group], function(i){
      cosSim(medians[[ctrl.group]], i)
   })
   cos_sims <- c(NA, cos_sims)
   names(cos_sims)[1] <- ctrl.group
   
   cos_sims <- data.frame(
      exp_group=factor(names(cos_sims), rev(names(cos_sims))), 
      cos_sim=cos_sims,
      row.names=NULL
   )
   
   p2 <- ggplot(cos_sims, aes(x=1, y=exp_group)) +
      geom_tile(aes(fill=cos_sim), color='black') +
      geom_text(aes(label=round(cos_sim,2))) +
      scale_fill_distiller(palette='YlGnBu', direction=-1, limits=c(0,1),  na.value='white', guide=F) +
      scale_y_discrete(expand=c(0,0), paste0('Cosine similarity vs. ', ctrl.group), position='right') +
      scale_x_discrete(expand=c(0,0)) +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         #axis.title.y=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank()
      )
   
   p_combined <- plot_grid(p1,p2, align='h', axis='tblr', nrow=1, rel_widths=rel.widths)
   return(p_combined)
}

## Export plots ================================
if(WRITE_OUTPUT){
   p_contexts <- (function(){
      
      df <- contexts_melt
      
      df <- df[order(df$exp_group),]
      df$exp_group <- as.character(df$exp_group)
      df$exp_group[df$sample_group=='HCC_multibiopsy.trunk'] <- 'HCC_multibiopsy\ntrunk'
      df$exp_group[df$sample_group=='HCC_multibiopsy.non_trunk'] <- 'HCC_multibiopsy\nnon_trunk'
      df$exp_group <- factor(df$exp_group, unique(df$exp_group))
      
      l <- lapply(levels(df$mut_type), function(mut.type){
         #mut.type='dbs'
         pd <- subset(df, mut_type==mut.type)
         plotContexts(
            pd, value.col='count_rel',
            ctrl.group='Healthy', exp.group.order=levels(df$exp_group)
         )
      })
      return(l)
   })()
   
   pdf(paste0(wd,'/plots/mut_contexts_rel.pdf'),13,10)
   for(i in p_contexts){ plot(i) }
   dev.off()
}


## Signatures ######################################################################################
##  LSQ fit --------------------------------
sig_metadata <- read.delim(mutSigExtractor::SIG_METADATA_PATH)
sig_freqs <- read.delim(paste0(base_dir,'/analysis3/select_sigs/sig_freqs.txt'))

SEL_SIGS <- unique(sig_freqs$sig_name[sig_freqs$frac_with_contrib>0.1])
#sig_freqs[sig_freqs$frac_with_contrib>0.1,]

fitToSignaturesWrapper <- function(contexts){
   #mut.context.counts=contexts$snv

   mut_types <- c('snv','dbs','indel')

   sig_profiles <- list(
      snv=SBS_SIGNATURE_PROFILES_V3,
      dbs=DBS_SIGNATURE_PROFILES,
      indel=INDEL_SIGNATURE_PROFILES
   )
   sig_profiles <- lapply(sig_profiles, function(i){ i[,colnames(i) %in% SEL_SIGS] })

   sig_contribs <- lapply(mut_types, function(i){
      fitToSignatures(
         mut.context.counts=contexts[[i]],
         signature.profiles=sig_profiles[[i]],
         #max.delta=0.001,
         use.lsq.r=F
      )
   })
   names(sig_contribs) <- mut_types

   return(sig_contribs)
}

sig_contribs <- fitToSignaturesWrapper(contexts)

hclustCustom <- function(m){
   #m=sig_contribs$snv
   
   sample_info <- data.frame(
      sample=rownames(m),
      getMetadata(rownames(m))
   )
   
   m_rel <- m/rowSums(m)
   m_rel[is.na(m_rel)] <- 0
   m_rel <- as.data.frame(m_rel)
   
   l <- split(m_rel, sample_info$exp_group)
   out <- lapply(l, function(i){
      hc <- hclust(dist(i))
      rownames(i)[hc$order]
   })
   unlist(out, use.names=F)
}

plotSigContrib <- function(
   m, sample.order=NULL, show.etiology=F, rel.contrib=T, 
   facet.space='fixed', strip.text.angle=0, strip.text.hjust=0.5,
   show.sample.names=T, use.short.sample.names=T, group.max.n.samples=20, fill.palette='Set3'
){
   #m=sig_contribs$snv

   m_rel <- m/rowSums(m)
   m_rel[is.na(m_rel)] <- 0

   ## To long form df
   df <- melt(as.matrix(m_rel))
   colnames(df) <- c('sample','sig','contrib_rel')

   ## Filter low contrib sigs
   df$contrib_abs <- melt(as.matrix(m))$value
   # df$sig <- as.character(df$sig)
   # df$sig[df$contrib_abs<50] <- 'Unassigned'
   # df$sig <- factor(df$sig, c(colnames(m),'Unassigned'))

   ## Sample order
   if(is.null(sample.order)){
      hc <- hclust(dist(m_rel))
      df$sample <- factor(df$sample, rownames(m_rel)[hc$order])
   } else {
      if(!all(sample.order %in% df$sample)){
         stop('Sample names in `df` are missing in `sample.order`')
      }
      df$sample <- factor(df$sample, sample.order)
   }

   ## Add metadata
   df <- cbind(df, getMetadata(df$sample))
   df$exp_group <- factor(df$exp_group, DISEASED_LIVER_GROUPS)
   df <- df[order(df$exp_group),]

   #subset(df, cohort=='Footprint' & sig=='SBS4')
   #subset(df, cohort=='Footprint' & sig=='ID3')

   if(show.etiology){
      df$etiology <- sig_metadata$sig_etiology[ match(df$sig, sig_metadata$sig_name) ]
      df$sig <- paste0(df$sig,'.',df$etiology)
      df$sig <- factor(df$sig, unique(df$sig))
   }

   ## Plot params
   which_contrib <- if(rel.contrib){ 'contrib_rel' } else { 'contrib_abs' }

   exp_group_counts <- unlist(lapply(split(df$sample, df$exp_group), function(i){
      length(unique(i))
   }))
   df$exp_group_counts <- unname(exp_group_counts[ as.character(df$exp_group) ])
   
   if(!use.short.sample.names){
      df$xlabs <- as.character(df$sample)
   } else {
      df$xlabs <- as.character(df$sample_name_short)
   }
   
   df$xlabs[df$exp_group_counts>group.max.n.samples] <- ''

   ## Plot main
   p <- ggplot(df, aes_string(x='sample', y=which_contrib)) +

      facet_grid(.~exp_group, scales='free_x', space=facet.space) +
      geom_bar(aes(fill=sig),stat='identity', position='stack', width=1) +

      scale_fill_brewer(palette=fill.palette, name='Signature') +
      scale_y_continuous(expand=c(0,0), name='Relative contribution') +
      scale_x_discrete(expand=c(0,0), name='Sample', labels=df$xlabs, breaks=df$sample) +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         panel.spacing.x=unit(2,'pt'),
         axis.text.x=element_blank(),
         #axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         axis.ticks.x=element_blank(),
         strip.text.x=element_text(angle=strip.text.angle, hjust=strip.text.hjust)
      )
   
   
   if(show.sample.names){
      p <- p + 
         theme(
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
            axis.ticks.x=element_line()
         )
   }
   
   return(p)
}

if(WRITE_OUTPUT){
   (function(){
      m <- sig_contribs$snv
      p1 <- plotSigContrib(m, sample.order=hclustCustom(m))
      
      m <- sig_contribs$indel
      p2 <- plotSigContrib(m, sample.order=hclustCustom(m), fill.palette='Pastel1')
      
      pdf(paste0(wd,'/plots/sig_contribs.pdf'), 11, 5)
      plot(p1)
      plot(p2)
      dev.off()
   })()
}


## Write supp table
if(WRITE_OUTPUT){
   (function(){
      l_contexts <- contexts
      l_sigs <- sig_contribs
      l_sigs$dbs <- NULL

      names(l_contexts) <- paste0('contexts.',names(l_contexts))
      names(l_sigs) <- paste0('sigs.',names(l_sigs))

      l <- c(l_contexts, l_sigs)
      l <- lapply(l, function(i){
         sample_names <- metadata$sample_name_2[ match(rownames(i), metadata$sample_name) ]
         data.frame(
            sample=sample_names, i,
            row.names=NULL, check.names=F
         )
      })
      openxlsx::write.xlsx(
         l, paste0(base_dir,'/analysis3/compare_mut_profiles/tables/S2_contexts_sigs_contribs.xlsx')
      )
   })()
}


# ## NMF --------------------------------
# nmfWrapper <- function(m, out.path, rank=2:15, ...){
#    if(!file.exists(out.path)){
#       require(NMF)
#       nmf_out <- nmf(m, rank=rank, seed=1, .options='v', ...)
#       saveRDS(nmf_out, out.path)
#    } else {
#       nmf_out <- readRDS(out.path)
#    }
#    return(nmf_out)
# }
# 
# calcProfileCosSim <- function(sigs_denovo, sigs_ref, contribs_denovo, context_counts){
# 
#    if(F){
#       sigs_denovo=nmf_out$fit$`9`@fit@H
#       sigs_ref=SBS_SIGNATURE_PROFILES_V3
#       contribs_denovo=nmf_out$fit$`9`@fit@W
#       context_counts=as.matrix(contexts$snv)
#    }
# 
#    sigs_denovo <- sigs_denovo / rowSums(sigs_denovo)
#    sigs_denovo[is.na(sigs_denovo)] <- 0
#    sigs_denovo <- t(sigs_denovo)
# 
#    ## Link de novo and reference sig names
#    cosSimMatrix <- function(profiles.de.novo, profiles.ref){
#       cosSim <- function(x, y) { sum(x*y)/sqrt(sum(x^2)*sum(y^2)) }
# 
#       out <- apply(profiles.de.novo, 2, function(i){
#          #i=profiles.de_novo[,1]
#          apply(profiles.ref, 2, function(j){
#             cosSim(i, j)
#          })
#       })
# 
#       t(out)
#    }
# 
#    cos_sim <- cosSimMatrix(sigs_denovo, sigs_ref)
#    #cos_sim[is.na(cos_sim)] <- 0
# 
#    max_sim <- round(apply(cos_sim,1,max),2) * 100
#    max_sim_ref_sig <- colnames(cos_sim)[ max.col(cos_sim) ]
#    sig_etiology <- sig_metadata$sig_etiology[ match(max_sim_ref_sig, sig_metadata$sig_name) ]
# 
#    sig_names <- paste0(
#       max_sim_ref_sig,'.',
#       sig_etiology,'.',
#       max_sim,'.',
#       'denovo_',1:ncol(sigs_denovo)
#    )
#    #sig_names <- gsub('^NA[.]','NOMATCH.',sig_names)
# 
#    colnames(sigs_denovo) <- sig_names
#    rownames(cos_sim) <- sig_names
#    cos_sim[is.na(cos_sim)] <- 0
# 
#    ## Annotate sig names to contribution matrix
#    #contribs_denovo <- i$W
#    colnames(contribs_denovo) <- sig_names
# 
#    ## Scale contribs
#    contexts_ss <- context_counts[rownames(contribs_denovo),]
#    mut_load_scale_factor <- rowSums(contexts_ss) / rowSums(contribs_denovo)
#    contribs_denovo <- contribs_denovo * mut_load_scale_factor
# 
#    ## Order sigs based on reference sigs
#    sig_order <- naturalsort::naturalorder(colnames(sigs_denovo))
#    sigs_denovo <- sigs_denovo[,sig_order]
#    cos_sim <- cos_sim[sig_order,]
#    contribs_denovo <- contribs_denovo[,sig_order]
# 
#    contribs_lsq <- mutSigExtractor::fitToSignatures(
#       context_counts[rownames(contribs_denovo),],
#       sigs_denovo,
#       use.r.implementation=F
#    )
# 
#   # plotContribRel(contribs_denovo)
# 
#    ## Output
#    list(
#       profiles=sigs_denovo,
#       cos_sim=cos_sim,
#       contribs_nmf=contribs_denovo,
#       contribs_lsq=contribs_lsq
#    )
# }
# 
# ## 
# nmf_out <- nmfWrapper(
#    as.matrix(contexts$snv),
#    out.path=paste0(wd,'/nmf/snv_allSamples.rds')
# )
# #plot(nmf_out)
# 
# profile_cos_sim <- calcProfileCosSim(
#    nmf_out$fit$`8`@fit@H, 
#    SBS_SIGNATURE_PROFILES_V3,
#    nmf_out$fit$`8`@fit@W,
#    as.matrix(contexts$snv)
# )
# 
# plotSigContrib(profile_cos_sim$contribs_nmf, show.etiology=F)
# 
# ##
# m <- as.matrix(contexts$snv)
# m <- m[subset(metadata, cohort=='Footprint', sample_name, drop=T),]
# 
# nmf_out <- nmfWrapper(
#    m,
#    out.path=paste0(wd,'/nmf/snv_footprintSamples.rds'),
#    rank=1:10
# )
# #plot(nmf_out)
# 
# profile_cos_sim <- calcProfileCosSim(
#    nmf_out$fit$`2`@fit@H, 
#    SBS_SIGNATURE_PROFILES_V3,
#    nmf_out$fit$`2`@fit@W,
#    as.matrix(contexts$snv)
# )
# 
# plotSigContrib(profile_cos_sim$contribs_nmf, show.etiology=F)

## Drivers ######################################################################################
## Oncoprint ================================
muts_som <- cacheAndReadData(paste0(base_dir,'/scripts/annotate_som_variants/mut_profile_som.txt.gz'))
muts_som <- muts_som[muts_som$sample %in% metadata$sample_name,]
muts_som$sample <- factor(muts_som$sample, metadata$sample_name)

gene_list <- openxlsx::read.xlsx(
   paste0(base_dir,'/analysis3/compare_mut_profiles/Compendium_Cancer_Genes_20200201.xlsx'),
   sheet='sel_genes'
)

prepDataOncoprint <- function(muts){
   
   #muts=subset(muts_som, ensembl_gene_id %in% gene_list$ensembl_gene_id)
   
   ## Remove irrelevant cols
   df <- muts
   df <- subset(df, select=-c(chrom, pos, ref, alt, hgvs_c))
   df$sample <- as.factor(df$sample)
   
   ## Select most impactful variant per sample per gene
   df <- df[order(df$max_score, decreasing=T),]
   df <- df[
      !duplicated(paste0(df$sample,':',df$ensembl_gene_id))
   ,]
   
   #df <- df[order(df$sample),]
   
   ## Rename eff
   df$snpeff_eff[df$snpeff_eff=='frameshift_variant'] <- 'frameshift'
   df$snpeff_eff[df$snpeff_eff=='missense_variant'] <- 'missense'
   df$snpeff_eff[df$snpeff_eff %in% c('stop_gained','stop_lost','start_lost')] <- 'nonsense'
   df$snpeff_eff[df$snpeff_eff %in% c('splice_acceptor_variant','splice_donor_variant')] <- 'essential_splice'
   df$snpeff_eff[df$snpeff_eff %in% c('disruptive_inframe_insertion','disruptive_inframe_deletion')] <- 'frameshift_inframe_disruptive'
   df$snpeff_eff[df$snpeff_eff %in% c('conservative_inframe_insertion','conservative_inframe_deletion')] <- 'frameshift_inframe_conservative'
   
   df$snpeff_eff[df$max_score<3] <- 'none'
   
   df$snpeff_eff <- factor(
      df$snpeff_eff,
      unique(c(
         'frameshift','frameshift_inframe_disruptive','frameshift_inframe_conservative','nonsense','missense','essential_splice',
         df$snpeff_eff[df$snpeff_eff!='none'],
         'none'
      ))
   )
   
   ## Order genes by mut freq
   df$snpeff_gene <- factor(df$snpeff_gene, unique(df$snpeff_gene))
   gene_mut_freq <- table(df$snpeff_gene[df$max_score>=3])
   gene_mut_freq <- sort(gene_mut_freq, decreasing=T)
   df$snpeff_gene <- factor(df$snpeff_gene, names(gene_mut_freq))
   
   df <- df[order(df$snpeff_gene),]
   df$sample <- factor(
      df$sample,
      unique(c(as.character(df$sample), levels(df$sample)))
   )
   
   ## Other cols
   df$supp_ann <- df$max_score_origin
   df$supp_ann[df$max_score<3 | df$supp_ann=='snpeff'] <- 'none'
   
   out <- df[,c('sample','snpeff_eff','snpeff_gene','max_score','supp_ann')]
   
   ## Add rows for missing samples
   missing_samples <- levels(df$sample)[ !(levels(df$sample) %in% df$sample) ]
   out_missing <- out[0,]
   out_missing <- out_missing[1:length(missing_samples),]; rownames(out_missing) <- NULL
   out_missing$sample <- missing_samples
   out_missing$snpeff_gene <- out$snpeff_gene[1]
   out_missing <- within(out_missing,{
      snpeff_eff <- 'none'
      max_score <- 0
      supp_ann <- 'none'
   })
   
   out <- rbind(out, out_missing)
   
   return(out)
}

##
p_oncoprint <- (function(){
   data <- prepDataOncoprint(
      subset(muts_som, ensembl_gene_id %in% gene_list$ensembl_gene_id)
   )
   data <- cbind(data, getMetadata(data$sample))
   
   ##
   m <- reshape2::dcast(
      data=subset(data, sample_group=='PCAWG_HCC'),
      formula=sample~snpeff_gene, value.var='snpeff_eff', fill='none'
   )
   rownames(m) <- m[,1]; m[,1] <- NULL
   m <- as.matrix(m)
   
   sum(apply(m,1,function(i){
      all(i=='none')
   }))
   
   data$exp_group <- factor(data$exp_group, DISEASED_LIVER_GROUPS)
   
   data <- data[order(data$snpeff_gene),]
   data$index_gene <- as.integer(data$snpeff_gene)
   data$snpeff_eff[data$snpeff_eff=='none'] <- NA
   
   data$index_sample <- as.integer(data$sample)

   grid_gene <- 1:(max(data$index_gene)-1) + 0.5
   grid_sample <- 1:(max(data$index_sample)-1) + 0.5
   
   ggplot(data, aes(x=index_gene, y=sample)) +
      facet_grid(exp_group~., scales='free',space='free') +
      
      #geom_hline(yintercept=grid_sample, color='darkgrey', size=0.05) +
      geom_vline(xintercept=grid_gene, color='darkgrey', size=0.05) +
      
      geom_tile(
         aes(fill=snpeff_eff), 
         color=ifelse(!is.na(data$snpeff_eff),'darkgrey',NA),
         size=0.05
      ) +
      
      scale_fill_brewer(name='Variant type', palette='Set3', na.translate=F) +
      
      scale_x_continuous(
         name='Gene name',
         breaks=unique(data$index_gene), labels=unique(data$snpeff_gene),
         expand=c(0,0),
         sec.axis = dup_axis()
      ) +
      
      scale_y_discrete(name='Sample', limits=rev) +
      
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         strip.text.y=element_text(angle=0, hjust=0),
         axis.text.x.top=element_text(angle=45, hjust=0),
         axis.text.x.bottom=element_text(angle=45, hjust=1),
         axis.title.x=element_blank()
      )
})()

if(WRITE_OUTPUT){
   pdf(paste0(wd,'/plots/oncoprint.pdf'), 10, 9)
   plot(p_oncoprint)
   dev.off()
}

## dnds ================================
if(WRITE_OUTPUT){
   library(dndscv)
   muts_som_all_raw <- cacheAndReadData(paste0(base_dir,'/analysis3/compare_mut_profiles/muts_som_all.txt.gz'))
   #muts_som_all$exp_group <- getMetadata(muts_som_all$sample, 'exp_group')
   
   ## Subset and get metadata
   muts_som_all <- subset(muts_som_all_raw, sample %in% metadata$sample_name)
   muts_som_all <- cbind(
      getMetadata(muts_som_all$sample,c('cohort','exp_group','sample_group')),
      muts_som_all
   )
   muts_som_all <- subset(muts_som_all, exp_group %in% c('Healthy','ALC','NASH','PSC','PCAWG_HCC','PCAWG_CCA'))
   
   ## Flatten clones from same donors
   muts_som_all <- within(muts_som_all,{
      sample_group <- as.character(sample_group)
      sample_group[cohort=='PCAWG'] <- sample[cohort=='PCAWG']
   })
   
   muts_som_all <- muts_som_all[
      !duplicated( muts_som_all[,c('sample_group','chrom','pos','ref','alt')] )
      ,]
   muts_som_all$sample <- NULL
   
   ## Only keep SNVs and indels
   muts_som_all$ref_len <- nchar(muts_som_all$ref)
   muts_som_all$alt_len <- nchar(muts_som_all$alt)
   
   muts_som_all <- muts_som_all[
      with(
         muts_som_all, 
         (ref_len==1 & alt_len==1) | (ref_len==1 & alt_len>1) | (ref_len>1 & alt_len==1)
      )
      ,]
   
   ## dndscv per exp group
   exp_group_names <- levels(droplevels(muts_som_all$exp_group))
   dnds_out <- lapply(exp_group_names, function(i){
      message('## Performing dndscv for: ',i)
      dndscv(
         subset(
            muts_som_all, 
            exp_group==i, 
            select=c(sample_group, chrom, pos, ref, alt)
         )
      )
   })
   names(dnds_out) <- exp_group_names
   
   ## Gather the gene enrichment table
   dnds_out.sel_cv <- lapply(exp_group_names, function(i){
      cbind(
         exp_group=i,
         dnds_out[[i]]$sel_cv
      )
   })
   dnds_out.sel_cv <- do.call(rbind, dnds_out.sel_cv)
   rownames(dnds_out.sel_cv) <- NULL
   
   ## Export output
   write.table(
      dnds_out.sel_cv, 
      paste0(base_dir,'/analysis3/compare_mut_profiles/dnds_out.sel_cv.txt'),
      sep='\t',quote=F,row.names=F
   )

   #dnds_out.sel_cv <- read.delim(paste0(base_dir,'/analysis3/compare_mut_profiles/dnds_out.sel_cv.txt.gz'))
   dnds_out.sel_cv <- dnds_out.sel_cv[order(dnds_out.sel_cv$qglobal_cv),]
   openxlsx::write.xlsx(
      dnds_out.sel_cv[,c('exp_group','gene_name','n_syn','n_mis','n_non','n_spl','n_ind','pglobal_cv','qglobal_cv')],
      paste0(base_dir,'/analysis3/compare_mut_profiles/tables/S3_dnds_out.sel_cv.compact.xlsx')
   )
}

## Ploidies ########################################################################################
## Load data ================================
ploidy <- read.delim(paste0(base_dir,'/scripts/extract_aneuploidy/m_ploidy.txt.gz'), check.names=F)
ploidy <- as.matrix(ploidy)

DISEASED_LIVER_GROUPS_2 <- c("Healthy","ALC","NASH","PSC","PCAWG_HCC","PCAWG_CCA")
metadata2 <- subset(
   metadata_raw, 
   !(sample_name %in% SAMPLE_BLACKLIST) & exp_group %in% DISEASED_LIVER_GROUPS_2,
   -c(germ_vcf,som_vcf,sv_vcf,cnv_tsv)
)

ploidy <- ploidy[as.character(metadata2$sample_name),]

## Plot ================================
p_ploidy <- (function(){
   ploidy_melt <- melt(ploidy)
   colnames(ploidy_melt) <- c('sample','chrom_arm','value')
   
   ploidy_melt <- cbind(ploidy_melt, getMetadata(ploidy_melt$sample, df=metadata2))
   ploidy_melt$exp_group <- factor(ploidy_melt$exp_group, DISEASED_LIVER_GROUPS_2)
   ploidy_melt$chrom_arm <- factor(ploidy_melt$chrom_arm, unique(ploidy_melt$chrom_arm ))
   
   hc <- hclust(dist(ploidy))
   ploidy_melt$sample <- factor(ploidy_melt$sample, rownames(ploidy)[hc$order])
   
   ploidy_melt$value[ploidy_melt$value>4] <- '>=5'
   ploidy_values <- c(1:4, '>=5')
   ploidy_melt$value <- factor(ploidy_melt$value, ploidy_values)
   
   ploidy_colors <- c('#9CC4DB','white','#FAE09A','#EE9679','#9373A8')
   names(ploidy_colors) <- ploidy_values
   
   ggplot(ploidy_melt, aes(x=chrom_arm, y=sample)) +
      facet_grid(exp_group~., space='free_y', scales='free_y') +
      geom_tile(aes(fill=value), color='grey') +
      
      scale_fill_manual(values=ploidy_colors, name='Ploidy', na.translate=F) +
      scale_y_discrete(name='Samples') +
      scale_x_discrete(name='Region', expand=c(0,0)) +
      
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         panel.spacing.y=unit(2,'points'),
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         strip.text.y=element_text(angle=0, hjust=0)
      )
})()

pdf(paste0(wd,'/plots/ploidy.pdf'), 11, 9)
plot(p_ploidy)
dev.off()



