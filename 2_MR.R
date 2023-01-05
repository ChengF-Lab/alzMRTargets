start<-Sys.time()
library(TwoSampleMR)
library(dplyr)
library(stringr)
#set the working directory from which the files will be read from
QTL_type<-"pqtl"
result_dir<-"AD_MR_pQTL_Jan042023"
GWAS_dir<-"AD"
File_PATTERN<-".tsv"

PATH="/hdd/yuan/project/AD"
input=sprintf("MR/QTL/xQTL_F_stat/%s",QTL_type)
input_GWAS=sprintf("GWAS/%s",GWAS_dir)
output_IV=sprintf("MR/results/%s/IV",result_dir)
output_MR=sprintf("MR/results/%s/MR_%s",result_dir,QTL_type)
WM_times=10000 # iterations in weighted median

setwd(file.path(PATH))
PATH_I <- file.path(PATH,input)
PATH_G <- file.path(PATH,input_GWAS)
PATH_O_IV <- file.path(PATH,output_IV)
PATH_O_MR <- file.path(PATH,output_MR)

file_list <- list.files(path=PATH_I,pattern=File_PATTERN)
GWAS_list <- list.files(path=PATH_G,pattern=File_PATTERN)

for (f in file_list){
  print(f)
  input_file<-sprintf('./%s/%s',input, f)
  LD_xQTL<-read.delim(input_file,header = T, sep="\t")
  
  exposure_dat <-format_data(
    LD_xQTL,
    type = "exposure",
    snps = NULL,
    header = TRUE,
    phenotype_col = "symbol",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "se",
    #eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    gene_col = "symbol",
    chr_col = "chr",
    pos_col = "pos"
  )
  
  #LD clumping
  library(ieugwasr)
  library(genetics.binaRies)
  LD_clump<-ld_clump(tibble(rsid=exposure_dat$SNP, pval=exposure_dat$pval.exposure),
                     clump_r2=0.2,
                     clump_kb=500,
                     plink_bin = get_plink_binary(),
                     bfile = "/hdd/yuan/project/AD/GWAS/1kg.v3/EUR")
  LD_exposure_dat <-exposure_dat %>% filter(SNP %in% LD_clump$rsid)
  
  for (g in GWAS_list){
    print(g)
    GWAS_file<-sprintf('./%s/%s',input_GWAS, g)
    GWAS_name<-strsplit(g, ".",fixed=T)[[1]][1]
    outcome_dat <- read_outcome_data(
      snps = LD_exposure_dat$SNP,
      filename = GWAS_file,
      sep = "\t",
      phenotype_col = g,
      snp_col = "rsid",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "eaf",
      pval_col = "pval"
    )
    #harmonise the exposure and outcome data
    dat <- harmonise_data(LD_exposure_dat, outcome_dat)
    dat <- subset(dat, dat$mr_keep == TRUE)
    #save IVs file
    write.table(dat, file = file.path(PATH_O_IV,str_c('IV-500kb#',GWAS_name,'#',f)), sep="\t", row.names = F)
    
    #======run MR========
    #Wald ratio for all SNPs
    mr_single <- mr_singlesnp(dat)
    mr_single<-mr_single[mr_single$SNP %in% dat$SNP, ]
    if (plyr::empty(mr_single) == FALSE) {
      or_single<-generate_odds_ratios(mr_single)
      or_single$FDR<- p.adjust(or_single$p, method = "fdr")
      write.table(or_single, file = file.path(PATH_O_MR,str_c('MR-500kb_single#',GWAS_name,'#',f)), sep="\t", row.names = F)
    }
    
    #MR for each exposure
    rsid_table <- as.data.frame(table(dat$exposure))  #exposure table
    #1--Wald ratio, No. snps of exposure==1
    df_rsid_1 <- subset(dat, dat$exposure %in% rsid_table$Var1[which(rsid_table$Freq ==1)])
    mr_rsid_1 <- mr(df_rsid_1)
    if (plyr::empty(mr_rsid_1) == FALSE) {
      or_snp1<-generate_odds_ratios(mr_rsid_1)
      or_snp1$FDR<- p.adjust(or_snp1$p, method = "BH")
      write.table(or_snp1, file = file.path(PATH_O_MR,str_c('MR-500kb_snp1#',GWAS_name,'#',f)), sep="\t", row.names = F)
    }
    #2--No. snps of exposure==2
    df_rsid_2 <- subset(dat, dat$exposure %in% rsid_table$Var1[which(rsid_table$Freq ==2)])
    mr_rsid_2 <- mr(df_rsid_2)
    if (plyr::empty(mr_rsid_2) == FALSE) {
      or_snp2<-generate_odds_ratios(mr_rsid_2)
      or_snp2$FDR<- p.adjust(or_snp2$p, method = "BH")
      write.table(or_snp2, file = file.path(PATH_O_MR,str_c('MR-500kb_snp2#',GWAS_name,'#',f)), sep="\t", row.names = F)
    }
    
    #3--No. snps of exposure>2
    df_rsid_more <- subset(dat, dat$exposure %in% rsid_table$Var1[which(rsid_table$Freq >2)])
    if (plyr::empty(df_rsid_more) == TRUE) {
      print("no genes in this script have > 2 snps available")
    } else {
      IVW_results <- data.frame()
      IVW_robust_results<- data.frame()
      
      Weighted_Median_results <- data.frame()
      
      maxlik_results <- data.frame()
      
      MRpresso_results<-data.frame()
      
      unique_exposures0 <- unique(df_rsid_more$exposure)
      unique_exposures <- unique_exposures0[!is.na(sort(unique_exposures0))]
      
      for (G in unique_exposures) {
        df_gene <- subset(df_rsid_more, exposure == G) #subset data to keep only 1 exposure of interest
        print(G)
        #non correlated
        dat_gene <- dat_to_MRInput(df_gene)
        
        IVW<- MendelianRandomization::mr_ivw(dat_gene[[1]],penalized = FALSE)
        if (IVW@Heter.Stat[2] < 0.05){
          IVW_robust<- MendelianRandomization::mr_ivw(dat_gene[[1]],penalized = TRUE)
          try(IVW_robust_results <- rbind(IVW_robust_results, c(IVW_robust@Exposure, IVW_robust@Outcome, IVW_robust@SNPs, IVW_robust@class[1], IVW_robust@Estimate, IVW_robust@StdError, IVW_robust@CILower,IVW_robust@CIUpper,IVW_robust@Pvalue)))
        }else{
          try(IVW_results <- rbind(IVW_results, c(IVW@Exposure, IVW@Outcome, IVW@SNPs, IVW@class[1], IVW@Estimate, IVW@StdError, IVW@CILower,IVW@CIUpper,IVW@Pvalue,IVW@Heter.Stat[1],IVW@Heter.Stat[2])))
        }
        
        Weighted_Median<-MendelianRandomization::mr_median(dat_gene[[1]],weighting = "weighted",iterations = WM_times,seed = 666)
        try(Weighted_Median_results <- rbind(Weighted_Median_results, c(Weighted_Median@Exposure, Weighted_Median@Outcome, Weighted_Median@SNPs, Weighted_Median@class[1], Weighted_Median@Estimate, Weighted_Median@StdError, Weighted_Median@CILower,Weighted_Median@CIUpper,Weighted_Median@Pvalue)))
        
        maxlik <- MendelianRandomization::mr_maxlik(dat_gene[[1]])
        try(maxlik_results <- rbind(maxlik_results, c(maxlik@Exposure, maxlik@Outcome, maxlik@SNPs, maxlik@class[1], maxlik@Estimate, maxlik@StdError,maxlik@CILower,maxlik@CIUpper,maxlik@Pvalue,maxlik@Heter.Stat[1],maxlik@Heter.Stat[2])))
        
        if (nrow(df_gene)>1000) {
          mr_presso<-run_mr_presso(df_gene, NbDistribution = 5000, SignifThreshold = 0.05)
          df_mrpresso<-mr_presso[[1]]$"Main MR results"
          df_mrpresso$outcome<-g
          df_mrpresso$Exposure<-G
          df_mrpresso$xqtl<-f
          MRpresso_results<-rbind(MRpresso_results,df_mrpresso)
        } else if (nrow(df_gene)>3) {
          mr_presso<-run_mr_presso(df_gene, NbDistribution = 1000, SignifThreshold = 0.05)
          df_mrpresso<-mr_presso[[1]]$"Main MR results"
          df_mrpresso$outcome<-g
          df_mrpresso$Exposure<-G
          df_mrpresso$xqtl<-f
          MRpresso_results<-rbind(MRpresso_results,df_mrpresso)
          
        }
      }
      #===========results===========
      #IVW
      names(IVW_results)<-c('Exposure','Outcome','SNPs','method','Estimate','StdError','CILower','CIUpper','Pvalue','Heter.Stat[1]','Heter.p')
      IVW_results$method_type <- 'IVW'
      if (plyr::empty(IVW_robust_results) == TRUE) {
        IVW_results$FDR<- p.adjust(IVW_results$Pvalue, method = "BH")
        write.table(IVW_results, file = file.path(PATH_O_MR,str_c('MR-more-500kb_IVW#',GWAS_name,'#',f)), sep="\t", row.names = F)
      }else{
        names(IVW_robust_results)<-c('Exposure','Outcome','SNPs','method','Estimate','StdError','CILower','CIUpper','Pvalue')
        IVW_robust_results$method_type <- 'IVW_robust'
        IVW_all<-bind_rows(IVW_results,IVW_robust_results)
        IVW_all$FDR<- p.adjust(IVW_all$Pvalue, method = "BH")
        write.table(IVW_all, file = file.path(PATH_O_MR,str_c('MR-more-500kb_IVW#',GWAS_name,'#',f)), sep="\t", row.names = F)
      }
      
      #Weighted_Median
      names(Weighted_Median_results)<-c('Exposure','Outcome','SNPs','method','Estimate','StdError','CILower','CIUpper','Pvalue')
      Weighted_Median_results$FDR<- p.adjust(Weighted_Median_results$Pvalue, method = "BH")
      Weighted_Median_results$method_type <- 'Weighted_Median'
      Weighted_Median_results$iterations<-WM_times
      write.table(Weighted_Median_results, file = file.path(PATH_O_MR,str_c('MR-more-500kb_WM#',GWAS_name,'#',f)), sep="\t", row.names = F)
      
      #maxlik
      names(maxlik_results)<-c('Exposure','Outcome','SNPs','method','Estimate','StdError','CILower','CIUpper','Pvalue','Heter.Stat[1]','Heter.p')
      maxlik_results$FDR<- p.adjust(maxlik_results$Pvalue, method = "BH")
      maxlik_results$method_type <- 'maxlik'
      
      write.table(maxlik_results, file = file.path(PATH_O_MR,str_c('MR-more-500kb_Maxlik#',GWAS_name,'#',f)), sep="\t", row.names = F)
      
      
      #MRpresso
      if (plyr::empty(MRpresso_results) == F){
        write.table(MRpresso_results, file = file.path(PATH_O_MR,str_c('MR-more-500kb_MRpresso#',GWAS_name,'#',f)), sep="\t", row.names = F)
      }
    }
  }
  
}

print('=========')
print("end")
end<-Sys.time()
print(str_c("start:",start,"-->end:",end))
