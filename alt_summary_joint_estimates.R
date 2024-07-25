# an estimator of joint SNP associations. 'Joint' here means conditional on all SNPs you give this function
alt_summary_joint_df=function(df,Rhat,n=1000) {
  # df: dataframe containing columns 'betahat' and 'betahatse' (rows are SNPs); respectively marginal association estimates and their SEs
  # Rhat [m x m matrix]: estimated LD matrix for SNPs in betahat and betahatse
  # n [scalar]: number of people used to estimate LD matrix (ie size of reference panel)
  betahat = df$betahat
  betahatse = df$betahatse
  m=length(betahat)
  K=fastmatrix::commutation(m=m,n=m,matrix=TRUE)
  Im=diag(m)
  Im2=diag(m^2)
  Im2K=Im2+K
  Thetahat=solve(Rhat)
  dd=median(betahatse^2)
  D=diag(dd)
  bhat=Thetahat%*%betahat
  adj=(n-m-1)/n
  bhat=adj*Thetahat%*%betahat
  var1=dd*adj*Thetahat
  P=kronecker(t(betahat),Im)
  KTheta=adj^2*kronecker(Thetahat,Thetahat)
  d1=n^2/((n-m)*(n-m-1)^2*(n-m-3))
  d2=d1*(n-m-1)
  vv=c(adj*Thetahat)
  vv=d1*vv%*%t(vv)
  TT=d2*Im2K%*%KTheta
  var2=adj^2*P%*%(vv+TT)%*%t(P)
  SEs=sqrt(diag(var1+var2))
  Ps=pchisq(bhat^2/SEs^2,1,lower.tail=FALSE)
  list(JointEstimates=bhat,JointEstimateSEs=SEs,JointEstimateP=Ps)
}

library(dplyr)
library(data.table)
library(ggplot2)
setwd("/Users/houy2/Documents/projects/AD/AD_genetics/results/update_v9_Jul072023/Noah")
#load LD matrix
df_LD<-read.delim("./r2_19485_8rsid.txt", sep = "\t", header = T)
Rhat<-as.matrix(df_LD[2:5]) # matrix of squared LD coefficients

df_ephx2<-read.delim("./EPHX2_multi_GWAS_all.tsv", sep = "\t", header = T)
df_ephx2 <- df_ephx2 %>%
  rename(betahat = beta, betahatse = se)

# Get unique gwas
unique_gwas <- unique(df_ephx2$outcome)
df_new <- data.table()

for (G in unique_gwas) {
  subset_df <- subset(df_ephx2, outcome == G)
  
  df_ephx2_reordered <- df_LD %>%
    select(RS_number) %>% # Select only the RS_number column or the column that is used for ordering
    inner_join(subset_df, by = "RS_number")
  
  # Call the modified function with the dataframe as input
  results_df <- alt_summary_joint_df(df_ephx2_reordered,sqrt(Rhat)) # can use sqrt(Rhat) bc all LD is positive
  results_df
  
  #save to df
  se_df <- data.frame(SE = results_df$JointEstimateSEs)
  rownames(se_df) <- rownames(results_df$JointEstimates)
  final_df <- data.frame(
    JointEstimateBeta = results_df$JointEstimates[,1],
    JointEstimateSE = se_df$SE,
    JointEstimateP = results_df$JointEstimateP[,1]
  )
  final_df$GWAS <- G
  final_df$rsid <- rownames(final_df)
  df_new <- rbind(df_new, final_df)
  
}

write.table(df_new,"JointEstimate_EPHX2_eur.tsv",sep = "\t",row.names = F)

df_eur<-subset(df_new, subset = (GWAS != 'AD_Kunkle_AFA_M1_b37_2020N' & GWAS != 'AD_Kunkle_AFA_M2_b37_2020'))

df_eur$GWAS<- factor(df_eur$GWAS, levels = c('AD_Bellenguez_2022N',
                                             'AD_NMAFallNoUKB23',
                                             'AD_GCST90012877_b37_2021N',
                                             'AD_GCST90012878_b37_2021N',
                                             'AD_Jansen_2019N',
                                             'AD_Kunkle_Stage1_2019N',
                                             'AD_Kunkle_Stage2_2019N'))

line_color<-c('AD_Bellenguez_2022N'='#90a4ae',
              'AD_NMAFallNoUKB23'='#795548',
              'AD_GCST90012877_b37_2021N'='#1565c0',
              'AD_GCST90012878_b37_2021N'='#f44336',
              'AD_Jansen_2019N'='#ff9900',
              'AD_Kunkle_Stage1_2019N'='#8d24aa',
              'AD_Kunkle_Stage2_2019N'='#B39DDB')


df_eur$rsid<- factor(df_eur$rsid, levels = c('rs11787077',
                                             'rs4236673',
                                             'rs9331896',
                                             'rs2279590',
                                             'rs1532278',
                                             'rs73223431',
                                             'rs751141',
                                             'rs2741342'))

dp_C <- ggplot(df_eur, aes(x = rsid, y =-log10(JointEstimateP), fill = group)) +
  #geom_boxplot() +
  geom_point(aes(fill = GWAS), size = 5,shape = 21,alpha=0.6, position = position_dodge(0.2)) +
  geom_hline(yintercept=-log10(5e-8), linetype=2,col="black",lwd=0.8)+
  geom_hline(yintercept=-log10(1e-5), linetype=2,col="gray",lwd=0.8)+
  scale_color_manual(values = line_color) +
  scale_fill_manual(values = line_color) +
  labs(title="Joint Estimate", x="", y="-Log10(P value)")+
  theme_classic() +
  theme(
    panel.border = element_blank(),
    legend.position = "right",
    axis.line = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

dp_C


pdf('./Joint Estimate_eur.pdf',width=7, height=4)
print(dp_C) # need print plot first
dev.off()

png('./Joint Estimate_eur.png',width=16, height=12,units="cm",res=500)
print(dp_C)
dev.off()


