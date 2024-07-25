rm(list=ls(all=TRUE))
source('functions.R')
library(dplyr);library(mvnfast);library(ggpubr)
nff=read.csv('GWAS_overlapping_matrix.csv')[,-1] %>% as.matrix()
pff=cov2cor(nff) # proportion of overlapping subjects in AD GWAS
nAD=diag(nff) # AD GWAS sample sizes
N0=((1/sqrt(nAD))%*%t(1/sqrt(nAD)))*pff
q=length(nAD) # number of AD GWAS
p=12 # number of gene phenotype trait-cohort pairs
ngenes=10
# gene rank defined by theta and power to reject H0: theta=0
thetas=seq(0,1,length.out=ngenes)
h2s=seq(0.01,0.15,length.out=ngenes)
trueranks=seq(ngenes,1,-1)
COMP=outer(thetas,h2s)
nX=500
RhoBB=ar1(p,0.95) # genetic correlation between gene phenotypes
M=3 # number of causal cis SNPs for each gene phenotype
q0=qchisq(1-5e-5,1) # IV signficance threshold
niter=1000
zSCORES=pSCORES=matrix(nr=niter,nc=ngenes)
for(thisgene in 1:ngenes) {
  THETA=matrix(thetas[thisgene]/p,nr=p,nc=q)
  h2i=h2s[thisgene]
  SigmaBB=h2i*RhoBB/M
  for(iter in 1:niter) {
    B=rmvn(M,rep(0,p),SigmaBB)
    adj=h2i/colSums(B^2)
    B=B%*%diag(sqrt(adj))
    B=B+matrix(rnorm(M*p,0,1/sqrt(nX)),M,p)
    PJ=sapply(1:p,function(h) pchisq(q0,1,ncp=nX*B[,h]^2,lower.tail=FALSE))
    DJ=sapply(1:p,function(h) rbinom(M,1,PJ[,h]))
    A=B%*%THETA # rows are SNPs, columns are AD GWAS sumstats
    for(h in 1:q) A=A+rmvn(M,rep(0,q),N0)
    # perform IV selection
    THETAHAT=THETAHATSE=matrix(NA,nr=p,nc=q)
    for(i in 1:p) {
      ix=which(DJ[,i]==1)
      if(length(ix)==0) next
      bix=(B[ix,i])
      aix=as.matrix(A[ix,])
      aix=matrix(aix,nr=length(ix),nc=q)
      for(j in 1:q) {
        THETAHAT[i,j]=solve(t(bix)%*%bix)%*%t(bix)%*%aix[,j]
        if(length(ix)==1) {
          se=sqrt(1/nAD[j]/bix^2+aix[,j]^2/nX/bix^4)
        } else {
          rss=sum((aix[,j]-bix*c(THETAHAT[i,j]))^2)/length(ix)
          se=sqrt(rss*diag(solve(t(bix)%*%bix)))
        }
        THETAHATSE[i,j]=se
      }
    }
    # calculate score for this gene
    Zs=THETAHAT/THETAHATSE
    zSCORES[iter,thisgene]=sum(Zs,na.rm=T)/ngenes
    pSCORES[iter,thisgene]=sum(Zs^2>1.96^2,na.rm=T)
  }
  cat('gene ',thisgene,' done\n',sep='')
}
# check that var(zSCORES[,1])\approx p*q
# calculate ranks for each iteration
zranks=t(apply(zSCORES,1,function(h) rev(rank(h,ties.method='random'))))
pranks=t(apply(pSCORES,1,function(h) rev(rank(h,ties.method='random'))))
zRANKS=data.frame(guess=c(zranks),gene=rep(1:ngenes,each=niter),truerank=rep(trueranks,each=niter))
pRANKS=data.frame(guess=c(pranks),gene=rep(1:ngenes,each=niter),truerank=rep(trueranks,each=niter))
bothranks=bind_rows(
  zRANKS %>% mutate(type='using Z-statistics'),
  pRANKS %>% mutate(type='using significant P-value counts')
)
# plot of estimated ranks using Z-scores
p1=bothranks %>%
  filter(type=='using Z-statistics') %>%
  ggplot(aes(x=factor(truerank),fill=factor(guess))) +
  geom_bar(position='fill') +
  labs(x='true rank of gene',y='proportion of simulations',title='a) Performance of gene ranks using alzScores') +
  scale_fill_brewer('estimated rank of gene',palette='Spectral') +
  theme_classic() +
  # facet_wrap(~type) +
  theme(legend.position='bottom')
# plot distribution of spearman rank correlations
zcorr=apply(zranks,1,function(h) cor(trueranks,h,method='spearman'))
pcorr=apply(pranks,1,function(h) cor(trueranks,h,method='spearman'))
bdf=data.frame(scor=c(zcorr,pcorr),type=rep(c('using Z-statistics','using significant P-value counts'),each=niter))
p2=bdf %>%
  filter(type=='using Z-statistics') %>%
  ggplot(aes(x=type,y=scor)) +
  geom_violin(fill='gray90') +
  geom_boxplot(fill='white',width=1/4) +
  theme_classic() +
  lims(y=c(0,1)) +
  labs(title='b) Correlation between estimated and true ranks',y=expression('Spearman '*rho),x='') +
  theme(axis.text.x=element_blank())
# plot proportion of correctly guessed ranks
propcorrect=function(guess,truth) sum(guess==truth)/length(guess)
propcorrectse=function(guess,truth) {p=propcorrect(guess,truth);p*(1-p)/sqrt(length(guess))}
zp=zRANKS %>% group_by(truerank) %>% 
  summarise(x=propcorrect(guess,truerank[1]),xse=propcorrectse(guess,truerank[1]))
pp=pRANKS %>% group_by(truerank) %>% 
  summarise(x=propcorrect(guess,truerank[1]),xse=propcorrectse(guess,truerank[1]))
both=bind_rows(
  zp %>% mutate(type='using Z-statistics'),
  pp %>% mutate(type='using significant P-value counts')
)
p3=both %>%
  filter(type=='using Z-statistics') %>%
  ggplot(aes(factor(truerank),x)) +
  geom_errorbar(aes(ymin=x-2*xse,ymax=x+2*xse),width=0.2) +
  geom_point(size=2/3) +
  # facet_wrap(~type) +
  theme_classic() +
  labs(x='true rank of gene',y='proportion of correctly estimated ranks',
       title='c) Distribution of correctly estimated gene ranks') +
  scale_y_continuous(breaks=seq(0,1,0.25),labels=seq(0,1,0.25),limits=c(0,1))
# put all together
ggarrange(p1,p2,p3,nrow=1,ncol=3)



