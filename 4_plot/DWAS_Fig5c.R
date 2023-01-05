setwd("C:/Users/houy2/Documents/Yuan/projects/AD/AD_genetics/results/update_v6_Sep132022/DWAS/update_July18/image")

library(ggplot2)
library(ggpubr)
library(gridExtra)

RR_data<-read.delim("MCI_DWAS_24_forestplot.tsv",header=T,check.names=F)
summary(RR_data)

RR_data$method<- factor(RR_data$method, 
                        levels = c('EXP','Abs','Rel','MPR'))
RR_data$Drug<- factor(RR_data$Drug, 
                      levels = c('Apixaban',
                                 'Levetiracetam',
                                 'Baclofen',
                                 'Rosuvastatin',
                                 'Spironolactone',
                                 'Torsemide',
                                 'Lamotrigine',
                                 'Propranolol',
                                 'Olmesartan',
                                 'Oxycodone',
                                 'Latanoprost',
                                 'Gabapentin',
                                 'Fluticasone',
                                 'Mirabegron',
                                 'Rivaroxaban',
                                 'Pantoprazole',
                                 'Losartan',
                                 'Amlodipine',
                                 'Atorvastatin',
                                 'Bupropion',
                                 'Oxybutynin',
                                 'Trazodone',
                                 'Clopidogrel'
                      ))
#"#BAA286","#4FC1E9","#7DB1B1",'#AC92EC',
p <- ggplot(data=RR_data,
            aes(x = Drug,y = HR_MCI,color=method))+
  geom_hline(yintercept =1, linetype=2)+
  geom_errorbar(aes(ymin=Lower_CI, ymax=Upper_CI),width=0,cex=1,alpha =1)+ 
  #geom_point(alpha = 1,aes(size = log10P))+
  geom_point(alpha = 1,shape=15,size=6)+
  xlab('')+ ylab("Hazard ratio, 95% CI")+
  scale_color_manual(values=c('#e57373',"#7986CB","#4DB6AC","#FFB74D"))+
  #scale_color_gradientn(limits = c(-1,1),colors = c("#c51b7d","#f1b6da","#f7f7f7", "#b8e186", "#4d9221")) +
  facet_wrap(~method,strip.position="right",nrow=40,scales = "free_y") +
  #facet_grid(. ~ Drug)+
  #theme_bw()+
  theme_classic()+
  theme(plot.title=element_text(size=16,face="bold"),
        panel.grid.major.y = element_line(size = 0.3, colour = "grey80",linetype =1),
        legend.position="left",
        axis.text.y=element_text(size=16,face="bold"),
        axis.text.x=element_text(hjust=1,vjust = 1,size=16,angle=45,colour = "Black"),
        axis.title=element_text(size=12),
        strip.text.y = element_text(hjust=0,vjust = 1,size=16,angle=0,colour = "Black"),
        strip.text.x = element_text(hjust=0,vjust = 1,size=16,angle=0,colour = "Black"))
p
pdf("forest_24drugs_v3.pdf",width=12, height=6)
print(p) # need print plot first
dev.off()

png("forest_24drugs_v3.png",width=15, height=12,units="cm",res=300)
print(p)
dev.off()

