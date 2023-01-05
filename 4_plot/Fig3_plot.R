library(ggplot2)
library(cowplot)

setwd("PATH")

df<-read.delim("EUR_MR_target_IVW_hp.tsv",header=T,check.names=F)
summary(df)
df$symbol <- factor(df$symbol, levels = rev(c('EPHX2','PRKCB','FOLH1','BACE2','EGFR','PLG','CD38','PARP10','MARK4','DRD2','CHRNE','DAGLB','APH1B','VKORC1','APP','PSEN2','PDE2A','NTRK1','GRIN3A','CACNA1D','F10','PDE5A','TAS2R60','GABBR1','NR3C1','ADRA2A','FPR1')))

df$xQTL_type <- factor(df$xQTL_type, levels = c(
  'Cortex-ROSMAP.All_pQTL',
  'Cortex-ROSMAP.Ctl_pQTL',
  'Cortex-Banner_pQTL',
  'Cortex-ROSMAP_eQTL',
  'Cortex-Mayo_eQTL',
  'Cortex-Metabrain_eQTL',
  'Cortex-Meta_eQTL',
  'Cerebellum-Mayo_eQTL',
  'Cerebellum-Metabrain_eQTL',
  'BasalGanglia-Metabrain_eQTL',
  'Hippocampus-Metabrain_eQTL',
  'SpinalCord-Metabrain_eQTL'
))
df$Outcome <- factor(df$Outcome, levels = c('AD_Bellenguez_2022','LOAD-Wightman_2021','AD-Schwartzentruber_2021','Family.Htx.AD-Schwartzentruber_2021','LOAD-Jansenetal_2019','AFA-Kunkle.M1_2020','AFA-Kunkle.M2_2020'))

df$color <- factor(df$color, levels = c('>1','0.1','0.01','0.001','0','-0.001','-0.01','-0.1','<-1'))

pa<-ggplot(df, aes(x = xQTL_type, y = symbol, fill = color)) +
  geom_tile(color = "white",lwd = 0.5,linetype = 1)+
  xlab('')+ ylab("")+
  scale_fill_manual(values=c("#c51b7d","#de77ae","#f1b6da","#fde0ef","#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221"),na.value = "#FAFAFA")+
  scale_y_discrete(position = "right")+
  scale_x_discrete(position = "top")+
  facet_wrap(~Outcome,strip.position="top",nrow=1) +
  theme(legend.text = element_text(colour="black",size=12),
        strip.background = element_rect(fill=c("#FAFAFA")),
        strip.placement = "outside",
        axis.line = element_line(size=0.3, colour = "white"),
        axis.ticks.x = element_line(colour = "black", size = 0.3),
        axis.ticks.y = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 12,angle=70,hjust=0,vjust=0),
        axis.text.y=element_text(colour="black", size = 12))+
  theme(panel.spacing.x = unit(1, "mm"))
pa
pp<- pa +
  geom_text(aes(label=sig),size=3.5,na.rm=TRUE,vjust =0.8) # you can play with the size
pp


pdf("./image/MR_target_IVW_hp_H.pdf",width=12, height=7)
print(pp)
dev.off()

df_target<-read.delim("EUR_MR_target_hist.tsv",header=T,check.names=F)
summary(df_target)
df_target$symbol <- factor(df_target$symbol, levels = rev(c('EPHX2','PRKCB','FOLH1','BACE2','EGFR','PLG','CD38','PARP10','MARK4','DRD2','CHRNE','DAGLB','APH1B','VKORC1','APP','PSEN2','PDE2A','NTRK1','GRIN3A','CACNA1D','F10','PDE5A','TAS2R60','GABBR1','NR3C1','ADRA2A','FPR1')))

df_target$group <- factor(df_target$group, levels = c('sig_posi#_p',
                                                      'sig_posi#_fdr',
                                                      'sig_neg#_p',
                                                      'sig_neg#_fdr'))


a1 <- ggplot(df_target, aes(x = count, y = symbol)) +
  geom_bar(stat="identity",aes(fill =group)) +
  ylab("")+xlab("")+
  scale_x_continuous(position = 'top')+
  scale_y_discrete(position = "left")+
  geom_vline(xintercept =0)+
  scale_fill_manual(values=c("#f1b6da","#c51b7d","#b8e186","#4d9221"))+
  theme(legend.text = element_text(colour="black",size=12),
        panel.grid.major = element_line(colour = "#E0E0E0"),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.3),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.3),
        axis.ticks.y = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_blank())+
  theme(panel.spacing.x = unit(1, "mm"))
a1

df_type<-read.delim("EUR_MR_target_bubble.tsv",header=T,check.names=F)
summary(df_type)
df_type$symbol <- factor(df_type$symbol, levels = rev(c('EPHX2','PRKCB','FOLH1','BACE2','EGFR','PLG','CD38','PARP10','MARK4','DRD2','CHRNE','DAGLB','APH1B','VKORC1','APP','PSEN2','PDE2A','NTRK1','GRIN3A','CACNA1D','F10','PDE5A','TAS2R60','GABBR1','NR3C1','ADRA2A','FPR1')))

df_type$group <- factor(df_type$group, levels = c('Amyloid','Tau','Inflammation','Neuromodulation','Microglia','Metabolism',"Vasculature",'Other'))


b1 <- ggplot(df_type, aes(x = group, y = symbol)) +
  geom_point(shape=15,aes(colour =type),size=5) +
  ylab("")+xlab("")+
  scale_y_discrete(position = "left")+
  scale_x_discrete(position = "top")+
  scale_color_manual(values=c("#42A5F5","#9575CD","#C0CA33","#A1887F","#26A69A","#78909C","#fdbf6f","#e57373","#FAFAFA"))+
  theme(legend.text = element_text(colour="black",size=12),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 12,angle=60,hjust=0,vjust=0),
        axis.text.y=element_text(colour="black", size = 12,angle=0,hjust=0,vjust=0.5))+
  theme(panel.spacing.x = unit(1, "mm"))
b1

pdf("./image/MR_target_IVW_bubble_plot-H15.pdf",width=6, height=8)
print(b1)
dev.off()

legend <- plot_grid(get_legend(pp), ncol = 1)
legend

a1 <- a1 + theme(legend.position = "none")
b1 <- b1 + theme(legend.position = "none")
pp <- pp + theme(legend.position = "none")

##combine b1 with a1
bb_plot <- plot_grid(b1, a1, align = "h", ncol = 2, axis = "t", rel_widths = c(6,4))
bb_plot

pdf("./image/MR_target_IVW_barbubble_plot.pdf",width=5, height=7.5)
print(bb_plot)
dev.off()
