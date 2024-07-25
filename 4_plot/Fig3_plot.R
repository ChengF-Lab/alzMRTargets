library(ggplot2)
library(cowplot)
library(dplyr)

setwd("/Users/houy2/Documents/projects/AD/AD_genetics/results/update_v9_Jul072023/MR_results/DT")

filename<-"AD#LD02_EUR_MR19_target_IVW_hp_Mar282024"
df_raw<-read.delim(sprintf("./heatmap/%s.tsv",filename),header=T,check.names=F)
summary(df_raw)

df <- filter(df_raw, Outcome != "AD.AAs_Kunkle_2020" & Outcome != "AD.AAsAPOE_Kunkle_2020")
target_list<-c('EPHX2','EGFR','PRKCB','FOLH1','BACE2','CD38','PLG','APH1B','CHRNE','PARP10','VKORC1','PTK2B','DRD2','ADRA1A','PDE5A','CACNA1D','NR3C1','GABBR1','FPR1')
df$symbol <- factor(df$symbol, levels = rev(target_list))

df$xQTL_type <- factor(df$xQTL_type, levels = c('pQTL_ROSMAP_Cortex.All',
                                                'pQTL_ROSMAP_Cortex.Ctl',
                                                'pQTL_Banner_Cortex',
                                                'eQTL_ROSMAP_Cortex',
                                                'eQTL_Mayo_Cortex',
                                                'eQTL_Metabrain_Cortex',
                                                'eQTL_Mayo_MetaAnalysis',
                                                'eQTL_Mayo_Cerebellum',
                                                'eQTL_Metabrain_Cerebellum',
                                                'eQTL_Metabrain_BasalGanglia',
                                                'eQTL_Metabrain_Hippocampus',
                                                'eQTL_Metabrain_SpinalCord'))
df$Outcome <- factor(df$Outcome, levels = c('AD_Bellenguez_2022',
                                            'AD.NoUKB23_Wightman_2021',
                                            'AD_Schwartzentruber_2021',
                                            'AD.FamilyHtx_Schwartzentruber_2021',
                                            'AD_Jansen_2019',
                                            'AD.Stage1&2_Kunkle_2019'))

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
  #axis.text.y=element_text(colour="black", size = 12))+
  theme(panel.spacing.x = unit(1, "mm"))
#coord_flip()
pa

# Adding the significance level stars using geom_text 
pp<- pa +
  geom_text(aes(label=sig),size=3.5,na.rm=TRUE,vjust =0.8) # you can play with the size
pp

pdf(sprintf("./heatmap/image/%s.pdf", filename),width=18, height=7)
#pdf("viral-critical.pdf",width=4.5, height=4)
print(pp) # need print plot first
dev.off()

#===================
df_pip<-read.csv("MR_target_EUR#19_PIPplot_Mar282024.csv",header=T,check.names=F)
summary(df_pip)
df_pip$symbol <- factor(df_pip$symbol, levels = rev(target_list))

pip<-ggplot(df_pip, aes(x = group, y = symbol,size=PIP, shape =group)) +
  geom_point(aes(colour = group),stroke = 1) +
  ylab("")+xlab("")+
  scale_y_discrete(position = "left")+
  scale_x_discrete(position = "top")+
  scale_shape_manual(values=c(7, 7))+
  scale_color_manual(values = c('FOCUS_DT' = "#C0CA33", 'coloc_DT' = "#fdbf6f"),na.value = "#FAFAFA") +
  theme(legend.text = element_text(colour="black",size=12),
        panel.grid.major = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.3),
        axis.ticks.y = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 12,angle=70,hjust=0,vjust=0),
        axis.text.y=element_blank())+
  theme(panel.spacing.x = unit(1, "mm"))
pip

pp <- pp + theme(legend.position = "none")
pip <- pip + theme(legend.position = "none")

##combine b1 with a1
pip_plot <- plot_grid(pp,pip, align = "h", ncol = 2, axis = "t", rel_widths = c(14,2))
#bb_plot<-plot_grid(bb_plot, legend, nrow = 1, rel_widths = c(10, 1.2))
pip_plot

pdf("./heatmap/image/MR_target_IVW_heatmapPIP_plot_Mar292024.pdf",width=12, height=6)
#pdf("viral-critical.pdf",width=4.5, height=4)
print(pip_plot) # need print plot first
dev.off()


#=============

df_target<-read.csv("MR_target_EUR#19_PIP_Mar282024.csv",header=T,check.names=F)
summary(df_target)


df_target$symbol <- factor(df_target$symbol, levels = rev(target_list))

df_target <- df_target %>%
  mutate(color_g = case_when(
    scoreSumZscores > 100 ~ ">100",
    scoreSumZscores > 50 & scoreSumZscores <= 100 ~ "50",
    scoreSumZscores > 0 & scoreSumZscores <= 50 ~ "0",
    scoreSumZscores > -20 & scoreSumZscores <= 0 ~ "-20",
    scoreSumZscores <= -20 ~ "-50"
  ))

df_target$color_g <- factor(df_target$color_g, levels = c('>100',
                                                          '50',
                                                          '0',
                                                          '-20',
                                                          '-50'))

a1 <- ggplot(df_target, aes(x = scoreSumZscores, y = symbol)) +
  geom_bar(stat="identity",aes(fill = color_g)) +
  ylab("")+xlab("")+
  scale_x_continuous(position = 'top')+
  scale_y_discrete(position = "left")+
  geom_vline(xintercept =0)+
  scale_fill_manual(values=c(">100"="#c51b7d",
                             "50"="#de77ae",
                             "0"= "#fde0ef",
                             "-20" = "#b8e186",
                             "-50" = "#4d9221"))+
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

df_target <- df_target %>%
  mutate(color_g2 = case_when(
    scoreMeanZscores > 2 ~ ">2",
    scoreMeanZscores > 1.5 & scoreMeanZscores <= 2 ~ "1.5",
    scoreMeanZscores > 1 & scoreMeanZscores <= 1.5 ~ "1",
    scoreMeanZscores > 0.5 & scoreMeanZscores <= 1 ~ "0.5",
    scoreMeanZscores > 0 & scoreMeanZscores <= 0.5 ~ "0",
    scoreMeanZscores > -0.5 & scoreMeanZscores <= 0 ~ "-0.5",
    scoreMeanZscores > -1 & scoreMeanZscores <= -0.5 ~ "-1",
    scoreMeanZscores > -1.5 & scoreMeanZscores <= -1 ~ "-1.5",
    scoreMeanZscores <= -1.5 ~ "-2"
  ))

df_target$color_g2 <- factor(df_target$color_g2, levels = c('>2',
                                                            '1.5',
                                                            '1',
                                                            '0.5',
                                                            '0',
                                                            '-0.5',
                                                            '-1',
                                                            "-1.5",
                                                            '-2'))

a2_s <- ggplot(df_target, aes(x = scoreMeanZscores, y = symbol)) +
  geom_bar(stat="identity",aes(fill = color_g2)) +
  ylab("")+xlab("")+
  scale_x_continuous(position = 'top')+
  scale_y_discrete(position = "left")+
  geom_vline(xintercept =0)+
  scale_fill_manual(values=c(">2"="#ec7014",
                             "1.5"= "#fe9929",
                             "1"= "#fec44f",
                             "0.5"= "#fee391",
                             "0"= "white",
                             "-0.5" = "#c6dbef",
                             "-1" = "#9ecae1",
                             "-1.5" = "#6baed6",
                             "-2" = "#4292c6"))+
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
a2_s


##pdf("./image/MR_target_IVW_bar_plot.pdf",width=6, height=8)
###pdf("viral-critical.pdf",width=4.5, height=4)
#print(a1) # need print plot first
#dev.off()

df_type<-read.delim("target_EUR#19_MOA_hist_Mar282024.tsv",header=T,check.names=F)
summary(df_type)
df_type$symbol <- factor(df_type$symbol, levels = rev(target_list))

df_type$group <- factor(df_type$group, levels = c('Amyloid','Tau','Inflammation','Neuromodulation','Microglia','Metabolism',"Vasculature"))


b1 <- ggplot(df_type, aes(x = group, y = symbol)) +
  geom_point(shape=15,aes(colour =group),size=5) +
  ylab("")+xlab("")+
  scale_y_discrete(position = "left")+
  scale_x_discrete(position = "top")+
  scale_color_manual(values = c('Amyloid' = "#42A5F5", 
                                'Inflammation' = "#9575CD", 
                                'Metabolism' = "#C0CA33", 
                                'Microglia' = "#A1887F", 
                                'Neuromodulation' = "#26A69A", 
                                'Tau' = "#fdbf6f", 
                                'Vasculature' = "#e57373"),
                     na.value = "#FAFAFA") +
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

pdf("./heatmap/image/MR_target_IVW_bubble_plot-H15.pdf",width=6, height=8)
#pdf("viral-critical.pdf",width=4.5, height=4)
print(b1) # need print plot first
dev.off()

legend <- plot_grid(get_legend(pp), ncol = 1)
legend

a1 <- a1 + theme(legend.position = "none")
a2_s <- a2_s + theme(legend.position = "none")
b1 <- b1 + theme(legend.position = "none")

##combine b1 with a1
bb_plot <- plot_grid(b1, a1,a2_s, align = "h", ncol = 3, axis = "t", rel_widths = c(5,2,2))
#bb_plot<-plot_grid(bb_plot, legend, nrow = 1, rel_widths = c(10, 1.2))
bb_plot

pdf("./heatmap/image/MR_target_IVW_barbubble_plot_Mar292024.pdf",width=6, height=6)
#pdf("viral-critical.pdf",width=4.5, height=4)
print(bb_plot) # need print plot first
dev.off()

