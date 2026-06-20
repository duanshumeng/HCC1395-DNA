#均匀性、稳定性评估
library(stringr)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(gridExtra)
source('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/Quartet_Multiomics_Ratio_Code/utils/theme_nature.r')
theme_nature_border <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.x = element_text(color = "black",face = "bold",size = 16),
        axis.title.y = element_text(color = "black",face = "bold",size = 16),
        #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
        legend.title = element_text(face = "bold",size = 16),
        legend.text = element_text(face = "bold",size = 14,colour = "black"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid.minor  = element_line(colour = NA),
        panel.grid.major.x   = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA)
  )+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.y = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))


library(ggplot2)
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Stability_and_homogenity/')
# SNV/Indel
#'#8E1025', '#F193B2', '#B4CB8F', '#598E44'
colors <- c('#598E44','#8E1025')
#基于已有的数据评估管间均匀性-------------------------------

rep.df <- read.csv('Hcc1395_Reproducibility_stats_total.csv')
rep.df.wes <- read.csv('Hcc1395_Reproducibility_stats_total.WES.csv')
rep.df.wes$rep_type <- 'Inter_vial'
rep.df$platform <- lapply(rep.df$source,function(x){
    if (str_count(x,'BGI_T10') == 2){
      return("BGI_T10")
    } else if (str_count(x,'BGI_T20') == 2){
      return("BGI_T20")
    } else if (str_count(x,'WUX') == 2){
      return("WUX")
    } else {
      return("")
    }
}) %>% unlist()

df_sub <- subset(rep.df,platform!="" & type %in% c('indels','SNVs'))

df_sub.p <- ggplot(df_sub, aes(x=type, y=Reproducibility, fill=type)) + 
  geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
  scale_fill_manual(values = colors) + 
  theme_nature_border + 
  theme(legend.position ='none'
  )+facet_grid(.~platform)+
  labs(x = NULL , y = 'Reproducibility')

ggsave('WGS.homogenity.BGIvsWUX.pdf',df_sub.p,width=3.15, height=2.7)

plot_box <- function(df){
  ggplot(df, aes(x=type, y=Reproducibility, fill=type)) + 
    geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
    scale_fill_manual(values = colors) + 
    theme_nature_border + 
    theme(legend.position ='none'
      )+
    labs(x = NULL , y = 'Reproducibility')
}

rep.df.sub <- subset(rep.df,type %in% c('indels','SNVs'),c('type','Reproducibility'))

p1 = plot_box(rep.df.sub)

ggsave('WGS.homogenity.interbottle.pdf',p1,width=5.15, height=4.7)

#计算标准偏差：
sd((subset(rep.df.sub,type=='indels')$Reproducibility)) %>% round(.,5)
sd((subset(rep.df.sub,type=='SNVs')$Reproducibility)) %>% round(.,5)
#计算均值：
median((subset(rep.df.sub,type=='indels')$Reproducibility)) %>% round(.,2)
median((subset(rep.df.sub,type=='SNVs')$Reproducibility)) %>% round(.,2)
#计算t检验的显著性差异：
t.test(subset(rep.df.sub,type=='indels')$Reproducibility,mu=0.92)
t.test(subset(rep.df.sub,type=='SNVs')$Reproducibility,mu=0.981)

rep.df.wes.sub <- subset(rep.df.wes,type %in% c('indels','SNVs'),c('type','Reproducibility'))
p2 = plot_box(rep.df.wes.sub)

ggsave('WES.homogenity.interbottle.pdf',p2,width=5.15, height=4.7)
#计算标准偏差：
sd((subset(rep.df.wes.sub,type=='indels')$Reproducibility)) %>% round(.,3)
sd((subset(rep.df.wes.sub,type=='SNVs')$Reproducibility)) %>% round(.,3)
#计算均值：
median((subset(rep.df.wes.sub,type=='indels')$Reproducibility))%>% round(.,2)
median((subset(rep.df.wes.sub,type=='SNVs')$Reproducibility))%>% round(.,2)

#基于序祯达产生的WES数据评估管间均匀性-------------------------------
rep.df.vial <- read.csv('SEQ_WES_Reproducibility_stats_total.csv')
rep.df.vial.2 <- read.csv('F1score_stats_total_rep_SEQ_2024-5-27.csv')
rep.df.vial.2.sub <- subset(rep.df.vial.2,!grepl('4-30|4-15',rep.df.vial.2$source))
rep.df.vial$rep_type <- lapply(rep.df.vial$source,function(x){
  x1 = str_split(str_split(x,'VS',simplify =T)[,1],'\\_',simplify = T)[,4]
  x2 = str_split(str_split(x,'VS',simplify =T)[,2],'\\_',simplify = T)[,4]
  if (x1 == x2){
    return('Inner_vial')
  } else {
    return('Inter_vial')
  }
}) %>% unlist()

rep.df.vial.2.sub$rep_type <- lapply(rep.df.vial.2.sub$source,function(x){
  x1 = str_split(str_split(x,'VS',simplify =T)[,1],'\\-',simplify = T)[,2]
  x2 = str_split(str_split(x,'VS',simplify =T)[,2],'\\-',simplify = T)[,2]
  if (x1 == x2){
    return('Inner_vial')
  } else {
    return('Inter_vial')
  }
}) %>% unlist()

rep.df.vial <- rbind(rep.df.vial,rep.df.vial.2.sub)

for (i in c('indels','SNVs')){
  aov_result <- aov(Reproducibility ~ rep_type, data = subset(rep.df.vial,type==i))
  summary(aov_result)
  print(i)
  F_value <- summary(aov_result)[[1]]$"F value"[1] %>% round(.,3)
  p_value <- summary(aov_result)[[1]]$"Pr(>F)"[1] %>% round(.,3)
  print(paste0('F_value:',F_value))
  print(paste0('p_value:',p_value))
}
#[1] "indels"
#[1] "F_value:0.006"
#[1] "p_value:0.938"
#[1] "SNVs"
#[1] "F_value:0.013"
#[1] "p_value:0.908"
#F值反映了组间变异与组内变异的比率。F值越小说明组间变异越小

plot_box_vial <- function(df){  
    p1 = ggplot(df, aes(x=rep_type, y=Reproducibility, fill=rep_type)) + 
    geom_boxplot(cex=0.5) +
    geom_signif(comparisons = list(c('Inner_vial','Inter_vial')),test.args =c(exact=FALSE),
                map_signif_level = T,textsize = 3)+
    scale_fill_manual(values = colors) + 
    theme_nature_border + 
    theme(legend.position ='none'
    ) +labs(x = NULL , y = 'Reproducibility')+facet_grid(.~type)
    return(p1)
    }

rep.df.vial.sub <- subset(rep.df.vial,type %in% c('indels','SNVs'),c('type','rep_type','Reproducibility'))

rep.df.vial.sub <- subset(rep.df.vial.sub,Reproducibility>0.92)
p.vial = plot_box_vial(rep.df.vial.sub)
ggsave('WES.homogenity.inter-inner.pdf',p.vial,width=8.15, height=5.7*0.7)

for (i in c('Inner_vial','Inter_vial')){
  print(i)
  #计算标准偏差：
  sd_i <- sd((subset(rep.df.vial.sub,type=='indels' & rep_type==i)$Reproducibility)) %>% round(.,3)
  sd_s <- sd((subset(rep.df.vial.sub,type=='SNVs' & rep_type==i)$Reproducibility))  %>% round(.,3)
  print(paste0('indel SD ',sd_i))
  print(paste0('snv SD ',sd_s))
  #计算均值：
  m_i <- median((subset(rep.df.vial.sub,type=='indels' & rep_type==i)$Reproducibility))  %>% round(.,2) 
  m_s <- median((subset(rep.df.vial.sub,type=='SNVs' & rep_type==i)$Reproducibility)) %>% round(.,2)
  print(paste0('indel Median ',m_i))
  print(paste0('snv Median ',m_s))
}

#[1] "Inner_vial"
#[1] "indel SD 0.028"
#[1] "snv SD 0.006"
#[1] "indel Median 0.95"
#[1] "snv Median 0.98"
#[1] "Inter_vial"
#[1] "indel SD 0.027"
#[1] "snv SD 0.006"
#[1] "indel Median 0.96"
#[1] "snv Median 0.98"
library(dplyr)
get_hom <- function(x,y){
  return(sqrt((x^2+y^2)))
}

#均匀性引入的不确定度-----
#SNV(0.009219544)
get_hom(0.006,0.007) 
#Indel(0.03201562)
get_hom(0.025,0.025) 


#基于DNA完整性评估短期稳定性-----------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Stability_and_homogenity')
ss_din <- read.csv("Stability_DIN_short-stage.csv",sep = '\t')
ss_din$Raw_sample <- ss_din$Sample
ss_din$Sample <- str_split(ss_din$Raw_sample,'\\-',simplify = T)[,1]
ss_din$Temperature <- str_split(ss_din$Raw_sample,'\\-',simplify = T)[,2] %>% paste0(.,'°C') %>% factor(.,levels = c("4°C", "25°C", "37°C"))
ss_din$Day <- str_split(ss_din$Raw_sample,'\\-',simplify = T)[,3] %>% as.numeric()
ss_din$DIN <- as.numeric(ss_din$DIN)



colors = c('#8d4891', '#f8e356', '#fe9536')
names(colors) = c("4°C", "25°C", "37°C")
p_bar <- function(df,v_type){
  df_sub <- subset(df,Sample==v_type)
  ebtop <- function(x){return(mean(x)+sd(x)/sqrt(length(x)))}
  ebbottom <- function(x){return(mean(x)-sd(x)/sqrt(length(x)))}
  p_bar_result <-
    ggplot(df_sub,aes(x = Day, y = DIN, fill = Temperature)) +
    theme_classic()+scale_fill_manual(values = colors,breaks = c("4°C", "25°C", "37°C")) +
    stat_summary(geom = "bar", fun = "mean", position = position_dodge(6), width=6)+
    stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(6), width=3) +
    #stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", method.args = list(alternative = "two.sided"))+
    labs(x=NULL)+ggtitle(paste0(v_type)) +
    theme_nature_border+
    scale_y_continuous(label=scientific,breaks=c(0,1,2,3,4,5,6,7,8,9,10))
  return(p_bar_result)
}
p_line <- function(df, v_type) {
  df_sub <- subset(df, Sample == v_type)
  
  ebtop <- function(x) { return(mean(x) + sd(x) / sqrt(length(x))) }
  ebbottom <- function(x) { return(mean(x) - sd(x) / sqrt(length(x))) }
  
  p_line_result <- ggplot(df_sub, aes(x = Day, y = DIN, group = Temperature,color =Temperature)) +
    #theme_classic() +
    scale_color_manual(values = colors) +
    stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5), size = 1) +
    stat_summary(geom = "point", fun = "mean", position = position_dodge(0.5), size = 3,alpha=0.6) +
    stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(0.5), width = 0.3) +
    labs(x = NULL) +
    ggtitle(paste0(v_type)) +
    theme_nature_border +
    scale_y_continuous(limits = c(0, 10), breaks = 0:10)+
    theme(
      legend.position = c(0.8, 0.2),  # 图例在图形内部的位置（x = 0.8, y = 0.2）
      legend.justification = c(1, 0))
    #ylim(c(0,10))
  
  return(p_line_result)
}

p_bar(ss_din,'HCC1395')
p_line(ss_din,'HCC1395')

hcc1395_pbar <- p_bar(ss_din,'HCC1395')
hcc1395_pline <- p_line(ss_din,'HCC1395')

#ggsave('HCC1395.Short-Stage-DIN.barplot.pdf',hcc1395_pbar,width=4.15, height=3.7)

hcc1395BL_pbar <- p_bar(ss_din,'HCC1395BL')
hcc1395BL_pline <- p_line(ss_din,'HCC1395BL')

#ggsave('HCC1395BL.Short-Stage-DIN.barplot.pdf',hcc1395BL_pbar,width=4.15, height=3.7)

F7_pbar <- p_bar(ss_din,'F7')
F7_pline <- p_line(ss_din,'F7')
#ggsave('F7.Short-Stage-DIN.barplot.pdf',F7_pbar,width=4.15, height=3.7)

get_box <- function(ss_din,smp,tmp){
  df = subset(ss_din,Sample==smp & Temperature==tmp)
  print(df)
  df$Day <- factor(df$Day,levels = c('0','10','20','30'))
  com_list = list(c('0','10'),c('0','20'),c('0','30'))
  p1 = ggplot(df, aes(x=Day, y=DIN,fill=Temperature)) + 
    geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
    geom_signif(comparisons = com_list,test.args =list(exact=TRUE),test = "t.test",
                map_signif_level = T,textsize = 3,y_position=c(9.5,9.6,9.7))+
    scale_fill_manual(values = colors) + 
    theme_nature_border + 
    theme(legend.position ='none'
    ) +ylim(8,10)+
    labs(x = 'Day', y = 'DIN')+ggtitle(paste0(smp,' ',tmp))
  return(p1)
}

#HCC1395
p4 <- get_box(ss_din,'HCC1395',"4°C")
p25 <- get_box(ss_din,'HCC1395',"25°C")
p37 <- get_box(ss_din,'HCC1395',"37°C")
hcc1395_pbox <- p4+p25+p37
#ggsave('HCC1395.Short-Stage-DIN.boxplot.pdf',hcc1395_pbox,width=7.15, height=4.7)

#HCC1395BL
p4 <- get_box(ss_din,'HCC1395BL',"4°C")
p25 <- get_box(ss_din,'HCC1395BL',"25°C")
p37 <- get_box(ss_din,'HCC1395BL',"37°C")
hcc1395BL_pbox <- p4+p25+p37
#ggsave('HCC1395BL.Short-Stage-DIN.boxplot.pdf',hcc1395BL_pbox,width=7.15, height=4.7)

#F7
get_box_F7 <- function(ss_din,smp,tmp){
  df = subset(ss_din,Sample==smp & Temperature==tmp)
  print(df)
  df$Day <- factor(df$Day,levels = c('0','10','20','30'))
  com_list = list(c('0','10'),c('0','20'),c('0','30'))
  p1 = ggplot(df, aes(x=Day, y=DIN,fill=Temperature)) + 
    geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
    geom_signif(comparisons = com_list,test.args =list(exact=FALSE),test = "t.test",
                map_signif_level = T,textsize = 3,y_position=c(8.7,8.8,8.9))+
    scale_fill_manual(values = colors) + 
    theme_nature_border + 
    theme(legend.position ='none'
    ) +ylim(7,10)+
    labs(x = 'Day', y = 'DIN')+ggtitle(paste0(smp,' ',tmp))
  return(p1)
}
p4 <- get_box_F7(ss_din,'F7',"4°C")
p25 <- get_box_F7(ss_din,'F7',"25°C")
p37 <- get_box_F7(ss_din,'F7',"37°C")
F7_pbox <- p4+p25+p37
#ggsave('F7.Short-Stage-DIN.boxplot.pdf',F7_pbox,width=7.15, height=4.7)

box_lst.ss = list()
line_lst.ss = list()

for (a in c('F7','HCC1395','HCC1395BL')){
  pbar <- p_line(ss_din,a)
  line_lst.ss[[a]] = pbar
  for (b in c("4°C", "25°C", "37°C")){
    n = paste0(a,' ', b)
    print(n)
    pbox <- get_box(ss_din,a,b)
    box_lst.ss[[n]] = pbox
  }
}

p_din.short.box <- plot_grid(plotlist = box_lst.ss, align = "h", 
                         nrow = 3)

p_din.short.line <- plot_grid(plotlist = line_lst.ss, align = "v", 
                             ncol = 3)
ggsave('DIN.short_stage.box.pdf',p_din.short.box,width=9, height=9)
ggsave('DIN.short_stage.line.pdf',p_din.short.line,width=8.15, height=5.7*0.7)

#线性回归模型
library(ggpmisc)
ss_din <- arrange(ss_din, (Temperature))
smp_tmp <- ss_din[,c('Sample','Temperature')] %>% unique()
p_lst = list()
df = data.frame(Sample=0,
                Temperature=0,
                T0.95=0,
                Beta1=0,
                Beta1.std = 0,
                Stability=0)

for (i in 1:dim(smp_tmp)[1]){
  smp = smp_tmp[i,'Sample']
  tmp = smp_tmp[i,'Temperature']
  test = subset(ss_din,Sample==smp & Temperature==tmp)
  fit <- lm(DIN~Day,test)
  fit.sum <- summary(fit)
  
  signature_level = 0.05
  degree_freedom = 12-2
  t0.95 <- qt(1-signature_level/2,degree_freedom)
  beta1 <- fit.sum$coefficients[2,1]
  beta1.s <- fit.sum$coefficients[2,2]
  sta <- ifelse(abs(beta1) < t0.95*beta1.s,'Yes','No')
  df.1 = data.frame(Sample=smp,
                  Temperature=tmp,
                  T0.95=t0.95,
                  Beta1=beta1,
                  Beta1.std = beta1.s,
                  Stability=sta)
  df = rbind(df,df.1)
  p <- ggplot(test, aes(x=Day, y=DIN)) +
    geom_point(color = "grey20",size = 1, alpha = 0.8)+ylim(8,10)
  #回归线
  #添加回归曲线
  p2 <- p + geom_smooth(formula = y ~ x, color = "red",
                        fill = "blue", method = "lm",se = T, level=0.95) +
    theme_bw() +
    stat_poly_eq(
      aes(label = paste(..eq.label..,..p.value.label..,sep = '~~~')),
      formula = y ~ x,  parse = TRUE,color="blue",
      size = 3.5, #公式字体大小
      label.x = 0.05,  #位置 ，0-1之间的比例
      label.y = 0.95) + 
    labs(title="test",x="Day" , y="DIN")+ggtitle(paste0(smp,' ',tmp))+theme_nature_border
  p_lst[[i]]=p2
  
}


p_din.short <- plot_grid(plotlist = p_lst, align = "h", 
          nrow = 3)
ggsave('DIN.short_stage.lm.pdf',p_din.short,width=9, height=9)

write.csv(df[-1,],file = 'DIN.short_stage.lm.csv',row.names = F)

#开瓶冻融稳定性评估--------------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Stability_and_homogenity')
tf_din <- read.csv("Thaw_Freeze_storage_DIN.csv",sep = '\t')
tf_din$Raw_sample <- tf_din$Sample %>% gsub('-N','_N',.)
tf_din$Sample <- str_split(tf_din$Raw_sample,'\\_',simplify = T)[,1]
tf_din$Temperature <- str_split(tf_din$Raw_sample,'\\_',simplify = T)[,2] %>% paste0(.,'°C')  %>% factor(.,levels = c("4°C", "-20°C", "-80°C"))
tf_din$Times <- str_extract(tf_din$Raw_sample, "(\\d*)$") %>% as.numeric()
tf_din$DIN <- as.numeric(tf_din$DIN)


tf_din <- dplyr::arrange(tf_din,Temperature)
colors.l = c('#8d4891', '#f8e356', '#fe9536')
names(colors.l) = c("4°C", "-20°C", "-80°C")
get_box.tf <- function(ss_din,smp,tmp){
  df = subset(ss_din,Sample==smp & Temperature==tmp)
  print(df)
  df$Times <- factor(df$Times,levels = c('0','5','10','15','20'))
  com_list = list(c('0','5'),c('0','10'),c('0','15'),c('0','20'))
  p1 = ggplot(df, aes(x=Times, y=DIN,fill=Temperature)) + 
    geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
    geom_signif(comparisons = com_list,test.args =list(exact=TRUE),test = "t.test",
                map_signif_level = T,textsize = 3,y_position=c(9.5,9.6,9.7,9.8))+
    scale_fill_manual(values = colors.l) + 
    theme_nature_border + 
    theme(legend.position ='none'
    ) +ylim(8,10)+
    labs(x = 'Times', y = 'DIN')+ggtitle(paste0(smp,' ',tmp))
  return(p1)
}

p_line.tf <- function(df, v_type) {
  df_sub <- subset(df, Sample == v_type)
  df_sub$Times <- factor(df_sub$Times,levels = c('0','5','10','15','20'))
  ebtop <- function(x) { return(mean(x) + sd(x) / sqrt(length(x))) }
  ebbottom <- function(x) { return(mean(x) - sd(x) / sqrt(length(x))) }
  
  p_line_result <- ggplot(df_sub, aes(x = Times, y = DIN, group = Temperature,color =Temperature)) +
    #theme_classic() +
    scale_color_manual(values = colors.l) +
    stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5), size = 1) +
    stat_summary(geom = "point", fun = "mean", position = position_dodge(0.5), size = 3,alpha=0.6) +
    stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(0.5), width = 0.3) +
    labs(x = NULL) +
    ggtitle(paste0(v_type)) +
    theme_nature_border +
    scale_y_continuous(limits = c(0, 10), breaks = 0:10)+
    theme(
      legend.position = c(0.8, 0.2),  # 图例在图形内部的位置（x = 0.8, y = 0.2）
      legend.justification = c(1, 0))
  #ylim(c(0,10))
  
  return(p_line_result)
}
box_lst = list()
line_lst = list()
for (a in c('F7','HCC1395','HCC1395BL')){
  pline = p_line.tf(tf_din,a)
  line_lst[[a]]=pline
  for (b in c("4°C", "-20°C", "-80°C")){
    n = paste0(a,' ', b)
    print(n)
    pbox <- get_box.tf(tf_din,a,b)
    box_lst[[n]] = pbox
  }
}

tf_plot <- grid.arrange(grobs = box_lst,ncol = 3)
tf.line.plot <- grid.arrange(grobs = line_lst,ncol = 3)

ggsave('DIN.thraw_fridge_stage.box.pdf',tf_plot,width=9, height=9)
ggsave('DIN.thraw_fridge_stage.line.pdf',tf.line.plot,width=8.15, height=5.7*0.7)

#线性回归---
library(ggpmisc)
smp_tmp <- tf_din[,c('Sample','Temperature')] %>% unique()
p_lst = list()
df = data.frame(Sample=0,
                Temperature=0,
                T0.95=0,
                Beta1=0,
                Beta1.std = 0,
                Stability=0)
for (i in 1:dim(smp_tmp)[1]){
  smp = smp_tmp[i,'Sample']
  tmp = smp_tmp[i,'Temperature']
  test = subset(tf_din,Sample==smp & Temperature==tmp)
  fit <- lm(DIN~Times,test)
  fit.sum <- summary(fit)
  
  signature_level = 0.05
  degree_freedom = 12-2
  t0.95 <- qt(1-signature_level/2,degree_freedom)
  beta1 <- fit.sum$coefficients[2,1]
  beta1.s <- fit.sum$coefficients[2,2]
  sta <- ifelse(abs(beta1) < t0.95*beta1.s,'Yes','No')
  df.1 = data.frame(Sample=smp,
                    Temperature=tmp,
                    T0.95=t0.95,
                    Beta1=beta1,
                    Beta1.std = beta1.s,
                    Stability=sta)
  df = rbind(df,df.1)
  p <- ggplot(test, aes(x=Times, y=DIN)) +
    geom_point(color = "grey20",size = 1, alpha = 0.8)+ylim(8,10)
  #回归线
  #添加回归曲线
  p2 <- p + geom_smooth(formula = y ~ x, color = "red",
                        fill = "blue", method = "lm",se = T, level=0.95) +
    theme_bw() +
    stat_poly_eq(
      aes(label = paste(..eq.label..,..p.value.label..,sep = '~~~')),
      formula = y ~ x,  parse = TRUE,color="blue",
      size = 3.5, #公式字体大小
      label.x = 0.05,  #位置 ，0-1之间的比例
      label.y = 0.95) + 
    labs(title="test",x="Times" , y="DIN")+ggtitle(paste0(smp,' ',tmp))+theme_nature_border
  p_lst[[i]]=p2
  
}


p_din.thraw_fridge <- plot_grid(plotlist = p_lst, align = "h", 
                         nrow = 3)
ggsave('DIN.thraw_fridge.lm.pdf',p_din.thraw_fridge,width=9, height=9)

write.csv(df[-1,],file = 'DIN.thraw_fridge.lm.csv',row.names = F)

#基于DNA完整性评估长期稳定性-----------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Stability_and_homogenity')
colors.l = c('#8d4891', '#f8e356', '#fe9536')
names(colors.l) = c("4°C", "-20°C", "-80°C")

ls_din <- read.csv("Stability_DIN_long-stage.csv")
ls_din$Raw_sample <- ls_din$Sample %>% gsub("","",.)
ls_din$Sample <- str_split(ls_din$Raw_sample,'\\_',simplify = T)[,1]
ls_din$Temperature <-  str_split(str_split(ls_din$Raw_sample,'\\_',simplify = T)[,2]
                                 ,'\\-N',simplify = T)[,1] %>% paste0(.,'°C') %>% factor(.,levels = c("4°C", "-20°C", "-80°C"))
ls_din$Day <- str_split(str_split(str_split(ls_din$Raw_sample,'\\_',simplify = T)[,2]
                                  ,'\\-N',simplify = T)[,2]
                        ,'\\-',simplify = T)[,2]

ls_din <- dplyr::arrange(ls_din,Temperature)
#ls_din$DIN <- as.numeric(ls_din$DIN)
m= c('0','1','3','5')
m_comp = list(c('0','1'),c('0','3'),c('0','5'))
ls_din$Day <- factor(ls_din$Day,levels = m)
get_box <- function(ss_din,smp,tmp,com_list){
  df = subset(ss_din,Sample==smp & Temperature==tmp)
  print(df)
  #com_list = list(c('0','1'))
  p1 = ggplot(df, aes(x=Day, y=DIN,fill=Temperature)) + 
    geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
    geom_signif(comparisons = com_list,test.args =list(exact=TRUE),test = "t.test",
                map_signif_level = T,textsize = 3,y_position=c(9.7))+
    scale_fill_manual(values = colors.l) + 
    theme_nature_border + 
    theme(legend.position ='none'
    ) +ylim(9,10)+
    labs(x = 'Month', y = 'DIN')+ggtitle(paste0(smp,' ',tmp))
  return(p1)
}

p4 <- get_box(ls_din,'HCC1395',"4°C",m_comp)
p20 <- get_box(ls_din,'HCC1395',"-20°C",m_comp)
p80 <- get_box(ls_din,'HCC1395',"-80°C",m_comp)

hcc1395.p.ls <- p4+p20+p80

p4 <- get_box(ls_din,'HCC1395BL',"4°C",m_comp)
p20 <- get_box(ls_din,'HCC1395BL',"-20°C",m_comp)
p80 <- get_box(ls_din,'HCC1395BL',"-80°C",m_comp)

hcc1395BL.p.ls <- p4+p20+p80

get_box_F7 <- function(ss_din,smp,tmp,com_list){
  df = subset(ss_din,Sample==smp & Temperature==tmp)
  print(df)
  #com_list = list(c('0','1'))
  p1 = ggplot(df, aes(x=Day, y=DIN,fill=Temperature)) + 
    geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
    geom_signif(comparisons = com_list,test.args =list(exact=TRUE),test = "t.test",
                map_signif_level = T,textsize = 3,y_position=c(8.7))+
    scale_fill_manual(values = colors.l) + 
    theme_nature_border + 
    theme(legend.position ='none'
    ) +ylim(8,10)+
    labs(x = 'Month', y = 'DIN')+ggtitle(paste0(smp,' ',tmp))
  return(p1)
}
p4 <- get_box_F7(ls_din,'F7',"4°C",m_comp)
p20 <- get_box_F7(ls_din,'F7',"-20°C",m_comp)
p80 <- get_box_F7(ls_din,'F7',"-80°C",m_comp)

F7.p.ls <- p4+p20+p80

p_ls.1 <- hcc1395.p.ls/hcc1395BL.p.ls+F7.p.ls
#ggsave('DIN.long_stage.M1.pdf',p_ls.1,width=9, height=9)
ggsave('DIN.long_stage.M5.pdf',p_ls.1,width=9, height=9)

p_line <- function(df, v_type) {
  df_sub <- subset(df, Sample == v_type)
  
  ebtop <- function(x) { return(mean(x) + sd(x) / sqrt(length(x))) }
  ebbottom <- function(x) { return(mean(x) - sd(x) / sqrt(length(x))) }
  
  p_line_result <- ggplot(df_sub, aes(x = Day, y = DIN, group = Temperature,color =Temperature)) +
    #theme_classic() +
    scale_color_manual(values = colors.l) +
    stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5), size = 1) +
    stat_summary(geom = "point", fun = "mean", position = position_dodge(0.5), size = 3,alpha=0.6) +
    stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(0.5), width = 0.3) +
    labs(x = NULL) +
    ggtitle(paste0(v_type)) +
    theme_nature_border +
    scale_y_continuous(limits = c(0, 10), breaks = 0:10)+
    theme(
      legend.position = c(0.8, 0.2),  # 图例在图形内部的位置（x = 0.8, y = 0.2）
      legend.justification = c(1, 0))
  #ylim(c(0,10))
  
  return(p_line_result)
}
line_lst.ls = list()

for (a in c('F7','HCC1395','HCC1395BL')){
  pbar <- p_line(ls_din,a)
  line_lst.ls[[a]] = pbar
}

p_din.long <- plot_grid(plotlist = line_lst.ls, align = "v", 
                                ncol = 3)
ggsave('DIN.long_stage.line.pdf',p_din.long,width=8.15, height=5.7*0.7)

#线性回归---
library(ggpmisc)
smp_tmp <- ls_din[,c('Sample','Temperature')] %>% unique()
ls_din$Day <- ls_din$Day %>% as.character() %>% as.numeric()
p_lst = list()
df = data.frame(Sample=0,
                Temperature=0,
                STDEV=0,
                T0.95=0,
                Beta1=0,
                Beta1.std = 0,
                Stability=0)
for (i in 1:dim(smp_tmp)[1]){
  smp = smp_tmp[i,'Sample']
  tmp = smp_tmp[i,'Temperature']
  test = subset(ls_din,Sample==smp & Temperature==tmp)
  test$Day <- as.numeric(test$Day)
  std=sd(test$DIN)*100
  fit <- lm(DIN~Day,test)
  fit.sum <- summary(fit)
  
  signature_level = 0.05
  degree_freedom = 4*(3-1)
  t0.95 <- qt(1-signature_level/2,degree_freedom)
  beta1 <- fit.sum$coefficients[2,1]
  beta1.s <- fit.sum$coefficients[2,2]
  sta <- ifelse(abs(beta1) < t0.95*beta1.s,'Yes','No')
  df.1 = data.frame(Sample=smp,
                    Temperature=tmp,
                    STDEV=std,
                    T0.95=t0.95,
                    Beta1=beta1,
                    Beta1.std = beta1.s,
                    Stability=sta)
  df = rbind(df,df.1)
  p <- ggplot(test, aes(x=Day, y=DIN)) +
    geom_point(color = "grey20",size = 1, alpha = 0.8)+ylim(8,10)
  #回归线
  #添加回归曲线
  p2 <- p + geom_smooth(formula = y ~ x, color = "red",
                        fill = "blue", method = "lm",se = T, level=0.95) +
    theme_bw() +
    stat_poly_eq(
      aes(label = paste(..eq.label..,..p.value.label..,sep = '~~~')),
      formula = y ~ x,  parse = TRUE,color="blue",
      size = 3.5, #公式字体大小
      label.x = 0.05,  #位置 ，0-1之间的比例
      label.y = 0.95) + 
    labs(title="test",x="Month" , y="DIN")+ggtitle(paste0(smp,' ',tmp))+theme_nature_border
  p_lst[[i]]=p2
  
}


p_din.long_stage <- plot_grid(plotlist = p_lst, align = "h", 
                                nrow = 3)

ggsave('DIN.long_stage.lm.M5.pdf',p_din.long_stage,width=9, height=9)

write.csv(df[-1,],file = 'DIN.long_stage.lm.csv',row.names = F)
#干冰运输稳定性评估---------------
ts_din <- read.csv("HCC1395.transport_stability.DIN.csv")
ts_din$Raw_sample <- ts_din$Sample
ts_din$Sample <- str_split(ts_din$Raw_sample,'\\_',simplify = T)[,4]
ts_din$Temperature <- str_split(ts_din$Raw_sample,'\\_',simplify = T)[,2]
ts_din$Day <- str_split(ts_din$Raw_sample,'\\_',simplify = T)[,3] %>% gsub('D','',.) %>% as.numeric()
ts_din$DIN <- as.numeric(ts_din$DIN)

colors = c('#fe9536')
names(colors) = c("CO2")
p_line <- function(df, v_type) {
  df_sub <- subset(df, Sample == v_type)
  
  ebtop <- function(x) { return(mean(x) + sd(x) / sqrt(length(x))) }
  ebbottom <- function(x) { return(mean(x) - sd(x) / sqrt(length(x))) }
  
  p_line_result <- ggplot(df_sub, aes(x = Day, y = DIN, group = Temperature,color =Temperature)) +
    #theme_classic() +
    scale_color_manual(values = colors) +
    stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5), size = 1) +
    stat_summary(geom = "point", fun = "mean", position = position_dodge(0.5), size = 3,alpha=0.6) +
    stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(0.5), width = 0.3) +
    labs(x = NULL) +
    ggtitle(paste0(v_type)) +
    theme_nature_border +
    scale_y_continuous(limits = c(0, 10), breaks = 0:10)+
    theme(
      legend.position = c(0.8, 0.2),  # 图例在图形内部的位置（x = 0.8, y = 0.2）
      legend.justification = c(1, 0))
  #ylim(c(0,10))
  
  return(p_line_result)
}
hcc1395_pline <- p_line(ts_din,'HCC1395')
hcc1395BL_pline <- p_line(ts_din,'HCC1395BL')
F7_pline <- p_line(ts_din,'D6')

p_din.transport <- F7_pline+hcc1395_pline+hcc1395BL_pline

ggsave('DIN.transport.line.pdf',p_din.transport,width=8.15, height=5.7*0.7)

#线性回归
smp_tmp <- ts_din[,c('Sample','Temperature')] %>% unique()
ts_din$Day <- ts_din$Day %>% as.character() %>% as.numeric()
p_lst = list()
df = data.frame(Sample=0,
                Temperature=0,
                STDEV=0,
                T0.95=0,
                Beta1=0,
                Beta1.std = 0,
                Stability=0)
for (i in 1:dim(smp_tmp)[1]){
  smp = smp_tmp[i,'Sample']
  tmp = smp_tmp[i,'Temperature']
  test = subset(ts_din,Sample==smp & Temperature==tmp)
  test$Day <- as.numeric(test$Day)
  std=sd(test$DIN)*100
  fit <- lm(DIN~Day,test)
  fit.sum <- summary(fit)
  
  signature_level = 0.05
  degree_freedom = 4*(3-1)
  t0.95 <- qt(1-signature_level/2,degree_freedom)
  beta1 <- fit.sum$coefficients[2,1]
  beta1.s <- fit.sum$coefficients[2,2]
  sta <- ifelse(abs(beta1) < t0.95*beta1.s,'Yes','No')
  df.1 = data.frame(Sample=smp,
                    Temperature=tmp,
                    STDEV=std,
                    T0.95=t0.95,
                    Beta1=beta1,
                    Beta1.std = beta1.s,
                    Stability=sta)
  df = rbind(df,df.1)
}

write.csv(df[-1,],file = 'DIN.transport.lm.csv',row.names = F)

#基于特征值准确性评估标准物质稳定性方案----------
p_line.v <- function(df, v_type) {
  df_sub <- subset(df, type == v_type)
  
  ebtop <- function(x) { return(mean(x) + sd(x) / sqrt(length(x))) }
  ebbottom <- function(x) { return(mean(x) - sd(x) / sqrt(length(x))) }
  
  p_line_result <- ggplot(df_sub, aes(x = Day, y = F1.score, group = Temperature,color =Temperature)) +
    #theme_classic() +
    scale_color_manual(values = colors.l) +
    stat_summary(geom = "line", fun = "mean", position = position_dodge(0.5), size = 1) +
    stat_summary(geom = "point", fun = "mean", position = position_dodge(0.5), size = 3,alpha=0.6) +
    stat_summary(geom = "errorbar", fun.min = ebbottom, fun.max = ebtop, position = position_dodge(0.5), width = 0.3) +
    labs(x = NULL) +
    ggtitle(paste0(v_type)) +
    theme_nature_border +
    scale_y_continuous(limits = c(0, 1))+
    theme(
      legend.position = c(0.8, 0.2),  # 图例在图形内部的位置（x = 0.8, y = 0.2）
      legend.justification = c(1, 0))
  #ylim(c(0,10))
  
  return(p_line_result)
}
#短期稳定性---------
f1_wgs.seq <- read.csv('F1score_stats_total_SEQ_2024-5-27.csv')
f1_wgs.s <- subset(f1_wgs.seq,grepl('4-15|4-30',f1_wgs.seq$source))

f1_wgs.s$Day <- str_split(f1_wgs.s$source,'[\\-,\\.]',simplify = T)[,4]
f1_wgs.s$Temperature <- '4°C'

f1_wes_seq.s <- read.csv('HCC1395_SEQ_WES_F1score_stats_total.csv')
f1_wes_seq.s$Day <- 0
f1_wes_seq.s$Temperature <- '4°C'

f1_wgs.s <- rbind(f1_wgs.s,f1_wes_seq.s)
f1_wgs.s <- subset(f1_wgs.s,type != 'records')
f1_wgs.s$type <- gsub('indels','InDels',f1_wgs.s$type)

s.indel <- p_line.v(f1_wgs.s,'InDels')
s.snv <- p_line.v(f1_wgs.s,'SNVs')
f1_wgs.s.p <- s.indel+s.snv
ggsave('F1.short_stage.line.pdf',f1_wgs.s.p,width=8.15, height=5.7*0.7)
#线性回归模型
library(ggpmisc)
smp_tmp <- f1_wgs.s[,c('type','Temperature')] %>% unique()
p_lst = list()
df = data.frame(type=0,
                Temperature=0,
                T0.95=0,
                STDEV=0,
                Beta1=0,
                Beta1.std = 0,
                Stability=0)

for (i in 1:dim(smp_tmp)[1]){
  smp = smp_tmp[i,'type']
  tmp = smp_tmp[i,'Temperature']
  test = subset(f1_wgs.s,type==smp & Temperature==tmp)
  std=sd(test$F1.score)*100
  fit <- lm(F1.score~Day,test)
  fit.sum <- summary(fit)
  
  signature_level = 0.05
  degree_freedom = 9-2
  t0.95 <- qt(1-signature_level/2,degree_freedom)
  beta1 <- fit.sum$coefficients[2,1]
  beta1.s <- fit.sum$coefficients[2,2]
  sta <- ifelse(abs(beta1) < t0.95*beta1.s,'Yes','No')
  df.1 = data.frame(type=smp,
                    Temperature=tmp,
                    T0.95=t0.95,
                    STDEV=std,
                    Beta1=beta1,
                    Beta1.std = beta1.s,
                    Stability=sta)
  df = rbind(df,df.1)
  p <- ggplot(test, aes(x=Day, y=F1.score)) +
    geom_point(color = "grey20",size = 1, alpha = 0.8)+ylim(0.6,1)
  #回归线
  #添加回归曲线
  p2 <- p + geom_smooth(formula = y ~ x, color = "red",
                        fill = "blue", method = "lm",se = T, level=0.95) +
    theme_bw() +
    stat_poly_eq(
      aes(label = paste(..eq.label..,..p.value.label..,sep = '~~~')),
      formula = y ~ x,  parse = TRUE,color="blue",
      size = 3.5, #公式字体大小
      label.x = 0.05,  #位置 ，0-1之间的比例
      label.y = 0.95) + 
    labs(title="test",x="Day" , y="F1.score")+ggtitle(paste0(smp,' ',tmp))+theme_nature_border
  p_lst[[i]]=p2
  
}


p_f1.short <- plot_grid(plotlist = p_lst, align = "h", 
                         nrow = 1)
ggsave('F1.short_stage.lm.pdf',p_f1.short,width=9*0.8, height=9*0.4)

write.csv(df[-1,],file = 'F1.short_stage.lm.csv',row.names = F)


#长期稳定性-----------
#2023-2月 (WUX、BGI) -> 0 
#2023-3月 (ARD、NVG、WGE) -> 1
f1_wgs <- read.csv('F1score_stats_total.deepsomatic-long-stability.csv')
f1_wgs <- subset(f1_wgs,!grepl('merge',f1_wgs$source))
f1_wgs$source <- gsub('WGS_WUX','WUX_WGS',f1_wgs$source)
f1_wgs$Batch <- str_split(f1_wgs$source,'\\_WGS',simplify = T)[,1]
f1_wgs <- subset(f1_wgs,Batch!='SL')
#f1_wgs <- subset(f1_wgs,grepl('TNseq.norm',f1_wgs$source) & grepl('2023',f1_wgs$source))


#2023-8月 (ELE) -> 6
f1_wes_ele <- read.csv('HCC1395_WES_F1score_stats_ELE.csv')
f1_wes_ele$date <- 6
#f1_wes_ele['F1.score'] <- 2 * (f1_wes_ele['precision'] * f1_wes_ele['recall']) / (f1_wes_ele['precision'] + f1_wes_ele['recall'])

#2024-1月 (SEQ) -> 11
f1_wes_seq <- read.csv('HCC1395_SEQ_WES_F1score_stats_total.csv')
f1_wes_seq$date <- 11

#2024-5 (SEQ) -> 15
f1_wgs.l <- subset(f1_wgs.seq,!grepl('4-15|4-30',f1_wgs.seq$source))
f1_wgs.l$date <- 15

#f1_wgs_2 <- rbind(subset(f1_wgs,type=='indels' & date=='2023-2',c('type','F1.score','date')) %>% sample_n(3),
#                 subset(f1_wgs,type=='SNVs'& date=='2023-2',c('type','F1.score','date')) %>% sample_n(3))
#f1_wgs_3 <- rbind(subset(f1_wgs,type=='indels' & date=='2023-3',c('type','F1.score','date')) %>% sample_n(3),
#                 subset(f1_wgs,type=='SNVs'& date=='2023-3',c('type','F1.score','date')) %>% sample_n(3))

#f1_wes_ele <- rbind(subset(f1_wes_ele,type=='indels',c('type','F1.score','date')) %>% sample_n(3),
#                    subset(f1_wes_ele,type=='SNVs',c('type','F1.score','date')) %>% sample_n(3)) 

#f1_wes_seq <- rbind(subset(f1_wes_seq,type=='indels',c('type','F1.score','date')) %>% sample_n(3),
#                    subset(f1_wes_seq,type=='SNVs',c('type','F1.score','date')) %>% sample_n(3)) 

#f1_long <- rbind(f1_wgs_2,f1_wgs_3) %>% rbind(.,f1_wes_ele) %>% rbind(.,f1_wes_seq)

f1_long <- f1_wgs

f1_long$Month <- as.character(f1_long$Month)
com_list=list(c('0','1'),c('0','32'),c('0','33'),c('0','34'),c('0','35'))
f1_long <- subset(f1_long,type!='records')
f1_long$F1.score <- as.numeric(f1_long$F1.score)
p1 = ggplot(f1_long, aes(x=Month, y=F1.score,fill=Month)) + 
  geom_boxplot(outlier.size = 0.8, lwd=0.3) + 
  geom_signif(comparisons = com_list,test.args =list(exact=FALSE),test = "t.test",
              map_signif_level = T,textsize = 3)+
  #scale_fill_manual(values = colors.l) + 
  theme_nature_border + 
  theme(legend.position ='none'
  ) +facet_grid(.~type)+xlim(c('0','1','32','33','34','35'))+
  scale_fill_brewer(palette = 'Dark2')+
  labs(x = 'Month', y = 'F1.score')

p1

f1_long$type <- gsub('indels','InDels',f1_long$type)
plm = list()
df = data.frame(Sample=0,
                  Temperature=0,
                  T0.95=0,
                  STDEV=0,
                  Beta1=0,
                  Beta1.std = 0,
                  Stability=0)
library(ggpmisc)
for (i in c('InDels','SNVs')){
  test=subset(f1_long,type==i)
  n=nrow(test)
  print(n)
  test$Month <- as.numeric(test$Month)
  std <- sd(test$F1.score)
  fit <- lm(F1.score~Month,test)
  fit.sum <- summary(fit)
  signature_level = 0.05
  degree_freedom = n-2
  t0.95 <- qt(1-signature_level/2,degree_freedom)
  beta1 <- fit.sum$coefficients[2,1]
  beta1.s <- fit.sum$coefficients[2,2]
  sta <- ifelse(abs(beta1) < t0.95*beta1.s,'Yes','No')
  smp=i
  tmp='-80°C'
  df.1 = data.frame(Sample=smp,
                    Temperature=tmp,
                    T0.95=t0.95,
                    STDEV=std,
                    Beta1=beta1,
                    Beta1.std = beta1.s,
                    Stability=sta)
  df <- rbind(df,df.1)
  p <- ggplot(test, aes(x=date, y=F1.score)) +
    geom_point(color = "grey20",size = 1, alpha = 0.8)+ylim(0.5,1)
  #回归线
  #添加回归曲线
  #p2 <- p + geom_smooth(formula = y ~ x, color = "red",
  #                      fill = "blue", method = "lm",se = T, level=0.95) +
  #  theme_bw() +
  #  stat_poly_eq(
  #    aes(label = paste(..eq.label..,..p.value.label..,sep = '~~~')),
  #    formula = y ~ x,  parse = TRUE,color="blue",
  #    size = 3.5, #公式字体大小
  #    label.x = 0.05,  #位置 ，0-1之间的比例
  #    label.y = 0.95)+scale_x_continuous(breaks =seq(0,15,3))+
  #  labs(x="Month" , y="F1.score")+ggtitle(paste0(smp,' ',tmp))+theme_nature_border
  #plm[[i]] = p2
}

#p_f1.long <- plot_grid(plotlist = plm, align = "h", 
#                        nrow = 1)
#ggsave('F1.long_stage.lm.pdf',p_f1.long,width=9*0.8, height=9*0.4)

write.csv(df[-1,],file = 'F1.long_stage.WGS.lm.csv',row.names = F)

f1_long$Day <- f1_long$Month %>% as.numeric()
f1_long$Temperature <- '-80°C'
s.indel <- p_line.v(f1_long,'InDels')
s.snv <- p_line.v(f1_long,'SNVs')
f1_long.p <- s.indel+s.snv
ggsave('F1.long_stage.WGS.line.pdf',f1_long.p,width=8.15, height=5.7*0.7)
#F1 score参考范围-----

get_95conf <- function(data){
  # 计算标准误差
  mean_data <- mean(data)
  se_data <- sd(data) / sqrt(length(data))
  
  # 计算t值 (95%的置信水平，双尾检验)
  alpha <- 0.05
  t_value <- qt(1 - alpha/2, df = length(data) - 1)
  
  # 计算置信区间
  lower_bound <- mean_data - t_value * se_data 
  upper_bound <- mean_data + t_value * se_data 
  
  conf_interval <- c(lower_bound %>% round(.,3), upper_bound %>% round(.,3))
  return(conf_interval)
}

data <- c(5, 8, 6, 7, 5, 10, 6, 8, 9, 7)


#WGS---
f1_wgs <- read.csv('HCC1395_WGS_F1score_stats_total.csv')
f1_wgs <- subset(f1_wgs,!grepl('merge',f1_wgs$source))

get_95conf(subset(f1_wgs,type=='indels')$F1.score) %>% print()
get_95conf(subset(f1_wgs,type=='SNVs')$F1.score) %>% print()

#median(subset(f1_wgs,type=='indels')$F1.score)
#median(subset(f1_wgs,type=='SNVs')$F1.score)

#WES---
f1_wes_ele <- read.csv('HCC1395_WES_F1score_stats_ELE.csv')

#f1_wes_ele['F1.score'] <- 2 * (f1_wes_ele['precision'] * f1_wes_ele['recall']) / (f1_wes_ele['precision'] + f1_wes_ele['recall'])

#2024-1月 (SEQ) -> 11
f1_wes_seq <- read.csv('HCC1395_SEQ_WES_F1score_stats_total.csv')


#2024-5 (SEQ) -> 15
f1_wes_seq.1 <- read.csv('F1score_stats_total_SEQ_2024-5-27.csv')

f1_wes <- rbind(f1_wes_ele,f1_wes_seq[,-11]) %>% rbind(.,f1_wes_seq.1[,-11])


get_95conf(subset(f1_wes,type=='indels')$F1.score) %>% print()
get_95conf(subset(f1_wes,type=='SNVs')$F1.score) %>% print()




