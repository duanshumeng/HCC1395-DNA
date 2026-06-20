setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Stability_and_homogenity')
library(dplyr)
library(stringr)
library(patchwork)
df_ji <- read.csv('All.Jacard_index_bedtools_isec.csv')

df_ji$source <- gsub('WGS_WUX','WUX_WGS',df_ji$source)
df_ji$smp1 = str_split(df_ji$source,'VS',simplify = T)[,1]
df_ji$batch1 = ifelse(grepl('2021',df_ji$smp1),'2021','2023')
df_ji$smp2 = str_split(df_ji$source,'VS',simplify = T)[,2]
df_ji$batch2 = str_split(df_ji$smp2,'_',simplify = T)[,3]
df_ji$batch2 = ifelse(grepl('2021',df_ji$smp2),'2021','2023')
df_ji <- subset(df_ji,batch1=='2023' & batch2=='2023')
df_ji$type <- gsub('indels','Indel',df_ji$type) %>% gsub('snps','SNV',.) 

get_batch_compre <- function(df_ji){
  df_ji$Library1 =  str_split(df_ji$smp1,'_WGS',simplify = T)[,1]
  df_ji$Library2 =  str_split(df_ji$smp2,'_WGS',simplify = T)[,1]
  
  #WUX:
  df_ji.wux = subset(df_ji,(df_ji$Library1 == 'WUX' | df_ji$Library2 == 'WUX') & (!grepl('BGI_T10|BGI_T20',df_ji$Library1) & !grepl('BGI_T10|BGI_T20',df_ji$Library2)))
  df_ji.wux$compare <- ifelse(df_ji.wux$Library1==df_ji.wux$Library2,'inner_batch','intra_batch')
  
  
  #BGI_T10
  df_ji.t10 = subset(df_ji,(df_ji$Library1 == 'BGI_T10' | df_ji$Library2 == 'BGI_T10') & (!grepl('WUX|BGI_T20',df_ji$Library1) & !grepl('WUX|BGI_T20',df_ji$Library2)))
  df_ji.t10$compare <- ifelse(df_ji.t10$Library1==df_ji.t10$Library2,'inner_batch','intra_batch')
  
  #BGI_T20
  df_ji.t20 = subset(df_ji,(df_ji$Library1 == 'BGI_T20' | df_ji$Library2 == 'BGI_T20') & (!grepl('WUX|BGI_T10',df_ji$Library1) & !grepl('WUX|BGI_T10',df_ji$Library2)))
  df_ji.t20$compare <- ifelse(df_ji.t20$Library1==df_ji.t20$Library2,'inner_batch','intra_batch')
  
  
  
  get_violin <- function(df_ji,lib){
    p1 <- ggplot(df_ji,aes(x=compare,y=JI))+
      geom_violin(aes(fill = compare), width = 0.9,alpha=0.5)+
      geom_boxplot(aes(fill = compare), width = 0.2,alpha=0.7) +  # 箱线图展示分布
      geom_jitter(aes(color = compare), width = 0.1, size = 3) +  # 散点展示个体值
      labs(
        title = paste("Jaccard by",lib),
        x = '',
        y = "Jaccard"
      ) +
      my_theme+ guides(fill=FALSE,color=FALSE)+
      stat_compare_means(method = "wilcox.test",
                         #aes(label = "p.format"), #显示方式
                         label = "p.signif",
                         size = 5,label.y=c(0.85),
                         comparisons = list(c('inner_batch','intra_batch')))+facet_grid(.~type)
    return(p1)
  }
  
  p1=get_violin(df_ji.wux,'WUX')
  p2=get_violin(df_ji.t10,'BGI_T10')
  p3=get_violin(df_ji.t20,'BGI_T20')
  
  return(p1+p2+p3)
  
}
df_ji.all <- subset(df_ji,!grepl('HCR',df_ji$source))
ji_all <- get_batch_compre(df_ji.all)
ggsave('Batch_reproducibility_all.pdf',ji_all,width =8.15*2.5 ,height =5.7,dpi = 300)


df_ji.hcr <- subset(df_ji,str_count(df_ji$source, "HCR")==2)
ji_hc <- get_batch_compre(df_ji.hcr)
ggsave('Batch_reproducibility_HC.pdf',ji_hc,width =8.15*2.5 ,height =5.7,dpi = 300)

#batch saturation------
JI_df <- read_csv('JI_cutdown_batch.csv')
JI_df$Type <- str_split(JI_df$sample,'\\.|\\_',simplify = T)[,5] 
JI_df$Datasize <- str_split(JI_df$sample,'\\.|\\_',simplify = T)[,4]  %>% as.numeric()
JI_df$recall <- JI_df$n_intersections/(JI_df$n_intersections+JI_df$FN)
JI_df$precision <- JI_df$n_intersections/(JI_df$n_intersections+JI_df$FP)
JI_df$F1 <- 2*JI_df$recall*JI_df$precision/(JI_df$precision+JI_df$recall)
JI_df[,c('Type','recall','precision','F1','Datasize')]
JI_df_sorted <- JI_df[order(JI_df$Datasize), c('Type','recall','precision','F1','Datasize')]

data_long <- pivot_longer(
  JI_df_sorted,
  cols = c(recall, precision, F1),
  names_to = "Metric",
  values_to = "Value"
)

# 绘制趋势图
f1_p <- ggplot(data_long, aes(x = Datasize, y = Value, color = Metric, group = Metric)) +
  geom_line(linewidth = 1) +  # 绘制线条
  geom_point(size = 2,alpha=0.7) +        # 添加数据点
  scale_color_manual(values = c("recall" = "#E41A1C", 
                                "precision" = "#377EB8", 
                                "F1" = "#4DAF4A")) +
  labs(title = "",
       x = "Datasize",
       y = "Value") +
  scale_x_continuous(breaks = JI_df_sorted$Datasize)+facet_grid(.~Type)+my_theme

ggsave('Batch_cutdown_F1.png',f1_p,width =8.15*1.5 ,height =5.7,dpi = 300)


