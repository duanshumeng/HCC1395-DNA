
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggsignif)
library(patchwork)
library(ggpubr)
library(cowplot)
library(ggsci)
options(scipen=200)
source('/Users/duanshumeng/生物信息/PGx_lab/毕业论文/Scripts/my_theme.R')
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/Pipeline_validation')



f1_s = read.csv("stats_total_compares_SEQC2.csv")
f1_s$source <- gsub('.bwa','',f1_s$source) %>% gsub('_1','',.) %>% paste0(.,'.SEQC2')
#f1 = read.csv("F1score_stats_total.csv")
f1 = read.csv("SEQC2_F1score_stats_total.validation.csv")
f1$source <- gsub('_stralka','.strelka',f1$source) %>% 
  gsub('strelka','strelka',.) %>% 
  gsub('.PASS','.PGx',.) %>% gsub('_deepsomatic','.deepsomatic',.) %>% gsub('_mutect2','.mutect2',.)

#f1$source <- str_split(f1$source,'\\.',simplify = TRUE)[,2]

f1_s['F1.score'] <- 2 * (f1_s['precision'] * f1_s['recall']) / (f1_s['precision'] + f1_s['recall'])
#f1['F1.score'] <- 2 * (f1['precision'] * f1['recall']) / (f1['precision'] + f1['recall'])

#合并SEQC2和PGx的数据绘制散点图
f1_s.sub <- f1_s[,c('type','source','F1.score')]
colnames(f1_s.sub) <- c('type','source','SEQC2.F1.score')
f1_s.sub$Sample <- str_split(f1_s.sub$source,'\\.',simplify = T)[,1]
f1_s.sub$Caller <- str_split(f1_s.sub$source,'\\.',simplify = T)[,2]
f1_s.clean <- subset(f1_s.sub,f1_s.sub$type!='records',c("type","SEQC2.F1.score","Sample","Caller"))

#f1_s.sub$source <- sub('.strelka.SEQC2','_1',f1_s.sub$source) %>% sub('.muTect2.SEQC2','_2',.)  %>% sub('.somaticSniper.SEQC2','_3',.) 
f1.sub <- f1[,c('type','source','F1.score')] 
colnames(f1.sub) <- c('type','source','PGx.F1.score')
f1.sub$Sample <- str_split(f1.sub$source,'\\.',simplify = T)[,1]
f1.sub$Caller <- str_split(f1.sub$source,'\\.',simplify = T)[,2]

# 1. 载入核心依赖包
library(tidyverse)
library(ggsci) # 提供 Nature 风格配色 (pal_npg)

# ==========================================
# 2. 载入并清洗数据
# ==========================================
f1_clean <- subset(f1.sub,f1.sub$type!='records',c("type","PGx.F1.score","Sample","Caller"))

f1_clean$type <- gsub('indels','InDel',f1_clean$type) %>% gsub('SNVs','SNV',.)

f1_s.clean$type <- gsub('indels','InDel',f1_s.clean$type) %>% gsub('SNVs','SNV',.)
# 过滤保留 indels 和 SNVs 面板，并统一规范大写以确保图表专业性
plot_data <- f1_clean %>%
  filter(type %in% c("InDel", "SNV")) %>%
  mutate(type = factor(type, levels = c("InDel", "SNV"), labels = c("InDel", "SNV")))

plot_data.s <- f1_s.clean %>%
  filter(type %in% c("InDel", "SNV")) %>%
  mutate(type = factor(type, levels = c("InDel", "SNV"), labels = c("InDel", "SNV")))


plot_data$Sample <- gsub('_T_1','',plot_data$Sample)

plot_data$Caller <- gsub('strelka','Strelka2',plot_data$Caller) %>% gsub('deepsomatic','Deepsomatic',.) %>% gsub('mutect2','Mutect2',.)

plot_data.s$Caller <- gsub('strelka','Strelka2',plot_data.s$Caller) %>% gsub('somaticSniper','SomaticSniper',.) %>% gsub('muTect2','Mutect2',.)

# ==========================================
# 3. 核心步骤：构建不同 Panel 的个性化阈值数据
# ==========================================
thresholds <- data.frame(
  type = factor(c("InDel", "SNV"), levels = c("InDel", "SNV")),
  cutoff = c(0.70, 0.85)
)

# ==========================================
# 5. ggplot2 绘图执行
# ==========================================
p1 <- ggplot(plot_data, aes(x = Sample, y = PGx.F1.score, fill = Caller)) +
  # 绘制分组条形图 (position_dodge 让同一小样本的各软件并排显示)
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.2) +
  
  # 按 type 进行分面 (每个变异类型拥有独立的 Panel)
  facet_wrap(~ type, ncol = 2, scales = "fixed") +
  
  # 关键步骤：映射不同面板对应的阈值线
  geom_hline(data = thresholds, aes(yintercept = cutoff), 
             color = "black", linetype = "dashed", linewidth = 0.5) +
  
  # 在阈值线上方添加标签提示
  geom_text(data = thresholds, aes(x = 0.6, y = cutoff + 0.02, label = paste("Cutoff:", cutoff)),
            inherit.aes = FALSE, size = 3, color = "black", fontface = "italic", hjust = 0) +
  
  # 视觉修饰
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_fill_npg() +  # 调用 Nature Publishing Group 官方美学配色
  labs(
    x = "Sample",
    y = "F1 Score",
    fill = "Caller",
  ) +
  my_theme

# 打印图表
print(p1)


p2 <- ggplot(plot_data.s, aes(x = Sample, y = SEQC2.F1.score, fill = Caller)) +
  # 绘制分组条形图 (position_dodge 让同一小样本的各软件并排显示)
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.2) +
  
  # 按 type 进行分面 (每个变异类型拥有独立的 Panel)
  facet_wrap(~ type, ncol = 2, scales = "fixed") +
  
  # 关键步骤：映射不同面板对应的阈值线
  geom_hline(data = thresholds, aes(yintercept = cutoff), 
             color = "black", linetype = "dashed", linewidth = 0.5) +
  
  # 在阈值线上方添加标签提示
  geom_text(data = thresholds, aes(x = 0.6, y = cutoff + 0.02, label = paste("Cutoff:", cutoff)),
            inherit.aes = FALSE, size = 3, color = "black", fontface = "italic", hjust = 0) +
  
  # 视觉修饰
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  scale_fill_npg() +  # 调用 Nature Publishing Group 官方美学配色
  labs(
    x = "Sample",
    y = "F1 Score",
    fill = "Caller",
  ) +
  my_theme

print(p2)
# 保存符合期刊分辨率要求的 PDF 格式文件
ggsave("F1_Barplot.pipeline_validation.pdf", plot = p, width = 6.5, height = 4.5, units = "in", dpi = 300)

