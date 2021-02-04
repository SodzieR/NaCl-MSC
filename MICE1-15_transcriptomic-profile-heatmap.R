
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(readxl)
library(ComplexHeatmap)

# Setup -------------------------------------------------------------------
setwd('D:/!ZMRiI/!MSC_PAPER')

data <- read.delim('all_data_normalized_TPM.txt', header = TRUE) %>% 
  select(id,MICE1,MICE2,MICE3,MICE4,MICE5,MICE6,MICE7,MICE8,MICE9,MICE10,
         MICE11,MICE12,MICE13,MICE14,MICE15) %>% 
  column_to_rownames(var = 'id')

metadata <- data.frame(sample=c('MICE1','MICE2','MICE3','MICE4','MICE5',
                                'MICE6','MICE7','MICE8','MICE9','MICE10',
                                'MICE11','MICE12','MICE13','MICE14','MICE15'),
                       group = c(1,1,1,1,1,
                                 2,2,2,2,2,
                                 3,3,3,3,3))

# Load --------------------------------------------------------------------

# To properly calculate the z-score, the data has to be will be approximately 
# normally distributed and suitable for calculating z-scores. Z-score are 
# used to enhance the visualisation of the data (it helps to even the values out).
TPM_heatmap_in_log2 <- data %>% 
  as.matrix() %>% 
  log2()

TPM_heatmap_in_log2[is.infinite(TPM_heatmap_in_log2)] <- 0

column_names <- colnames(TPM_heatmap_in_log2)

TPM_heatmap_in_log2_scaled <- base::apply(TPM_heatmap_in_log2, 1, scale) %>% 
  t()

colnames(TPM_heatmap_in_log2_scaled) <- column_names

TPM_heatmap_in_log2_scaled[is.nan(TPM_heatmap_in_log2_scaled)] <- 0

col_fun1 = circlize::colorRamp2(c(min(TPM_heatmap_in_log2_scaled), 0, max(TPM_heatmap_in_log2_scaled)), 
                                c("green", "black", "red"))
#
Heatmap(TPM_heatmap_in_log2_scaled, name = 'MICE 1:15', col = col_fun1,
        cluster_rows = F, cluster_columns = FALSE, show_row_names = F,
        border = T, show_column_names = TRUE, column_title = 'Z-scored log2(TPM)',
        top_annotation = HeatmapAnnotation(
          empty = anno_empty(border = FALSE),
          foo = anno_block(gp = gpar(fill = 'black'), labels = c('Group 1', 'Group 2', 'Group 3'),
                           labels_gp = gpar(col = "white", fontsize = 10))),
        column_split = rep(1:3, each = 5)
        )


#Clustering problem - not enough ram

col_fun2 = circlize::colorRamp2(c(min(TPM_heatmap_in_log2_scaled), 0, max(TPM_heatmap_in_log2_scaled)),
                             c('blue','white','red'))

ann <- data.frame(metadata$group)
colnames(ann) <- c('group')
colours <- list('group' = c('1' = 'red', 
                            '2' = 'magenta',
                            '3' = 'blue'))

colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_width = unit(c(1, 4), 'cm'),
                            gap = unit(1, 'mm'))

temp2 <- Heatmap(TPM_heatmap_in_log2_scaled, name = 'expression', col = col_fun2,
        cluster_rows = F, cluster_columns = FALSE, show_row_names = F,
        border = T, show_column_names = TRUE, column_title = 'Z-scored log2(TPM)', 
        top_annotation=colAnn, row_dend_width = unit(4, 'cm'))

pdf("bwr_Defauult-clustering_z-scored_log2TPM_7x7.pdf", width=7, height=7)
draw(temp2, heatmap_legend_side="left", annotation_legend_side="right")
dev.off()


















