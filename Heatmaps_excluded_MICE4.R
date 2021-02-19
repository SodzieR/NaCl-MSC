
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ComplexHeatmap)
library(readxl)

# Setup -------------------------------------------------------------------
setwd('D:\\!ZMRiI\\!MSC_PAPER\\!!_excluded_MICE4')

data_TPM <- read.delim('all_data_normalized_TPM.txt', header = TRUE) %>% 
  dplyr::select(id,MICE1,MICE2,MICE3,MICE5,
                MICE6,MICE7,MICE8,MICE9,MICE10,
                MICE11,MICE12,MICE13,MICE14,MICE15)

data_1vs2_IPA <- read_excel('1(-MICE4)vs2_IPA.xls', sheet = 'analysis_ready_1vs2_IPA')

data_1vs3_IPA <- read_excel('1(-MICE4)vs3_IPA.xls', sheet = 'analysis_ready_1vs3_IPA')

wszystkie_geny_IPA <- data_1vs2_IPA %>% 
  full_join(data_1vs3_IPA) %>% 
  distinct(id, .keep_all = TRUE)

data <- data_TPM %>% 
  full_join(wszystkie_geny_IPA) %>% 
  arrange(symbol) %>% 
  distinct(symbol, .keep_all = TRUE) %>% 
  filter(!is.na(symbol)) %>% 
  column_to_rownames(var = 'symbol') %>% 
  select(-c(id))

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

#

col_fun2 = circlize::colorRamp2(c(-1.5, 0, 1.5),
                                c('blue','white','red'))


temp2 <- Heatmap(TPM_heatmap_in_log2_scaled, name = 'expression', col = col_fun2,
                 cluster_rows = T, cluster_columns = FALSE, show_row_names = F,
                 border = F, show_column_names = TRUE,column_title = 'Z-scored log2(TPM)',
                 row_dend_width = unit(3, 'cm')
)

pdf("XDDDD5_bwr_Defauult-clustering_z-scored_log2TPM_3x4.pdf", width=5, height=7)
draw(temp2, heatmap_legend_side="left")
dev.off()
#

gplots::heatmap.2(x=TPM_heatmap_in_log2_scaled, 
                  Colv = FALSE, 
                  Rowv = TRUE,
                  scale="none",
                  col="bluered", breaks=c(-1.5,-0.75,0.75,1.5),
                  trace="none",
                  dendrogram='none')











