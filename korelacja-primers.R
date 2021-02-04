
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(readxl)
# rstatix is for kendall's corr matrix
library(rstatix)
# library(Hmisc) 
library(corrplot)
library(circlize)

# Setup -------------------------------------------------------------------
setwd('D:/!ZMRiI/!MSC_PAPER/')

#####################################################
#####################################################

#raw_read_counts <- read_excel('raw_read_counts.xlsx')

raw_read_counts <- read.delim('all_data_normalized_TPM.txt', header = TRUE)

primers <- toupper(c('EF1a',
                     'EF1a',
                     'GAPDH',
                     'GAPDH1',
                     'GAPDH2',
                     'GAPDH2',
                     'ZO1',
                     'Clca1',
                     'Gob5',
                     'Cldn1',
                     'Cldn4',
                     'Cldn7',
                     'Cldn18',
                     'Ocln',
                     'MMP2',
                     'MMP3',
                     'MMP7',
                     'MMP9' ,
                     'MMP12',
                     'TIMP1',
                     'TIMP2',
                     'IL1B',
                     'IL3',
                     'IL4',
                     'IL5',
                     'IL6',
                     'Cxcl1',
                     'Cxcl2',
                     'IL13',
                     'IL10',
                     'IL17A',
                     'IL17',
                     'IL22',
                     'IL23a',
                     'IL25',
                     'IL33',
                     'TSLP',
                     'TGFB1',
                     'IFNY',
                     'TNFa',
                     'Csf2',
                     'Csf3',
                     'Col1a2',
                     'Col3a2',
                     'Col5a1',
                     'Col5a2',
                     'Col5a3',
                     'Muc5a',
                     'Muc5b',
                     'Tnfsf13b',
                     'Tnfsf13',
                     'Tnfrsf17',
                     'Tnfsf13c'
))

######################################################
######################################################

read_counts_to_HGNC <- (AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db,
                                              key=raw_read_counts$id, 
                                              columns="SYMBOL",
                                              keytype="ENSEMBL")) %>% 
  rename('id' = ENSEMBL) %>% 
  mutate('SYMBOL' = toupper(SYMBOL))

raw_read_HGNC <- raw_read_counts %>% 
  full_join(read_counts_to_HGNC) %>% 
  select(id, SYMBOL, everything()) %>% 
  select(-c('X','X.1'))
# Load --------------------------------------------------------------------

genes_from_hypothesis <- c('TGFB1','TGFB2','TGFB3','FOXP3','SMAD7','IL2RB',
                           'IL2RG','TGFBR1','TGFBR2','TGFBR3','PTGES2','CCL2',
                           'IL6','IL15','IL15RA','IL10','IFNG','IL4','IL5',
                           'IL13','JAK3','STAT1','STAT5A','STAT3','JAK1','MTOR',
                           'AKTIP','AKT3', primers) %>% 
  sort()

genes_from_hypothesis_no_primers <- c('TGFB1','TGFB2','TGFB3','FOXP3','SMAD7','IL2RB',
                           'IL2RG','TGFBR1','TGFBR2','TGFBR3','PTGES2','CCL2',
                           'IL6','IL15','IL15RA','IL10','IFNG','IL4','IL5',
                           'IL13','JAK3','STAT1','STAT5A','STAT3','JAK1','MTOR',
                           'AKTIP','AKT3') %>% 
  sort()


data_1 <- raw_read_HGNC %>% 
  select(SYMBOL,MICE1,MICE2,MICE3,MICE4,MICE5) %>% 
  filter(SYMBOL %in% genes_from_hypothesis) %>% 
  column_to_rownames(var = 'SYMBOL') %>% t()

data_2 <- raw_read_HGNC %>% 
  select(SYMBOL,MICE6,MICE7,MICE8,MICE9,MICE10) %>% 
  filter(SYMBOL %in% genes_from_hypothesis) %>% 
  column_to_rownames(var = 'SYMBOL') %>% t()

data_3 <- raw_read_HGNC %>% 
  select(SYMBOL,MICE11,MICE12,MICE13,MICE14,MICE15) %>% 
  filter(SYMBOL %in% genes_from_hypothesis) %>% 
  column_to_rownames(var = 'SYMBOL') %>% t()


# Correlation matrices ----------------------------------------------------
# Data_1 ------------------------------------------------------------------
data_1_cor_r <- data.frame(data_1) %>% 
  cor_mat(method = 'kendall')

rownames(data_1_cor_r) <- NULL

data_1_cor_p <- data_1_cor_r %>% 
  cor_get_pval()

rownames(data_1_cor_p) <- NULL


data_1_cor_r_found_primers <- names(data_1_cor_r)[(names(data_1_cor_r) %in% primers)]

data_1_cor_r2 <- data_1_cor_r %>% 
  select('rowname', data_1_cor_r_found_primers) %>% 
  dplyr::filter(rowname %in% genes_from_hypothesis_no_primers)

data_1_cor_p_found_primers <- names(data_1_cor_p)[(names(data_1_cor_p) %in% primers)]

data_1_cor_p2 <- data_1_cor_p %>% 
  select('rowname', data_1_cor_p_found_primers) %>% 
  dplyr::filter(rowname %in% genes_from_hypothesis_no_primers)

data_1_cor_r_corrplot <- data_1_cor_r2 %>%
  column_to_rownames(var = 'rowname')

data_1_cor_p_corrplot <- data_1_cor_p2 %>%
  column_to_rownames(var = 'rowname')

#
pdf("group1-primers.pdf", width=7, height=7)
corrplot(as.matrix(data_1_cor_r_corrplot), type="lower", order="original", tl.col = 'black',
         p.mat = as.matrix(data_1_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'Primers - group1', method = 'square', tl.srt = 45,mar=c(0,0,1,0))
dev.off()

# Data 2 ------------------------------------------------------------------
data_2_cor_r <- data.frame(data_2) %>% 
  cor_mat(method = 'kendall')

rownames(data_2_cor_r) <- NULL

data_2_cor_p <- data_2_cor_r %>% 
  cor_get_pval()

rownames(data_2_cor_p) <- NULL


data_2_cor_r_found_primers <- names(data_2_cor_r)[(names(data_2_cor_r) %in% primers)]

data_2_cor_r2 <- data_2_cor_r %>% 
  select('rowname', data_2_cor_r_found_primers) %>% 
  dplyr::filter(rowname %in% genes_from_hypothesis_no_primers)

data_2_cor_p_found_primers <- names(data_2_cor_p)[(names(data_2_cor_p) %in% primers)]

data_2_cor_p2 <- data_2_cor_p %>% 
  select('rowname', data_2_cor_p_found_primers) %>% 
  dplyr::filter(rowname %in% genes_from_hypothesis_no_primers)

data_2_cor_r_corrplot <- data_2_cor_r2 %>%
  column_to_rownames(var = 'rowname')

data_2_cor_p_corrplot <- data_2_cor_p2 %>%
  column_to_rownames(var = 'rowname')

#
pdf("group2-primers.pdf", width=7, height=7)
corrplot(as.matrix(data_2_cor_r_corrplot), type="lower", order="original", tl.col = 'black',
         p.mat = as.matrix(data_2_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'Primers - group2', method = 'square', tl.srt = 45,mar=c(0,0,1,0))
dev.off()


# Data 3 ------------------------------------------------------------------
data_3_cor_r <- data.frame(data_3) %>% 
  cor_mat(method = 'kendall')

rownames(data_3_cor_r) <- NULL

data_3_cor_p <- data_3_cor_r %>% 
  cor_get_pval()

rownames(data_3_cor_p) <- NULL


data_3_cor_r_found_primers <- names(data_3_cor_r)[(names(data_3_cor_r) %in% primers)]

data_3_cor_r2 <- data_3_cor_r %>% 
  select('rowname', data_3_cor_r_found_primers) %>% 
  dplyr::filter(rowname %in% genes_from_hypothesis_no_primers)

data_3_cor_p_found_primers <- names(data_3_cor_p)[(names(data_3_cor_p) %in% primers)]

data_3_cor_p2 <- data_3_cor_p %>% 
  select('rowname', data_3_cor_p_found_primers) %>% 
  dplyr::filter(rowname %in% genes_from_hypothesis_no_primers)

data_3_cor_r_corrplot <- data_3_cor_r2 %>%
  column_to_rownames(var = 'rowname')

data_3_cor_p_corrplot <- data_3_cor_p2 %>%
  column_to_rownames(var = 'rowname')

#
pdf("group3-primers.pdf", width=7, height=7)
corrplot(as.matrix(data_3_cor_r_corrplot), type="lower", order="original", tl.col = 'black',
         p.mat = as.matrix(data_3_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'Primers - group3', method = 'square', tl.srt = 45,mar=c(0,0,1,0))
dev.off()
