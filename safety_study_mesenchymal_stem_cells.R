
# wd ----------------------------------------------------------------------
setwd('D:\\!ZMRiI\\!MSC_PAPER')

# Dependencies ------------------------------------------------------------

library(clusterProfiler)
library(AnnotationHub)
library(readxl)
library(ggplot2)
library(RColorBrewer)

# Pre-analysed GSEA datas -------------------------------------------------

preanalyzded_1vs2 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\!GSEA_FIGRUES\\FULL-GO-GSEA-1VS2.xlsx') %>% 
  dplyr::mutate('-log10(padj)' = -log10(p.adjust))

preanalyzded_1vs3 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\!GSEA_FIGRUES\\FULL-GO-GSEA-1VS3.xlsx') %>% 
  dplyr::mutate('-log10(padj)' = -log10(p.adjust))

preanalyzded_2vs3 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\!GSEA_FIGRUES\\FULL-GO-GSEA-2VS3.xlsx') %>% 
  dplyr::mutate('-log10(padj)' = -log10(p.adjust))


# FIGURE PLOTs
###############################################
# , hjust = 0.0001 - put it after position argument in geom text
###############################################
# 1vs2
pdf('GSEA_plots.pdf', width = 10)
(GO_GSEA_1vs2_df_plot <- ggplot(preanalyzded_1vs2, aes(factor(factor(Description, levels = Description), levels = Description), `-log10(padj)`, fill = enrichmentScore)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
    geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
    coord_flip() +
    ggtitle('GO_GSEA 1vs2') +
    labs(y= '-log10(padj)', x = "Description") +
    theme_classic())
(GO_GSEA_1vs2_df_plot_NES <- ggplot(preanalyzded_1vs2, aes(factor(factor(Description, levels = Description), levels = Description), `-log10(padj)`, fill = NES)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
  coord_flip() +
  ggtitle('GO_GSEA 1vs2') +
  labs(y= '-log10(padj)', x = "Description") +
  theme_classic())

# 1vs3
(GO_GSEA_1vs3_df_plot <- ggplot(preanalyzded_1vs3, aes(factor(Description, levels = Description), `-log10(padj)`, fill = enrichmentScore)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
  coord_flip() +
  ggtitle('GO_GSEA 1vs3') +
  labs(y= '-log10(padj)', x = "Description") +
  theme_classic())
(GO_GSEA_1vs3_df_plot_NES <- ggplot(preanalyzded_1vs3, aes(factor(factor(Description, levels = Description), levels = Description), `-log10(padj)`, fill = NES)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
  coord_flip() +
  ggtitle('GO_GSEA 1vs3') +
  labs(y= '-log10(padj)', x = "Description") +
  theme_classic())
# 2vs3 
(GO_GSEA_2vs3_df_plot <- ggplot(preanalyzded_2vs3, aes(factor(Description, levels = Description), `-log10(padj)`, fill = enrichmentScore)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
  coord_flip() +
  ggtitle('GO_GSEA 2vs3') +
  labs(y= '-log10(padj)', x = "Description") +
  theme_classic())
(GO_GSEA_2vs3_df_plot_NES <- ggplot(preanalyzded_2vs3, aes(factor(factor(Description, levels = Description), levels = Description), `-log10(padj)`, fill = NES)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
  coord_flip() +
  ggtitle('GO_GSEA 2vs3') +
  labs(y= '-log10(padj)', x = "Description") +
  theme_classic())
dev.off()
# Setup -------------------------------------------------------------------

Treg_Macrophages_M2 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\makrofagi-i-inne.xlsx')

data_1vs2 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\1vs2.xlsx') %>% 
  dplyr::filter(!is.na(id))
data_1vs2$log2FoldChange <- as.numeric(data_1vs2$log2FoldChange)
data_1vs2$padj <- as.numeric(data_1vs2$padj)

data_1vs2 <- data_1vs2 %>% 
  arrange(padj) %>% 
  dplyr::distinct(gene_ID, .keep_all = TRUE)
  
# -

data_1vs3 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\1vs3.xlsx') %>% 
  dplyr::filter(!is.na(id))
data_1vs3$log2FoldChange <- as.numeric(data_1vs3$log2FoldChange)
data_1vs3$padj <- as.numeric(data_1vs3$padj)

data_1vs3 <- data_1vs3 %>% 
  arrange(padj) %>% 
  dplyr::distinct(gene_ID, .keep_all = TRUE)

# -

data_2vs3 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\2vs3.xlsx') %>% 
  dplyr::filter(!is.na(id))
data_2vs3$log2FoldChange <- as.numeric(data_2vs3$log2FoldChange)
data_2vs3$padj <- as.numeric(data_2vs3$padj)

data_2vs3 <- data_2vs3 %>% 
  arrange(padj) %>% 
  dplyr::distinct(gene_ID, .keep_all = TRUE)

# Load --------------------------------------------------------------------

# Grab Mice Genome Annotation
hub <- AnnotationHub()
q <- query(hub, 'org.Mm.eg.db')
id <- q$ah_id[length(q)]
miceGenomeAnnotation <- hub[[id]]

# Non-biased GSEA ---------------------------------------------------------

# GO Gene Set Enrichment Analysis - 1 vs 2

set.seed(424242)

geneList_1vs2 <- data_1vs2$log2FoldChange
names(geneList_1vs2) = as.character(data_1vs2$id)
geneList_1vs2 <- sort(geneList_1vs2, decreasing = TRUE)

GO_GSEA_1vs2 <- gseGO(geneList     = geneList_1vs2,
              OrgDb        = miceGenomeAnnotation,
              ont          = "BP",
              minGSSize    = 15,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = TRUE,
              keyType = 'ENSEMBL')

GO_GSEA_1vs2_df <- data.frame(GO_GSEA_1vs2) %>% 
  mutate('-log10(padj)' = -log10(p.adjust))
GO_GSEA_1vs2_df <- GO_GSEA_1vs2_df %>% 
  dplyr::arrange(`-log10(padj)`) %>% 
  filter(p.adjust < 0.05)

# writexl::write_xlsx(data.frame(setReadable(GO_GSEA_1vs2, OrgDb = miceGenomeAnnotation)), 'GO_GSEA_1vs2_df.xlsx')

GO_GSEA_1vs2_df_plot <- ggplot(GO_GSEA_1vs2_df, aes(factor(Description, levels = Description), `-log10(padj)`, fill = enrichmentScore)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), position = position_dodge(width = 0.9), hjust = 0.0001) +
  coord_flip() +
  ggtitle('GO_GSEA 1vs2')

# GO ORA - 1 vs 2

set.seed(424242)

genes_1vs2_grepl <- data_1vs2 %>% dplyr::filter(base::grepl('CCL|CXCL|XCL|CX3C', toupper(gene_ID)))
genes_1vs2_filter <- data_1vs2 %>% dplyr::filter(toupper(gene_ID) %in% c('TLR2','TLR4','IFN','TLR7','TLR9','IL12','IL6','TNF',
                                                                         'AZGP1',
                                                                         'B2M',
                                                                         'CALR',
                                                                         'CANX',
                                                                         'CD1A',
                                                                         'CD1B',
                                                                         'CD1C',
                                                                         'CD1D',
                                                                         'CD1E',
                                                                         'CD4',
                                                                         'CD8A',
                                                                         'CD8B',
                                                                         'CD74',
                                                                         'CREB1',
                                                                         'CTSB',
                                                                         'CTSE',
                                                                         'CTSL1',
                                                                         'CTSS',
                                                                         'FCER1G',
                                                                         'FCGRT',
                                                                         'PDIA3',
                                                                         'HFE',
                                                                         'HLA-A',
                                                                         'HLA-B',
                                                                         'HLA-C',
                                                                         'HLA-DMA',
                                                                         'HLA-DMB',
                                                                         'HLA-DOA',
                                                                         'HLA-DOB',
                                                                         'HLA-DPA1',
                                                                         'HLA-DPB1',
                                                                         'HLA-DQA1',
                                                                         'HLA-DQA2',
                                                                         'HLA-DQB1',
                                                                         'HLA-DRA',
                                                                         'HLA-DRB1',
                                                                         'HLA-DRB3',
                                                                         'HLA-DRB4',
                                                                         'HLA-DRB5',
                                                                         'HLA-E',
                                                                         'HLA-F',
                                                                         'HLA-G',
                                                                         'HLA-H',
                                                                         'MR1',
                                                                         'HSPA1A',
                                                                         'HSPA1B',
                                                                         'HSPA1L',
                                                                         'HSPA2',
                                                                         'HSPA4',
                                                                         'HSPA5',
                                                                         'HSPA6',
                                                                         'HSPA8',
                                                                         'HSP90AA1',
                                                                         'HSP90AB1',
                                                                         'ICAM1',
                                                                         'IFNA1',
                                                                         'IFNA2',
                                                                         'IFNA4',
                                                                         'IFNA5',
                                                                         'IFNA6',
                                                                         'IFNA7',
                                                                         'IFNA8',
                                                                         'IFNA10',
                                                                         'IFNA13',
                                                                         'IFNA14',
                                                                         'IFNA16',
                                                                         'IFNA17',
                                                                         'IFNA21',
                                                                         'IFNG',
                                                                         'KIR2DL1',
                                                                         'KIR2DL2',
                                                                         'KIR2DL3',
                                                                         'KIR2DL4',
                                                                         'KIR2DS1',
                                                                         'KIR2DS3',
                                                                         'KIR2DS4',
                                                                         'KIR2DS5',
                                                                         'KIR3DL1',
                                                                         'KIR3DL2',
                                                                         'KLRC1',
                                                                         'KLRC2',
                                                                         'KLRC3',
                                                                         'KLRD1',
                                                                         'LTA',
                                                                         'CIITA',
                                                                         'MICA',
                                                                         'MICB',
                                                                         'NFYA',
                                                                         'NFYB',
                                                                         'NFYC',
                                                                         'LGMN',
                                                                         'PSMB8',
                                                                         'PSMC1',
                                                                         'PSMC2',
                                                                         'PSMC3',
                                                                         'PSMC4',
                                                                         'PSMC5',
                                                                         'PSMC6',
                                                                         'PSMD1',
                                                                         'PSMD2',
                                                                         'PSMD3',
                                                                         'PSMD4',
                                                                         'PSMD5',
                                                                         'PSMD7',
                                                                         'PSMD8',
                                                                         'PSMD10',
                                                                         'PSMD11',
                                                                         'PSMD13',
                                                                         'PSME1',
                                                                         'PSME1',
                                                                         'PSME2',
                                                                         'PSME2',
                                                                         'RELB',
                                                                         'RFX5',
                                                                         'RFXAP',
                                                                         'SLC10A2',
                                                                         'TAP1',
                                                                         'TAP2',
                                                                         'TAPBP',
                                                                         'THBS1',
                                                                         'SHFM1',
                                                                         'KLRC4',
                                                                         'AP3B1',
                                                                         'RFXANK',
                                                                         'PSMD6',
                                                                         'PSME3',
                                                                         'PSMD14',
                                                                         'CLEC4M',
                                                                         'IFI30',
                                                                         'PROCR',
                                                                         'ADRM1',
                                                                         'KIAA0368',
                                                                         'TRPC4AP',
                                                                         'CD209',
                                                                         'UBXN1',
                                                                         'ERAP1',
                                                                         'TAPBPL',
                                                                         'KIR2DL5A',
                                                                         'ERAP2',
                                                                         'ULBP3',
                                                                         'ULBP2',
                                                                         'ULBP1',
                                                                         'KIR3DL3',
                                                                         'RAET1E',
                                                                         'RAET1L',
                                                                         'UBR1',
                                                                         'RAET1G',
                                                                         'PDIA2'))
genes_1vs2_all <- genes_1vs2_filter %>% dplyr::full_join(genes_1vs2_grepl)

geneList_1vs2_ORA <- genes_1vs2_all$log2FoldChange
names(geneList_1vs2_ORA) = as.character(genes_1vs2_all$id)
geneList_1vs2_ORA <- sort(geneList_1vs2_ORA, decreasing = TRUE)
gene_1vs2_ORA <- names(geneList_1vs2_ORA)[abs(geneList_1vs2_ORA) > 0.5]

GO_ORA_1vs2 <- enrichGO(gene          = gene_1vs2_ORA,
                universe      = names(geneList_1vs2_ORA),
                OrgDb         = miceGenomeAnnotation,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE,
                keyType = 'ENSEMBL')

GO_ORA_1vs2_df <- data.frame(GO_ORA_1vs2) %>% 
  mutate('-log10(padj)' = -log10(p.adjust))

GO_ORA_1vs2_df <- GO_ORA_1vs2_df %>% 
  dplyr::arrange(`-log10(padj)`) %>% 
  filter(p.adjust < 0.1)

# writexl::write_xlsx(data.frame(setReadable(GO_GSEA_1vs2, OrgDb = miceGenomeAnnotation)), 'GO_GSEA_1vs2_df.xlsx')

GO_ORA_1vs2_df_plot <- ggplot(GO_ORA_1vs2_df, aes(factor(Description, levels = Description), `-log10(padj)`)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  coord_flip() +
  ggtitle('GO_GSEA 1vs2')


# writexl::write_xlsx(data.frame(GO_ORA_1vs2), 'ORA_dendritic-antigenPresentation/GO_ORA_1vs2_dendritic-antigenPresentation.xlsx')

# KEGG 1vs2 ---------------------------------------------------------------

data_1vs2_kegg_mm <- AnnotationDbi::select(org.Mm.eg.db,
                                        key=data_1vs2$id, 
                                        columns="ENTREZID",
                                        keytype="ENSEMBL") %>% 
  dplyr::rename('id' = ENSEMBL) %>% 
  dplyr::full_join(data_1vs2) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::filter(!is.na(ENTREZID)) %>% 
  dplyr::distinct(ENTREZID, .keep_all = TRUE)


geneList_1vs2_kegg_mm <- data_1vs2_kegg_mm$log2FoldChange
names(geneList_1vs2_kegg_mm) = as.character(data_1vs2_kegg_mm$ENTREZID)
geneList_1vs2_kegg_mm <- sort(geneList_1vs2_kegg_mm, decreasing = TRUE)

set.seed(424242)

kk2_1vs2 <- gseKEGG(geneList     = geneList_1vs2_kegg_mm,
                    organism     = 'mmu',
                    #nPerm        = 10000,
                    maxGSSize    = 500,
                    minGSSize    = 10,
                    pvalueCutoff = 1,
                    verbose      = FALSE)

library("pathview")
kk2_1vs2_nk_cytotoxicity <- pathview(gene.data  = geneList_1vs2_kegg_mm,
                     pathway.id = "mmu04650",
                     species    = "mmu",
                     #limit      = list(gene=max(abs(geneList_1vs2_kegg_mm)), cpd=1)
                     )

# GO Gene Set Enrichment Analysis - 1 vs 3

set.seed(424242)

geneList_1vs3 <- data_1vs3$log2FoldChange
names(geneList_1vs3) = as.character(data_1vs3$id)
geneList_1vs3 <- sort(geneList_1vs3, decreasing = TRUE)

GO_GSEA_1vs3 <- gseGO(geneList     = geneList_1vs3,
                      OrgDb        = miceGenomeAnnotation,
                      ont          = "BP",
                      minGSSize    = 15,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = TRUE,
                      keyType = 'ENSEMBL')

GO_GSEA_1vs3_df <- data.frame(GO_GSEA_1vs3) %>% 
  mutate('-log10(padj)' = -log10(p.adjust))
GO_GSEA_1vs3_df <- GO_GSEA_1vs3_df %>% 
  dplyr::arrange(`-log10(padj)`) %>% 
  filter(p.adjust < 0.05)


# writexl::write_xlsx(data.frame(setReadable(GO_GSEA_1vs3, OrgDb = miceGenomeAnnotation)), 'GO_GSEA_1vs3_df.xlsx')

GO_GSEA_1vs3_df_plot <- ggplot(GO_GSEA_1vs3_df, aes(factor(Description, levels = Description), `-log10(padj)`, fill = enrichmentScore)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), position = position_dodge(width = 0.9), hjust = 0.0001) +
  coord_flip() +
  ggtitle('GO_GSEA 1vs3')

ggplot(GO_GSEA_1vs3_df, aes(factor(Description, levels = Description), `-log10(padj)`, fill = NES)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), position = position_dodge(width = 0.9), hjust = 0.0001) +
  coord_flip() +
  ggtitle('GO_GSEA 1vs3')

# GO ORA - 1 vs 3

set.seed(424242)

genes_1vs3_grepl <- data_1vs3 %>% dplyr::filter(base::grepl('CCL|CXCL|XCL|CX3C', toupper(gene_ID)))
genes_1vs3_filter <- data_1vs3 %>% dplyr::filter(toupper(gene_ID) %in% c('TLR2','TLR4','IFN','TLR7','TLR9','IL12','IL6','TNF',
                                                                         'AZGP1',
                                                                         'B2M',
                                                                         'CALR',
                                                                         'CANX',
                                                                         'CD1A',
                                                                         'CD1B',
                                                                         'CD1C',
                                                                         'CD1D',
                                                                         'CD1E',
                                                                         'CD4',
                                                                         'CD8A',
                                                                         'CD8B',
                                                                         'CD74',
                                                                         'CREB1',
                                                                         'CTSB',
                                                                         'CTSE',
                                                                         'CTSL1',
                                                                         'CTSS',
                                                                         'FCER1G',
                                                                         'FCGRT',
                                                                         'PDIA3',
                                                                         'HFE',
                                                                         'HLA-A',
                                                                         'HLA-B',
                                                                         'HLA-C',
                                                                         'HLA-DMA',
                                                                         'HLA-DMB',
                                                                         'HLA-DOA',
                                                                         'HLA-DOB',
                                                                         'HLA-DPA1',
                                                                         'HLA-DPB1',
                                                                         'HLA-DQA1',
                                                                         'HLA-DQA2',
                                                                         'HLA-DQB1',
                                                                         'HLA-DRA',
                                                                         'HLA-DRB1',
                                                                         'HLA-DRB3',
                                                                         'HLA-DRB4',
                                                                         'HLA-DRB5',
                                                                         'HLA-E',
                                                                         'HLA-F',
                                                                         'HLA-G',
                                                                         'HLA-H',
                                                                         'MR1',
                                                                         'HSPA1A',
                                                                         'HSPA1B',
                                                                         'HSPA1L',
                                                                         'HSPA2',
                                                                         'HSPA4',
                                                                         'HSPA5',
                                                                         'HSPA6',
                                                                         'HSPA8',
                                                                         'HSP90AA1',
                                                                         'HSP90AB1',
                                                                         'ICAM1',
                                                                         'IFNA1',
                                                                         'IFNA2',
                                                                         'IFNA4',
                                                                         'IFNA5',
                                                                         'IFNA6',
                                                                         'IFNA7',
                                                                         'IFNA8',
                                                                         'IFNA10',
                                                                         'IFNA13',
                                                                         'IFNA14',
                                                                         'IFNA16',
                                                                         'IFNA17',
                                                                         'IFNA21',
                                                                         'IFNG',
                                                                         'KIR2DL1',
                                                                         'KIR2DL2',
                                                                         'KIR2DL3',
                                                                         'KIR2DL4',
                                                                         'KIR2DS1',
                                                                         'KIR2DS3',
                                                                         'KIR2DS4',
                                                                         'KIR2DS5',
                                                                         'KIR3DL1',
                                                                         'KIR3DL2',
                                                                         'KLRC1',
                                                                         'KLRC2',
                                                                         'KLRC3',
                                                                         'KLRD1',
                                                                         'LTA',
                                                                         'CIITA',
                                                                         'MICA',
                                                                         'MICB',
                                                                         'NFYA',
                                                                         'NFYB',
                                                                         'NFYC',
                                                                         'LGMN',
                                                                         'PSMB8',
                                                                         'PSMC1',
                                                                         'PSMC2',
                                                                         'PSMC3',
                                                                         'PSMC4',
                                                                         'PSMC5',
                                                                         'PSMC6',
                                                                         'PSMD1',
                                                                         'PSMD2',
                                                                         'PSMD3',
                                                                         'PSMD4',
                                                                         'PSMD5',
                                                                         'PSMD7',
                                                                         'PSMD8',
                                                                         'PSMD10',
                                                                         'PSMD11',
                                                                         'PSMD13',
                                                                         'PSME1',
                                                                         'PSME1',
                                                                         'PSME2',
                                                                         'PSME2',
                                                                         'RELB',
                                                                         'RFX5',
                                                                         'RFXAP',
                                                                         'SLC10A2',
                                                                         'TAP1',
                                                                         'TAP2',
                                                                         'TAPBP',
                                                                         'THBS1',
                                                                         'SHFM1',
                                                                         'KLRC4',
                                                                         'AP3B1',
                                                                         'RFXANK',
                                                                         'PSMD6',
                                                                         'PSME3',
                                                                         'PSMD14',
                                                                         'CLEC4M',
                                                                         'IFI30',
                                                                         'PROCR',
                                                                         'ADRM1',
                                                                         'KIAA0368',
                                                                         'TRPC4AP',
                                                                         'CD209',
                                                                         'UBXN1',
                                                                         'ERAP1',
                                                                         'TAPBPL',
                                                                         'KIR2DL5A',
                                                                         'ERAP2',
                                                                         'ULBP3',
                                                                         'ULBP2',
                                                                         'ULBP1',
                                                                         'KIR3DL3',
                                                                         'RAET1E',
                                                                         'RAET1L',
                                                                         'UBR1',
                                                                         'RAET1G',
                                                                         'PDIA2'))
genes_1vs3_all <- genes_1vs3_filter %>% dplyr::full_join(genes_1vs3_grepl)

geneList_1vs3_ORA <- genes_1vs3_all$log2FoldChange
names(geneList_1vs3_ORA) = as.character(genes_1vs3_all$id)
geneList_1vs3_ORA <- sort(geneList_1vs3_ORA, decreasing = TRUE)
gene_1vs3_ORA <- names(geneList_1vs3_ORA)[abs(geneList_1vs3_ORA) > 0.5]

GO_ORA_1vs3 <- enrichGO(gene          = gene_1vs3_ORA,
                        universe      = names(geneList_1vs3_ORA),
                        OrgDb         = miceGenomeAnnotation,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE,
                        keyType = 'ENSEMBL')

GO_ORA_1vs3_df <- data.frame(GO_ORA_1vs3) %>% 
  mutate('-log10(padj)' = -log10(p.adjust))

GO_ORA_1vs3_df <- GO_ORA_1vs3_df %>% 
  dplyr::arrange(`-log10(padj)`) %>% 
  filter(p.adjust < 0.1)

# writexl::write_xlsx(data.frame(setReadable(GO_GSEA_1vs3, OrgDb = miceGenomeAnnotation)), 'GO_GSEA_1vs3_df.xlsx')

GO_ORA_1vs3_df_plot <- ggplot(GO_ORA_1vs3_df, aes(factor(Description, levels = Description), `-log10(padj)`)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  coord_flip() +
  ggtitle('GO_GSEA 1vs3')

# KEGG 1vs3 ---------------------------------------------------------------

data_1vs3_kegg_mm <- AnnotationDbi::select(org.Mm.eg.db,
                                           key=data_1vs3$id, 
                                           columns="ENTREZID",
                                           keytype="ENSEMBL") %>% 
  dplyr::rename('id' = ENSEMBL) %>% 
  dplyr::full_join(data_1vs3) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::filter(!is.na(ENTREZID)) %>% 
  dplyr::distinct(ENTREZID, .keep_all = TRUE)


geneList_1vs3_kegg_mm <- data_1vs3_kegg_mm$log2FoldChange
names(geneList_1vs3_kegg_mm) = as.character(data_1vs3_kegg_mm$ENTREZID)
geneList_1vs3_kegg_mm <- sort(geneList_1vs3_kegg_mm, decreasing = TRUE)

set.seed(424242)

kk2_1vs3 <- gseKEGG(geneList     = geneList_1vs3_kegg_mm,
                    organism     = 'mmu',
                    #nPerm        = 10000,
                    maxGSSize    = 500,
                    minGSSize    = 10,
                    pvalueCutoff = 1,
                    verbose      = FALSE)

library("pathview")
kk2_1vs2_akt_pathway <- pathview(gene.data  = geneList_1vs3_kegg_mm,
                     pathway.id = "mmu04151",
                     species    = "mmu",
                     #limit      = list(gene=max(abs(geneList_1vs3_kegg_mm)), cpd=1)
                     )

# GO Gene Set Enrichment Analysis - 2 vs 3

set.seed(424242)

geneList_2vs3 <- data_2vs3$log2FoldChange
names(geneList_2vs3) = as.character(data_2vs3$id)
geneList_2vs3 <- sort(geneList_2vs3, decreasing = TRUE)

GO_GSEA_2vs3 <- gseGO(geneList     = geneList_2vs3,
                      OrgDb        = miceGenomeAnnotation,
                      ont          = "BP",
                      minGSSize    = 15,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = TRUE,
                      keyType = 'ENSEMBL')

GO_GSEA_2vs3_df <- data.frame(GO_GSEA_2vs3) %>% 
  mutate('-log10(padj)' = -log10(p.adjust))
GO_GSEA_2vs3_df <- GO_GSEA_2vs3_df %>% 
  dplyr::arrange(`-log10(padj)`) %>% 
  filter(p.adjust < 0.05)


GO_GSEA_2vs3_df_plot <- ggplot(GO_GSEA_2vs3_df, aes(factor(Description, levels = Description), `-log10(padj)`, fill = enrichmentScore)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  geom_text(aes(label = round(enrichmentScore, 3)), position = position_dodge(width = 0.9), hjust = 0.0001) +
  coord_flip() +
  ggtitle('GO_GSEA 2vs3')

# GO ORA - 2 vs 3

set.seed(424242)

genes_2vs3_grepl <- data_2vs3 %>% dplyr::filter(base::grepl('CCL|CXCL|XCL|CX3C', toupper(gene_ID)))
genes_2vs3_filter <- data_2vs3 %>% dplyr::filter(toupper(gene_ID) %in% c('TLR2','TLR4','IFN','TLR7','TLR9','IL12','IL6','TNF',
                                                                         'AZGP1',
                                                                         'B2M',
                                                                         'CALR',
                                                                         'CANX',
                                                                         'CD1A',
                                                                         'CD1B',
                                                                         'CD1C',
                                                                         'CD1D',
                                                                         'CD1E',
                                                                         'CD4',
                                                                         'CD8A',
                                                                         'CD8B',
                                                                         'CD74',
                                                                         'CREB1',
                                                                         'CTSB',
                                                                         'CTSE',
                                                                         'CTSL1',
                                                                         'CTSS',
                                                                         'FCER1G',
                                                                         'FCGRT',
                                                                         'PDIA3',
                                                                         'HFE',
                                                                         'HLA-A',
                                                                         'HLA-B',
                                                                         'HLA-C',
                                                                         'HLA-DMA',
                                                                         'HLA-DMB',
                                                                         'HLA-DOA',
                                                                         'HLA-DOB',
                                                                         'HLA-DPA1',
                                                                         'HLA-DPB1',
                                                                         'HLA-DQA1',
                                                                         'HLA-DQA2',
                                                                         'HLA-DQB1',
                                                                         'HLA-DRA',
                                                                         'HLA-DRB1',
                                                                         'HLA-DRB3',
                                                                         'HLA-DRB4',
                                                                         'HLA-DRB5',
                                                                         'HLA-E',
                                                                         'HLA-F',
                                                                         'HLA-G',
                                                                         'HLA-H',
                                                                         'MR1',
                                                                         'HSPA1A',
                                                                         'HSPA1B',
                                                                         'HSPA1L',
                                                                         'HSPA2',
                                                                         'HSPA4',
                                                                         'HSPA5',
                                                                         'HSPA6',
                                                                         'HSPA8',
                                                                         'HSP90AA1',
                                                                         'HSP90AB1',
                                                                         'ICAM1',
                                                                         'IFNA1',
                                                                         'IFNA2',
                                                                         'IFNA4',
                                                                         'IFNA5',
                                                                         'IFNA6',
                                                                         'IFNA7',
                                                                         'IFNA8',
                                                                         'IFNA10',
                                                                         'IFNA13',
                                                                         'IFNA14',
                                                                         'IFNA16',
                                                                         'IFNA17',
                                                                         'IFNA21',
                                                                         'IFNG',
                                                                         'KIR2DL1',
                                                                         'KIR2DL2',
                                                                         'KIR2DL3',
                                                                         'KIR2DL4',
                                                                         'KIR2DS1',
                                                                         'KIR2DS3',
                                                                         'KIR2DS4',
                                                                         'KIR2DS5',
                                                                         'KIR3DL1',
                                                                         'KIR3DL2',
                                                                         'KLRC1',
                                                                         'KLRC2',
                                                                         'KLRC3',
                                                                         'KLRD1',
                                                                         'LTA',
                                                                         'CIITA',
                                                                         'MICA',
                                                                         'MICB',
                                                                         'NFYA',
                                                                         'NFYB',
                                                                         'NFYC',
                                                                         'LGMN',
                                                                         'PSMB8',
                                                                         'PSMC1',
                                                                         'PSMC2',
                                                                         'PSMC3',
                                                                         'PSMC4',
                                                                         'PSMC5',
                                                                         'PSMC6',
                                                                         'PSMD1',
                                                                         'PSMD2',
                                                                         'PSMD3',
                                                                         'PSMD4',
                                                                         'PSMD5',
                                                                         'PSMD7',
                                                                         'PSMD8',
                                                                         'PSMD10',
                                                                         'PSMD11',
                                                                         'PSMD13',
                                                                         'PSME1',
                                                                         'PSME1',
                                                                         'PSME2',
                                                                         'PSME2',
                                                                         'RELB',
                                                                         'RFX5',
                                                                         'RFXAP',
                                                                         'SLC10A2',
                                                                         'TAP1',
                                                                         'TAP2',
                                                                         'TAPBP',
                                                                         'THBS1',
                                                                         'SHFM1',
                                                                         'KLRC4',
                                                                         'AP3B1',
                                                                         'RFXANK',
                                                                         'PSMD6',
                                                                         'PSME3',
                                                                         'PSMD14',
                                                                         'CLEC4M',
                                                                         'IFI30',
                                                                         'PROCR',
                                                                         'ADRM1',
                                                                         'KIAA0368',
                                                                         'TRPC4AP',
                                                                         'CD209',
                                                                         'UBXN1',
                                                                         'ERAP1',
                                                                         'TAPBPL',
                                                                         'KIR2DL5A',
                                                                         'ERAP2',
                                                                         'ULBP3',
                                                                         'ULBP2',
                                                                         'ULBP1',
                                                                         'KIR3DL3',
                                                                         'RAET1E',
                                                                         'RAET1L',
                                                                         'UBR1',
                                                                         'RAET1G',
                                                                         'PDIA2'))
genes_2vs3_all <- genes_2vs3_filter %>% dplyr::full_join(genes_2vs3_grepl)

geneList_2vs3_ORA <- genes_2vs3_all$log2FoldChange
names(geneList_2vs3_ORA) = as.character(genes_2vs3_all$id)
geneList_2vs3_ORA <- sort(geneList_2vs3_ORA, decreasing = TRUE)
gene_2vs3_ORA <- names(geneList_2vs3_ORA)[abs(geneList_2vs3_ORA) > 0.5]

GO_ORA_2vs3 <- enrichGO(gene          = gene_2vs3_ORA,
                        universe      = names(geneList_2vs3_ORA),
                        OrgDb         = miceGenomeAnnotation,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 1,
                        qvalueCutoff  = 1,
                        readable      = TRUE,
                        keyType = 'ENSEMBL')

GO_ORA_2vs3_df <- data.frame(GO_ORA_2vs3) %>% 
  mutate('-log10(padj)' = -log10(p.adjust))

GO_ORA_2vs3_df <- GO_ORA_2vs3_df %>% 
  dplyr::arrange(`-log10(padj)`) %>% 
  filter(p.adjust < 0.1)

GO_ORA_2vs3_df_plot <- ggplot(GO_ORA_2vs3_df, aes(factor(Description, levels = Description), `-log10(padj)`)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
  coord_flip() +
  ggtitle('GO_GSEA 2vs3')


# writexl::write_xlsx(data.frame(GO_ORA_2vs3), 'ORA_dendritic-antigenPresentation/GO_ORA_2vs3_dendritic-antigenPresentation.xlsx')

# tight-junction mucin ----------------------------------------------------

geny_tj_muc <- c('CLDN11',
                 'CLDN2',
                 'CLDN24',
                 'CLDND2',
                 'CLDN17',
                 'CLDN9',
                 'CLDN10',
                 'JAM3',
                 'CLDN6',
                 'MAGIX',
                 'CLDN16',
                 'CLDN13',
                 'TMEM1',
                 'CLDN15',
                 'CLDN20',
                 'CLDN22',
                 'MAGI2',
                 'MARVELD2',
                 'CLDN23',
                 'MAGI3',
                 'MAGI1',
                 'CLDN8',
                 'MPDZ',
                 'CGN',
                 'F11R',
                 'CLDN3',
                 'TJP2',
                 'CLDND1',
                 'CLDN12',
                 'MARVELD3',
                 'INADL',
                 'JAM2',
                 'TJP1',
                 'CLDN18',
                 'MARVELD1',
                 'CLDN5',
                 'CGNL1',
                 'OCLN',
                 'CLDN19',
                 'TMEM235',
                 'TJAP1',
                 'CLDN14',
                 'CLDN4',
                 'CLMP',
                 'TJP3',
                 'CLDN1',
                 'CLDN7',
                 'CLCA1',
                 'PROL1',
                 'MUC3',
                 'OVGP1',
                 'MUC13',
                 'MUC19',
                 'MUC2',
                 'MUC15',
                 'MUC6',
                 'MUC4',
                 'MUC5AC',
                 'MUC5B',
                 'MUC16',
                 'MUC20',
                 'MUC1',
                 'EMCN'
)

df_tj_muc_1vs2 <- data_1vs2 %>% 
  dplyr::filter(toupper(gene_ID) %in% geny_tj_muc) %>% 
  dplyr::select(gene_ID, log2FoldChange, pvalue, padj) %>% 
  dplyr::rename('Log2FC_1vs2' = log2FoldChange,
                'pvalue_1vs2' = pvalue,
                'padj1vs2' = padj)
df_tj_muc_1vs3 <- data_1vs3 %>% 
  dplyr::filter(toupper(gene_ID) %in% geny_tj_muc) %>% 
  dplyr::select(gene_ID, log2FoldChange, pvalue, padj) %>% 
  dplyr::rename('Log2FC_1vs3' = log2FoldChange,
                'pvalue_1vs3' = pvalue,
                'padj1vs3' = padj)
df_tj_muc_2vs3 <- data_2vs3 %>% 
  dplyr::filter(toupper(gene_ID) %in% geny_tj_muc) %>% 
  dplyr::select(gene_ID, log2FoldChange, pvalue, padj) %>% 
  dplyr::rename('Log2FC_2vs3' = log2FoldChange,
                'pvalue_2vs3' = pvalue,
                'padj2vs3' = padj)

df_tj_muc_combined <- df_tj_muc_1vs2 %>% 
  dplyr::full_join(df_tj_muc_1vs3) %>% 
  dplyr::full_join(df_tj_muc_2vs3)

# writexl::write_xlsx(df_tj_muc_combined, 'tightJunction_mucin-combined.xlsx')

df_tj_muc_combined_heatmap <- df_tj_muc_combined %>% 
  dplyr::select(gene_ID, Log2FC_1vs2, Log2FC_1vs3, Log2FC_2vs3) %>% 
  tibble::column_to_rownames(var = 'gene_ID')


col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ComplexHeatmap::Heatmap(as.matrix(df_tj_muc_combined_heatmap), name = 'Log2FoldChange',
                        cluster_rows = T, cluster_columns = F, show_row_names = T,
                        border = T, row_dend_width = unit(4, "cm"), col = col_fun)

temp1 <- genes_1vs2_grepl %>% 
  dplyr::full_join(genes_1vs3_grepl) %>% 
  dplyr::full_join(genes_2vs3_grepl) %>%
  dplyr::mutate(gene_ID = toupper(gene_ID)) %>% 
  dplyr::select(gene_ID) %>% 
  dplyr::distinct(gene_ID, .keep_all = T)



# cztery_geny <- temp1$gene_ID

# inflammasome related from tina.tan
cztery_geny <- c('IL12A','NLRP5','NLRP9A','NLRPX1','NAIP','NLRP4B','IL12B','NRLP4D',
                 'NLRP1A','NLRP9B','NLRP4C','IFNB1','NLRP4F','BCL2','NLRP1B','NLRP4E',
                 'POP1','NLRP4A','PYCARD','ASC','INFG','TNFSF4','NLRC4','NLRP9C','AIM2',
                 'CARD6','SKP1A','HSP90AA1','UBE2N','CHUK','MAP2K4','PELI1','TAB2','NLRP3K7',
                 'MAPK1','MAPK8','IL18','MAPK2K6','P2RX7','PEA15A','TAB1','MAPK9','TAB3','CASP8',
                 'BTRC','MAPK12','IKBKG','XIAP','HSP90B1','SUGT1','MEFV','NLRP3','MAPK13','TNF',
                 'CD40LG','PEA15B','TNFSF11','PELI3','CASP1','CCL2','PANX1','PSTPIP1','CCL7',
                 'NOD2','IL6','CXCL2','IL1RN','IL1R2','TNFSF14','CTSB','NFKBIA','NFKB1','NFKBIB',
                 'MAP2K1','MAPK11','BCL2L1','BIRC2','CCL5','NIRC5','CIITA','CITA','IRF1','RIPK2',
                 'IFIH1','TNIP2','IL1A','BIRC3','IL1B','MYD88','CXCL1','IRAK3','CASP4','IL33',
                 'IL1R1','MAPK3K8','IRAK4','TIRAP','IRF2','DDX58','RIG1','RIGL','TRAF6','IRAK2',
                 'RELA','NLRP6','PELI2','IRAK1','TMEM189','APP','MAPK3','MAPK3K3','AGER','TXNIP',
                 'IL1RAP','NOD1','RBX1','CUL1','CFLAR','TXN','IKBKBK','SQSTM1','FADD','HSP90AB1',
                 'TOLLIP')

cztery_geny <- c('FASL', 'TNFSF10', 'TGFB', 'IL10','TGFB1', 'TGFB2', 'TGFB3', 'IFNB',
                 'PTGES2')

cztery_geny <- c('IFNG','IDO','PGE2','IL1A','IL1B','IL6','IL10',
                 'IL17A','IL22','PDL1','PDL2','PTGES2','TNFSF10','FASL','TGFB3',
                 'TGFB2', 'TGFB1','FOXP3', 'SMAD7')

cztery_geny <- c('PTGES2','IL6','TGFB1','TGFB2','TGFB3','CCL2','CSF3')

cztery_geny <- c('NLRP3','AIM2','RIG1','CASP1','ASC','PYCARD', 'TLR2','TLR4',
                 'CASP18')

# sars cov 2 receptors
cztery_geny <- c('ACE2',
'ACE',
'TMPRSS2',
'SLC6A19',
'APH1A',
'APH1B',
'APOD',
'BSG',
'CD44',
'EGFR',
'GP6',
'ITGA3',
'ITGA6',
'ITGB1',
'JUP',
'LGALS3',
'NCSTN',
'NME1',
'NOD2',
'NXNL1',
'PPIA',
'PPIB',
'PSEN1',
'PSENEN',
'S100A9',
'SDC1',
'SLC16A1',
'SLC16A3',
'SLC16A4',
'SLC16A7',
'SLC16A8',
'SLC2A1',
'SLC3A2',
'SLC7A5',
'SPN',
'NFAT5',
'NFATC1',
'NFATC2',
'NFATC3',
'NFATC4',
'DPP4')

#sars cov 2 anti viral
# doi:10.3390/v12091039
cztery_geny <- c('TLR3','TLR7','TLR8','IFNG','IFNZ','INFA','IL6','MIP1A',
                 'MCP3','CSF2','IL2','IL10','IP10','CCL2','MCP1','CXCL1','CXCL5','
                 IL1A', 'IL1B')

cztery_geny1vs2 <- data_1vs2 %>% 
  dplyr::filter(toupper(gene_ID) %in% cztery_geny) %>% 
  dplyr::select(gene_ID, log2FoldChange, pvalue, padj) %>% 
  dplyr::rename('Log2FC_1vs2' = log2FoldChange,
                'pvalue_1vs2' = pvalue,
                'padj1vs2' = padj)
cztery_geny1vs3 <- data_1vs3 %>% 
  dplyr::filter(toupper(gene_ID) %in% cztery_geny) %>% 
  dplyr::select(gene_ID, log2FoldChange, pvalue, padj) %>% 
  dplyr::rename('Log2FC_1vs3' = log2FoldChange,
                'pvalue_1vs3' = pvalue,
                'padj1vs3' = padj)
cztery_geny2vs3 <- data_2vs3 %>% 
  dplyr::filter(toupper(gene_ID) %in% cztery_geny) %>% 
  dplyr::select(gene_ID, log2FoldChange, pvalue, padj) %>% 
  dplyr::rename('Log2FC_2vs3' = log2FoldChange,
                'pvalue_2vs3' = pvalue,
                'padj2vs3' = padj)

cztery_genycombined <- cztery_geny1vs2 %>% 
  dplyr::full_join(cztery_geny1vs3) %>% 
  dplyr::full_join(cztery_geny2vs3)

# writexl::write_xlsx(df_tj_muc_combined, 'tightJunction_mucin-combined.xlsx')

cztery_genyheatmap <- cztery_genycombined %>% 
  dplyr::select(gene_ID, Log2FC_1vs2, Log2FC_1vs3, Log2FC_2vs3) %>% 
  mutate(gene_ID = toupper(gene_ID)) %>% 
  tibble::column_to_rownames(var = 'gene_ID')


col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ComplexHeatmap::Heatmap(as.matrix(cztery_genyheatmap), name = 'Log2FoldChange',
                        cluster_rows = F, cluster_columns = F, show_row_names = T,
                        border = T, row_dend_width = unit(4, "cm"), col = col_fun, 
                        column_title = 'ACE-2-, CD147- and CD26- related genes (receptors)')

# writexl::write_xlsx(cztery_genycombined, 'cztery_genycombined-combined.xlsx')


# 2vs3 Kegg ---------------------------------------------------------------

data_2vs3_kegg_mm <- AnnotationDbi::select(org.Mm.eg.db,
                                           key=data_2vs3$id, 
                                           columns="ENTREZID",
                                           keytype="ENSEMBL") %>% 
  dplyr::rename('id' = ENSEMBL) %>% 
  dplyr::full_join(data_2vs3) %>% 
  dplyr::arrange(padj) %>% 
  dplyr::filter(!is.na(ENTREZID)) %>% 
  dplyr::distinct(ENTREZID, .keep_all = TRUE)


geneList_2vs3_kegg_mm <- data_2vs3_kegg_mm$log2FoldChange
names(geneList_2vs3_kegg_mm) = as.character(data_2vs3_kegg_mm$ENTREZID)
geneList_2vs3_kegg_mm <- sort(geneList_2vs3_kegg_mm, decreasing = TRUE)

set.seed(424242)

kk2_2vs3 <- gseKEGG(geneList     = geneList_2vs3_kegg_mm,
                    organism     = 'mmu',
                    #nPerm        = 10000,
                    maxGSSize    = 500,
                    minGSSize    = 10,
                    pvalueCutoff = 1,
                    verbose      = FALSE)
