
https://www.biostars.org/p/467197/#484439
  
# setwd()  ----------------------------------------------------------------
setwd('D:\\!ZMRiI\\!MSC_PAPER')

# Dependencies ------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)

# Setup -------------------------------------------------------------------
data_1vs2 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\1vs2.xlsx')
data_1vs2$log2FoldChange <- as.numeric(data_1vs2$log2FoldChange)
data_1vs2$padj <- as.numeric(data_1vs2$padj)
data_1vs2$stat <- as.numeric(data_1vs2$stat)

data_1vs3 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\1vs3.xlsx')
data_1vs3$log2FoldChange <- as.numeric(data_1vs3$log2FoldChange)
data_1vs3$padj <- as.numeric(data_1vs3$padj)
data_1vs3$stat <- as.numeric(data_1vs3$stat)

data_2vs3 <- read_excel('D:\\!ZMRiI\\!MSC_PAPER\\2vs3.xlsx')
data_2vs3$log2FoldChange <- as.numeric(data_2vs3$log2FoldChange)
data_2vs3$padj <- as.numeric(data_2vs3$padj)
data_2vs3$stat <- as.numeric(data_2vs3$stat)
#-
data_1vs2_filtred <- data_1vs2 %>% 
  full_join(AnnotationDbi::select(org.Mm.eg.db,
                                  key=data_1vs2$id, 
                                  columns="SYMBOL",
                                  keytype="ENSEMBL") %>% 
              dplyr::rename('id' = ENSEMBL)) %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  dplyr::arrange(stat) %>% 
  # I don't think that we need distinct() here, as we have group_by() piped to summarize(mean())
  #distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))

data_1vs2_ranks <- data_1vs2_filtred$stat
names(data_1vs2_ranks) <- data_1vs2_filtred$SYMBOL

# Load --------------------------------------------------------------------
GO_GSEA_1vs2 <- gseGO(geneList     = data_1vs2_ranks,
                      OrgDb        = org.Mm.eg.db,
                      ont          = "BP",
                      minGSSize    = 15,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      
                      verbose      = TRUE,
                      keyType = 'SYMBOL')
