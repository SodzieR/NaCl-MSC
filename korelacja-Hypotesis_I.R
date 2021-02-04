
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(readxl)
# rstatix is for kendalls corr matrix
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

data_1_cor_r_corrplot <- data_1_cor_r %>%
  column_to_rownames(var = 'rowname')

data_1_cor_p_corrplot <- data_1_cor_p %>%
  column_to_rownames(var = 'rowname')

#
#pdf("group1.pdf", width=7, height=7)
corrplot(as.matrix(data_1_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_1_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'Significant data 1', method = 'square', tl.srt = 45)
##dev.off()
#
#
#pdf("group1-all.pdf", width=7, height=7)
corrplot(as.matrix(data_1_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_1_cor_p_corrplot), sig.level = 0.05, insig = "pch", 
         title = 'Significant data 1', method = 'square', addrect = 3, tl.srt = 45)
##dev.off()
#
#pdf("group1-all-circular.pdf", width=7, height=7)
col_fun1 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
circos.heatmap(data_1_cor_r_corrplot, col = col_fun1, dend.side = "inside", 
               rownames.side = "outside", rownames.cex = 1)
lgd = ComplexHeatmap::Legend(title = "Group 1 - r coefficient", col_fun = col_fun1)
grid::grid.draw(lgd)
circos.clear()
##dev.off()

#
#pdf("group1-chord-corr.pdf", width=7, height=7)
set.seed(424242)
chordDiagram(data_1_cor_r_corrplot, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.01))
#dev.off()
#
#pdf("group1-chord-data.pdf", width=14, height=14)
set.seed(424242)
chordDiagram(data_1)
#dev.off()
#
data_1_cor_r_circular_significant <- data_1_cor_r %>% 
  cor_gather() %>% 
  dplyr::filter(p < 0.05)
rownames(data_1_cor_r_circular_significant) <- NULL

data_1_cor_r_circular_significant_matrix <- cor_spread(data_1_cor_r_circular_significant)

data_1_cor_r_circular_significant_matrix_noRownames <- data_1_cor_r_circular_significant_matrix %>%
  column_to_rownames(var = 'rowname')

data_1_cor_r_circular_significant_matrix_noRownames[is.na(data_1_cor_r_circular_significant_matrix_noRownames)] <- 0

#pdf("group1_corr_circular_singificant.pdf", width=7, height=7)
col_fun1 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
circos.heatmap(as.matrix(data_1_cor_r_circular_significant_matrix_noRownames), col = col_fun1, dend.side = "inside", 
               rownames.side = "outside", rownames.cex = 1, bg.border = "black")
lgd = ComplexHeatmap::Legend(title = "Group 1 - r coefficient", col_fun = col_fun1)
grid::grid.draw(lgd)
circos.clear()
#dev.off()
#
#pdf("group1-chord-significant.pdf", width=7, height=7)
chordDiagram(data_1_cor_r_circular_significant_matrix_noRownames)
#dev.off()

##########################################
mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]
data_1_cor_r_circular_significant_matrix_noRownames_sorted <- mat_ord(data_1_cor_r_circular_significant_matrix_noRownames)

diag(data_1_cor_r_circular_significant_matrix_noRownames_sorted) <- 0

data_1_cor_r_circular_significant_matrix_noRownames_sorted1 <- data_1_cor_r_circular_significant_matrix_noRownames_sorted %>% 
  arrange(rownames(data_1_cor_r_circular_significant_matrix_noRownames_sorted))

#colfunc <- colorRamps::magenta1green(18)

#colfunc <- colorRamps::primary.colors(18)

# get grid colors 
othercol1 = structure(rep("blue", length(data_1_cor_r_circular_significant_matrix_noRownames_sorted1)),
                     names = rownames(data_1_cor_r_circular_significant_matrix_noRownames_sorted1)) 
colfunc = c(othercol1)

#link colors
link_colors = colorRamp2(range(data_1_cor_r_circular_significant_matrix_noRownames_sorted1), c("red", "blue"), transparency = 0.25)

#pdf("SAMPLE-GROUP1-chord-significant.pdf", width=7, height=7)
colfunc = othercol1

col_fun_cor_2 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
lgd_links = ComplexHeatmap::Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = othercol, 
                                   title_position = "topleft", title = "Links")

chordDiagram(as.matrix(data_1_cor_r_circular_significant_matrix_noRownames_sorted1), grid.col = colfunc,
             annotationTrack = c("name", "grid"), col = link_colors)
title("Group 1 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
#dev.off()
# Data_2 ------------------------------------------------------------------

data_2_cor_r <- data.frame(data_2) %>% 
  cor_mat(method = 'kendall')

rownames(data_2_cor_r) <- NULL

data_2_cor_p <- data_2_cor_r %>% 
  cor_get_pval()

rownames(data_2_cor_p) <- NULL

data_2_cor_r_corrplot <- data_2_cor_r %>%
  column_to_rownames(var = 'rowname')

data_2_cor_p_corrplot <- data_2_cor_p %>%
  column_to_rownames(var = 'rowname')

#
#pdf("group2.pdf", width=7, height=7)
corrplot(as.matrix(data_2_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_2_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'Significant data 2', method = 'square', tl.srt = 45)
#dev.off()
#
#data_2_cor_r %>%
#  cor_reorder() %>%
#  pull_lower_triangle() %>%
#  cor_plot(label = TRUE)
#
#pdf("group2-all.pdf", width=7, height=7)
corrplot(as.matrix(data_2_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_2_cor_p_corrplot), sig.level = 0.05, insig = "pch", 
         title = 'Significant data 2', method = 'square', addrect = 3, tl.srt = 45)
#dev.off()
#

#pdf("group2-all-circular.pdf", width=7, height=7)
col_fun1 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
circos.heatmap(data_2_cor_r_corrplot, col = col_fun1, dend.side = "inside", 
               rownames.side = "outside", rownames.cex = 1)
lgd = ComplexHeatmap::Legend(title = "Group 2 - r coefficient", col_fun = col_fun1)
grid::grid.draw(lgd)
circos.clear()
#dev.off()
#
#pdf("group2-chord-corr.pdf", width=7, height=7)
set.seed(424242)
chordDiagram(data_2_cor_r_corrplot, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.01))
#dev.off()
#
#pdf("group2-chord-data.pdf", width=14, height=14)
set.seed(424242)
chordDiagram(data_2)
dev.off

#
data_2_cor_r_circular_significant <- data_2_cor_r %>% 
  cor_gather() %>% 
  dplyr::filter(p < 0.05)
rownames(data_2_cor_r_circular_significant) <- NULL

data_2_cor_r_circular_significant_matrix <- cor_spread(data_2_cor_r_circular_significant)

data_2_cor_r_circular_significant_matrix_noRownames <- data_2_cor_r_circular_significant_matrix %>%
  column_to_rownames(var = 'rowname')

data_2_cor_r_circular_significant_matrix_noRownames[is.na(data_2_cor_r_circular_significant_matrix_noRownames)] <- 0

#pdf("group2_corr_circular_singificant.pdf", width=7, height=7)
col_fun1 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
circos.heatmap(as.matrix(data_2_cor_r_circular_significant_matrix_noRownames), col = col_fun1, dend.side = "inside", 
               rownames.side = "outside", rownames.cex = 1, bg.border = "black")
lgd = ComplexHeatmap::Legend(title = "Group 2 - r coefficient", col_fun = col_fun1)
grid::grid.draw(lgd)
circos.clear()
#dev.off()
#
#pdf("group2-chord-significant.pdf", width=7, height=7)
chordDiagram(data_2_cor_r_circular_significant_matrix_noRownames)
#dev.off()

mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]
data_3_cor_r_circular_significant_matrix_noRownames_sorted <- mat_ord(data_3_cor_r_circular_significant_matrix_noRownames)

diag(data_3_cor_r_circular_significant_matrix_noRownames_sorted) <- 0

data_3_cor_r_circular_significant_matrix_noRownames_sorted2 <- data_3_cor_r_circular_significant_matrix_noRownames_sorted %>% 
  arrange(rownames(data_3_cor_r_circular_significant_matrix_noRownames_sorted))

#colfunc <- colorRamps::magenta2green(28)

#colfunc <- colorRamps::primary.colors(28)

# get grid colors 

#link colors
link_colors = colorRamp2(range(data_1_cor_r_circular_significant_matrix_noRownames_sorted1), c("red", "blue"), transparency = 0.25)

#pdf("TEMPTEMPTEMP-chord-significant.pdf", width=7, height=7)

col_fun_cor_3 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
lgd_links = ComplexHeatmap::Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_cor_3, 
                                   title_position = "topleft", title = "Link")

chordDiagram(as.matrix(data_1_cor_r_circular_significant_matrix_noRownames_sorted1), grid.col = colfunc,
             annotationTrack = c("name", "grid"), col = link_colors)
title("Group 3 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links, x = unit(3, "cm"), just = 'left')
#dev.off()

##################################################################################################
mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]
data_2_cor_r_circular_significant_matrix_noRownames_sorted <- mat_ord(data_2_cor_r_circular_significant_matrix_noRownames)

diag(data_2_cor_r_circular_significant_matrix_noRownames_sorted) <- 0

data_2_cor_r_circular_significant_matrix_noRownames_sorted2 <- data_2_cor_r_circular_significant_matrix_noRownames_sorted %>% 
  arrange(rownames(data_2_cor_r_circular_significant_matrix_noRownames_sorted))

#colfunc <- colorRamps::magenta2green(28)

#colfunc <- colorRamps::primary.colors(28)

# get grid colors 
othercol = structure(rep("blue", length(data_2_cor_r_circular_significant_matrix_noRownames_sorted2)),
                     names = rownames(data_2_cor_r_circular_significant_matrix_noRownames_sorted2)) 
remove <- c('IL10','IL15RA' ,'TGFB1','CCL2')  
othercol2 <- othercol[!othercol%in% remove]
colfunc = c('IL10' = 'red','IL15RA' = 'red','TGFB1' = 'red','CCL2' = 'red', othercol2)


#link colors
link_colors = colorRamp2(range(data_3_cor_r_circular_significant_matrix_noRownames_sorted2), c("red", "blue"), transparency = 0.25)

#pdf("SAMPLE-GROUP2-chord-significant.pdf", width=7, height=7)

col_fun_cor_3 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
lgd_links = ComplexHeatmap::Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_cor_3, 
                                   title_position = "topleft", title = "Link")

chordDiagram(as.matrix(data_2_cor_r_circular_significant_matrix_noRownames_sorted2), grid.col = colfunc,
             annotationTrack = c("name", "grid"), col = link_colors)
title("Group 2 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

#dev.off()
# Data_3 ------------------------------------------------------------------

data_3_cor_r <- data.frame(data_3) %>% 
  cor_mat(method = 'kendall')

rownames(data_3_cor_r) <- NULL

data_3_cor_p <- data_3_cor_r %>% 
  cor_get_pval()

rownames(data_3_cor_p) <- NULL

data_3_cor_r_corrplot <- data_3_cor_r %>%
  column_to_rownames(var = 'rowname')

data_3_cor_p_corrplot <- data_3_cor_p %>%
  column_to_rownames(var = 'rowname')

#
#pdf("group3.pdf", width=7, height=7)
corrplot(as.matrix(data_3_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_3_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'Significant data 3', method = 'square', tl.srt = 45)
#dev.off()
#
#data_3_cor_r %>%
#  cor_reorder() %>%
#  pull_lower_triangle() %>%
#  cor_plot(label = TRUE)

#
#pdf("group3-all.pdf", width=7, height=7)
corrplot(as.matrix(data_3_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_3_cor_p_corrplot), sig.level = 0.05, insig = "pch", 
         title = 'Significant data 3', method = 'square', addrect = 3, tl.srt = 45)
#dev.off()
#
# ---
# https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/
# ---
#pdf("group3-all-circular.pdf", width=7, height=7)
col_fun1 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
circos.heatmap(data_3_cor_r_corrplot, col = col_fun1, dend.side = "inside", 
               rownames.side = "outside", rownames.cex = 1)
lgd = ComplexHeatmap::Legend(title = "Group 3 - r coefficient", col_fun = col_fun1)
grid::grid.draw(lgd)
circos.clear()
#dev.off()

#
#pdf("group3-chord-corr.pdf", width=7, height=7)
set.seed(424242)
chordDiagram(data_3_cor_r_corrplot, annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.01))
#dev.off()
#
#pdf("group3-chord-data.pdf", width=14, height=14)
set.seed(424242)
chordDiagram(data_3)
#dev.off()
#
data_3_cor_r_circular_significant <- data_3_cor_r %>% 
  cor_gather() %>% 
  dplyr::filter(p < 0.05)
rownames(data_3_cor_r_circular_significant) <- NULL

data_3_cor_r_circular_significant_matrix <- cor_spread(data_3_cor_r_circular_significant)

data_3_cor_r_circular_significant_matrix_noRownames <- data_3_cor_r_circular_significant_matrix %>%
  column_to_rownames(var = 'rowname')

data_3_cor_r_circular_significant_matrix_noRownames[is.na(data_3_cor_r_circular_significant_matrix_noRownames)] <- 0

#pdf("group3_corr_circular_singificant.pdf", width=7, height=7)
col_fun1 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
circos.heatmap(as.matrix(data_3_cor_r_circular_significant_matrix_noRownames), col = col_fun1, dend.side = "inside", 
               rownames.side = "outside", rownames.cex = 1, bg.border = "black")
lgd = ComplexHeatmap::Legend(title = "Group 3 - r coefficient", col_fun = col_fun1)
grid::grid.draw(lgd)
circos.clear()
#dev.off()
#
#pdf("group3-chord-significant.#pdf", width=7, height=7)

#############################################################################
mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]
data_3_cor_r_circular_significant_matrix_noRownames_sorted <- mat_ord(data_3_cor_r_circular_significant_matrix_noRownames)

diag(data_3_cor_r_circular_significant_matrix_noRownames_sorted) <- 0

data_3_cor_r_circular_significant_matrix_noRownames_sorted2 <- data_3_cor_r_circular_significant_matrix_noRownames_sorted %>% 
  arrange(rownames(data_3_cor_r_circular_significant_matrix_noRownames_sorted))

#colfunc <- colorRamps::magenta2green(28)

#colfunc <- colorRamps::primary.colors(28)

# get grid colors 
othercol = structure(rep("blue", length(data_3_cor_r_circular_significant_matrix_noRownames_sorted2)),
                     names = rownames(data_3_cor_r_circular_significant_matrix_noRownames_sorted2)) 
remove <- c('JAK3','IL5')  
othercol2 <- othercol[!othercol%in% remove]
colfunc = c("IL5" = "red", "JAK3" = "red", othercol2)

#link colors
link_colors = colorRamp2(range(data_3_cor_r_circular_significant_matrix_noRownames_sorted2), c("red", "blue"), transparency = 0.25)

#pdf("SAMPLE-GROUP3-chord-significant.pdf", width=7, height=7)

col_fun_cor_3 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
lgd_links = ComplexHeatmap::Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_cor_3, 
                   title_position = "topleft", title = "Link")

chordDiagram(as.matrix(data_3_cor_r_circular_significant_matrix_noRownames_sorted2), grid.col = colfunc,
             annotationTrack = c("name", "grid"), col = link_colors)
title("Group 3 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
#dev.off()
# -------------------------------------------------------------------------
# ---
# https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html
##dev.off()
# ---

# Connection plot ---------------------------------------------------------

#===============================================================================================
#
# grey2red
#
#===============================================================================================

grey2red = function(n, base, gamma)
{
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = red = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, red, maxColorValue = 1); 
}


# Example of grey2red:

if (TRUE)
{
  par(mfrow = c(5,1))
  par(mar = c(1,3,1,1))
  n= 100
  barplot(rep(1, n), col = grey2red(n, 0, 1))
  barplot(rep(1, n), col = grey2red(n, 1, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 1))
  barplot(rep(1, n), col = grey2red(n, 0.5, 0.2))
  barplot(rep(1, n), col = grey2red(n, 0.5, 5.0))
}



#===============================================================================================
#
# Circle plot for generating VisANT-like plots
#
#===============================================================================================

circlePlot = function(
  adjacency,
  labels,
  order,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0), 
  radii = c(0.8, 0.8),
  startAngle = 0,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  lineColors = grey2red(50, 0.6, 1),
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
  xMargin = 1-radii[1],
  yMargin = 1-radii[2],
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  
  ...)
{
  
  if (startNewPlot)
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);
  
  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  #checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);
  
  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = apply(abs(adjx), 2, sum)-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);
  
  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);
  
  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];
  
  actualCexPts = rep(0, n);
  for (node in 1:n)
  {
    cex = min.cex.points;
    if (variable.cex.points)
      cex = min.cex.points + (max.cex.points - min.cex.points) * 
        (oConn[node] - minConn)/(maxConn - minConn)
    actualCexPts[node] = cex
  }
  
  diag(oAdj) = 0;
  maxA = max(abs(oAdj));
  if (sum(oAdj < 0) > 0)
  {
    adjCol = numbers2colors(oAdj, signed = TRUE, lim = c(-maxA, maxA), colorspace::diverge_hsv(2));
    
  } else {
    adjCol = numbers2colors(oAdj, signed = FALSE, lim = c(0, maxA), colorspace::diverge_hsv(2));
    
  }
  
  
  ltA = oAdj;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;
  
  adjOrder = order(c(abs(ltA)))
  rows = row(oAdj)[adjOrder];
  cols = col(oAdj)[adjOrder];
  
  nLines = n*(n-1)/2;
  for (line in 1:nLines)
  {
    n1 = rows[line];
    n2 = cols[line];
    a = oAdj[n1, n2];
    normA = abs(a)/maxA;
    
    w = min.line.width;
    if (variable.line.width)
      w = min.line.width + (max.line.width - min.line.width) * normA;
    
    #pRadius1 = par("cxy") * actualCexPts[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * actualCexPts[n2]/35;
    lineLen = sqrt( (x[n1] - x[n2])^2 + (y[n1] - y[n2])^2);
    x1 = x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen
    
    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n2, n1]);
  }
  
  for (node in 1:n)
    points(x[node], y[node], pch = oPch[node], cex = actualCexPts[node], bg = oPBg[node], col = oPColors[node]);
  
  for (node in 1:n)
  {
    cex = min.cex.labels;
    if (variable.cex.labels)
      cex = min.cex.labels + (max.cex.labels - min.cex.labels) *
        (oConn[node] - minConn)/(maxConn - minConn)
    textWidth = strwidth(oLabs[node], cex = cex);
    textHeight = strheight(oLabs[node], cex = cex);
    if (variableLabelAngle)
    {
      ang = angles[node]/pi * 180;
      if (ang < 180) 
      {
        dir = 1;
      } else {
        dir = -1;
        ang = ang - 180;
      }
      ang = (90 - ang)/2
      xDir = 1;
      yDir = 1;
      cosAng = cos(ang/180*pi);
      sinAng = sin(ang/180*pi);
    } else {
      ang = 0;
      xDir = x[node];
      yDir = y[node];
      cosAng = 1;
      sinAng = 1;
      dir = 1;
    }
    angRad = ang/180*pi;
    pRadius = par("cxy") * actualCexPts[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(sum(c(cosAng^2, sinAng^2) * pRadius^2));
    rotMat = matrix( c(cosAng, sinAng, -sinAng, cosAng), 2, 2);
    labelShift = rotMat %*% as.matrix(c(textWidth, textHeight));
    text(x[node] + dir * xDir * (labelShift[1]/2 + cosAng * effPointRadius + xLabelOffset[node]), 
         y[node] + dir * yDir * (labelShift[2]/2 + sinAng * effPointRadius + yLabelOffset[node]), 
         labels = oLabs[node], adj = c(0.5, 0.5), 
         cex = cex, col = oLColors[node], srt = ang, xpd = TRUE);
  }
  
}

# Example of circle plot:

if (FALSE)
{
  
  sizeGrWindow(8,8)
  par(mfrow = c(1,1));
  nS = 100;
  nn = 30;
  mod = simulateModule(rnorm(nS), nn);
  adjacency = abs(cor(mod))^3;
  
  order = c(1:nn);
  
  labels = paste("Gene", c(1:nn));
  circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, radii = c(0.6, 0.6));
  # Fix the overlaping labels - for now this requires semi-manual intervention.
  xOffset = rep(0.01, nn);
  xOffset[16] = -strwidth(labels[16])/3;
  circlePlot(adjacency, labels, order, xLabelOffset = xOffset);
  
  
  # Plot two circles in one plot
  
  circlePlot(adjacency, labels, order, variable.cex.labels = FALSE, center = c(-0.5, -0.5), 
             radii = c(0.35, 0.35));
  
  circlePlot(adjacency, labels, order, startNewPlot = FALSE, 
             variable.cex.labels = FALSE, center = c(0.5, 0.5), 
             radii = c(0.35, 0.35));
  
  
}


library(WGCNA)

# Data 1 - circle corr ----------------------------------------------------

data_1_adjacency <- mat_ord(data_1_cor_r_circular_significant_matrix_noRownames) %>% 
  as.matrix()

diag(data_1_adjacency) <- 0

# pdf('connection_diagram_corr.pdf', width = 10)

par(mar = c(0,0.5,3,0.5));
circlePlot(data_1_adjacency, 
           labels = rownames(data_1_adjacency),
           order = sort(rownames(data_1_adjacency), index.return=TRUE)$ix,
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 2)

par(mar = c(0,0.5,3,0.5));
circlePlot(data_1_adjacency, 
           labels = rownames(data_1_adjacency),
           order = c(1:nrow(data_1_adjacency)),
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 2)


# Data 2 - circle corr ----------------------------------------------------

data_2_adjacency <- mat_ord(data_2_cor_r_circular_significant_matrix_noRownames) %>% 
  as.matrix()

diag(data_2_adjacency) <- 0

par(mar = c(0,0.5,3,0.5));
circlePlot(data_2_adjacency, 
           labels = rownames(data_2_adjacency),
           order = c(1:nrow(data_2_adjacency)),
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 2)

par(mar = c(0,0.5,3,0.5));
circlePlot(data_2_adjacency, 
           labels = rownames(data_2_adjacency),
           order = sort(rownames(data_2_adjacency), index.return=TRUE)$ix,
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 2)

# Data 3 - circle corr ----------------------------------------------------

data_3_adjacency <- mat_ord(data_3_cor_r_circular_significant_matrix_noRownames) %>% 
  as.matrix()

diag(data_3_adjacency) <- 0

par(mar = c(0,0.5,3,0.5));
circlePlot(data_3_adjacency, 
           labels = rownames(data_3_adjacency),
           order = c(1:nrow(data_3_adjacency)),
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 2)

par(mar = c(0,0.5,3,0.5));
circlePlot(data_3_adjacency, 
           labels = rownames(data_3_adjacency),
           order = sort(rownames(data_3_adjacency), index.return=TRUE)$ix,
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 2)



# -------------------------------------------------------------------------

qgraph::qgraph(data_3_adjacency <- mat_ord(data_3_cor_r_circular_significant_matrix_noRownames) %>% 
                 as.matrix()
               , shape='circle', 
               posCol='darkgreen', negCol='darkred', layout='circle', vsize=10)






