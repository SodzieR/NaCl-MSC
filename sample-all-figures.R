
# wd ----------------------------------------------------------------------
setwd('D:\\!ZMRiI\\!MSC_PAPER')

# Dependencies ------------------------------------------------------------

library(clusterProfiler)
library(AnnotationHub)
library(readxl)
library(ggplot2)
library(RColorBrewer)


# Dependencies ------------------------------------------------------------
library(tidyverse)
library(readxl)
# rstatix is for kendall's corr matrix
library(rstatix)
# library(Hmisc) 
library(corrplot)
library(circlize)
setwd('D:/!ZMRiI/!MSC_PAPER/')

#####################################################
#####################################################

#raw_read_counts <- read_excel('raw_read_counts.xlsx')

raw_read_counts <- read.delim('all_data_normalized_TPM.txt', header = TRUE)
######################################################

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
pdf('all-sample.pdf', width = 16, height = 10)
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
    geom_bar(stat = 'identity', position = 'dodge', width = 0.25) +
    scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
    geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
    coord_flip() +
    ggtitle('GO_GSEA 2vs3') +
    labs(y= '-log10(padj)', x = "Description") +
    theme_classic())
(GO_GSEA_2vs3_df_plot_NES <- ggplot(preanalyzded_2vs3, aes(factor(factor(Description, levels = Description), levels = Description), `-log10(padj)`, fill = NES)) +
    geom_bar(stat = 'identity', position = 'dodge', width = 0.25) +
    scale_fill_gradient2(low = 'blue', mid = 'white', midpoint = 0, high = 'red') +
    geom_text(aes(label = round(enrichmentScore, 3)), hjust = -0.1) +
    coord_flip() +
    ggtitle('GO_GSEA 2vs3') +
    labs(y= '-log10(padj)', x = "Description") +
    theme_classic())

corrplot(as.matrix(data_1_cor_r_corrplot), type="lower", order="original", tl.col = 'black',
         p.mat = as.matrix(data_1_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'MICE1,2,3,4,5', method = 'square', tl.srt = 45,mar=c(0,0,1,0))

corrplot(as.matrix(data_2_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_2_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'MICIE6,7,8,9,10', method = 'square', tl.srt = 45, mar=c(0,0,1,0))

corrplot(as.matrix(data_3_cor_r_corrplot), type="lower", order="alphabet", tl.col = 'black',
         p.mat = as.matrix(data_3_cor_p_corrplot), sig.level = 0.05, insig = "blank", 
         title = 'MICE11,12,13,14,15', method = 'square', tl.srt = 45, mar=c(0,0,1,0))

##########################################
#
data_1_cor_r_circular_significant <- data_1_cor_r %>% 
  cor_gather() %>% 
  dplyr::filter(p < 0.05)
rownames(data_1_cor_r_circular_significant) <- NULL

data_1_cor_r_circular_significant_matrix <- cor_spread(data_1_cor_r_circular_significant)

data_1_cor_r_circular_significant_matrix_noRownames <- data_1_cor_r_circular_significant_matrix %>%
  column_to_rownames(var = 'rowname')

data_1_cor_r_circular_significant_matrix_noRownames[is.na(data_1_cor_r_circular_significant_matrix_noRownames)] <- 0


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

##################################################################################################

#
data_2_cor_r_circular_significant <- data_2_cor_r %>% 
  cor_gather() %>% 
  dplyr::filter(p < 0.05)
rownames(data_2_cor_r_circular_significant) <- NULL

data_2_cor_r_circular_significant_matrix <- cor_spread(data_2_cor_r_circular_significant)

data_2_cor_r_circular_significant_matrix_noRownames <- data_2_cor_r_circular_significant_matrix %>%
  column_to_rownames(var = 'rowname')

data_2_cor_r_circular_significant_matrix_noRownames[is.na(data_2_cor_r_circular_significant_matrix_noRownames)] <- 0


mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]
data_2_cor_r_circular_significant_matrix_noRownames_sorted <- mat_ord(data_2_cor_r_circular_significant_matrix_noRownames)

data_2_cor_r_circular_significant_matrix_noRownames_sorted[seq(from = 1, to = nrow(data_2_cor_r_circular_significant_matrix_noRownames_sorted), by = 2),] <- 0

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
link_colors = colorRamp2(range(data_2_cor_r_circular_significant_matrix_noRownames_sorted2), c("red", "blue"), transparency = 0.25)

#pdf("SAMPLE-GROUP2-chord-significant.pdf", width=7, height=7)

col_fun_cor_3 = colorRamp2(c(-1, 0, 1), c("red", "white", "blue"))
lgd_links = ComplexHeatmap::Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_cor_3, 
                                   title_position = "topleft", title = "Link")

chordDiagram(as.matrix(data_2_cor_r_circular_significant_matrix_noRownames_sorted2), grid.col = colfunc,
             annotationTrack = c("name", "grid"), col = link_colors)
title("Group 2 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

#############################################################################

data_3_cor_r_circular_significant <- data_3_cor_r %>% 
  cor_gather() %>% 
  dplyr::filter(p < 0.05)
rownames(data_3_cor_r_circular_significant) <- NULL

data_3_cor_r_circular_significant_matrix <- cor_spread(data_3_cor_r_circular_significant)

data_3_cor_r_circular_significant_matrix_noRownames <- data_3_cor_r_circular_significant_matrix %>%
  column_to_rownames(var = 'rowname')

data_3_cor_r_circular_significant_matrix_noRownames[is.na(data_3_cor_r_circular_significant_matrix_noRownames)] <- 0


mat_ord <- function(mx) mx[, c(rownames(mx), setdiff(colnames(mx), rownames(mx)))]
data_3_cor_r_circular_significant_matrix_noRownames_sorted <- mat_ord(data_3_cor_r_circular_significant_matrix_noRownames)

data_3_cor_r_circular_significant_matrix_noRownames_sorted[seq(from = 1, to = nrow(data_3_cor_r_circular_significant_matrix_noRownames_sorted), by = 2),] <- 0


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

data_1_adjacency <- mat_ord(data_1_cor_r_circular_significant_matrix_noRownames) %>% 
  as.matrix()

diag(data_1_adjacency) <- 0

# pdf('connection_diagram_corr.pdf', width = 10)
col_fun_cor_4 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd_links3 = ComplexHeatmap::Legend(at = c(-1, -0.5, 0, 0.5, 1), col_fun = col_fun_cor_4, 
                                   title_position = "topleft", title = "Links")

par(mar = c(0,0.5,3,0.5));
circlePlot(data_1_adjacency, 
           labels = rownames(data_1_adjacency),
           order = sort(rownames(data_1_adjacency), index.return=TRUE)$ix,
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 1.5)
title("Group 1 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links3, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

data_2_adjacency <- mat_ord(data_2_cor_r_circular_significant_matrix_noRownames) %>% 
  as.matrix()

diag(data_2_adjacency) <- 0

par(mar = c(0,0.5,3,0.5));
circlePlot(data_2_adjacency, 
           labels = rownames(data_2_adjacency),
           order = sort(rownames(data_2_adjacency), index.return=TRUE)$ix,
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 1.5)
title("Group 2 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links3, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))

data_3_adjacency <- mat_ord(data_3_cor_r_circular_significant_matrix_noRownames) %>% 
  as.matrix()

diag(data_3_adjacency) <- 0
par(mar = c(0,0.5,3,0.5));
circlePlot(data_3_adjacency, 
           labels = rownames(data_3_adjacency),
           order = sort(rownames(data_3_adjacency), index.return=TRUE)$ix,
           radii = c(0.9, 0.7), lineColors = greenWhiteRed(50),min.cex.labels = 1.1, max.cex.labels = 1.1,
           xLabelOffset = 0, yLabelOffset = 0,plotBox = c(-1.7, 1.7, -1, 1), variable.cex.points = F,
           min.cex.points = 1.5)
title("Group 3 - all significant correlations", cex = 0.8)
ComplexHeatmap::draw(lgd_links3, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))


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
# Heatmap(TPM_heatmap_in_log2_scaled, name = 'MICE 1:15', col = col_fun1,
#        cluster_rows = F, cluster_columns = FALSE, show_row_names = F,
#       border = T, show_column_names = TRUE, column_title = 'Z-scored log2(TPM)',
#       top_annotation = HeatmapAnnotation(
#          empty = anno_empty(border = FALSE),
#          foo = anno_block(gp = gpar(fill = 'black'), labels = c('Group 1', 'Group 2', 'Group 3'),
#                           labels_gp = gpar(col = "white", fontsize = 10))),
#        column_split = rep(1:3, each = 5)
#        )


#Clustering problem - not enough ram

col_fun2 = circlize::colorRamp2(c(min(TPM_heatmap_in_log2_scaled), 0, max(TPM_heatmap_in_log2_scaled)),
                                c('blue','white','red'))

ann <- data.frame(metadata$group)
colnames(ann) <- c('group')
colours <- list('group' = c('1' = 'red', 
                            '2' = 'green',
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

#pdf("bwr_Defauult-clustering_z-scored_log2TPM_7x7.pdf", width=7, height=7)
draw(temp2, heatmap_legend_side="left", annotation_legend_side="right")
#dev.off()


dev.off()
# Setup --------------------