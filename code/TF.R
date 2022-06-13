
# Load R libraries ####
library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(ggplotify)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(data.table)
library(Seurat)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)
library(tidyquant)
library(ggrepel)
library(ggsci)
library(scales)
#library(ComplexHeatmap)

# Load Aging Seurat object ####
obj <- readRDS('C:/BGI/×ªÂ¼ÔëÉù/object.rds')

cluster.annotation <- c("¦Â_¢ò", "¦Â_¢ò", "¦Â_¢ò", "¦Â_¢ò", "¦Â_¢ñ",       
                        "IC", "¦Á", "AC", "AC", "IC",          
                        "¦Â_¢ñ", "IC", "IC", "AC", "¦Â_¢ñ",    
                        "¦Ä", "AC", "EC", "EC", "DC", 
                        "IC", "¦Á",  "AC", "DC", "AC", "DC")

names(cluster.annotation) <- levels(obj)
obj <- RenameIdents(obj, cluster.annotation)
# add cell types to meta data
cell.data <- data.table(barcode = colnames(obj),
                        celltype = Idents(obj))
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
obj <- AddMetaData(obj, cell.data, col.name = "celltypes")

## REMOVE unknown ###
#obj <- subset(obj, subset = celltype != 'unknown')

## save annotated data ###
saveRDS(obj, file = 'C:/BGI/×ªÂ¼ÔëÉù/annotated_object.rds')


### add age metadata to seurat object ###
cell.data <- data.table(barcode = colnames(seu.ica),
                        origin = seu.ica$orig.ident)
cell.data$group[cell.data$origin == 'DY210605' | cell.data$origin == 'DY210604'] <- '6Months'
cell.data$group[cell.data$origin == 'XY21060501' | cell.data$origin == 'XY21060502'] <- '10Days'
### delete the origin column ###
cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
cell.data$origin <- NULL
seu.ica <- AddMetaData(seu.ica, cell.data, col.name = "grouping")

### save data ###
saveRDS(seu.ica, file = 'C:/BGI/×ªÂ¼ÔëÉù/annotated_object.rds')


#### load annotated data ###
seu.ica <- readRDS('C:/BGI/×ªÂ¼ÔëÉù/annotated_object.rds')

# Define the celltypes ####
celltypes <- unique(seu.ica@meta.data$celltype)
celltypes <- celltypes[which(!is.na(celltypes))]
celltypes <- setdiff(celltypes, c("red_blood_cells", "low_quality_cells", "Gamma-Delta_T_cells"))

# Define function to calculate euclidean distances accounting for cell number and nUMI ####
getEuclideanDistance <- function(celltype, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(seu.ica, subset = celltypes == celltype)
  expr <- GetAssayData(object = tmp, slot = "counts")
  
  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])
  
  # Down_Sample_Matrix <-function (expr_mat) {
  #   min_lib_size <- min(colSums(expr_mat))
  #   down_sample <- function(x) {
  #     prob <- min_lib_size/sum(x)
  #     return(unlist(lapply(x, function(y) {
  #       rbinom(1, y, prob)
  #     })))
  #   }
  #   down_sampled_mat <- apply(expr_mat, 2, down_sample)
  #   return(down_sampled_mat)
  # }
  # ds_expr <- Down_Sample_Matrix(expr)
  ds_expr <- SampleUMI(data = expr, max.umi =  100000)
  
  nsample <- min(table(tmp@meta.data$grouping)[c("10Days", "6Months")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "6Months")], nsample)
  young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "10Days")], nsample)
  ds_expr_r <- ds_expr[, c(young_r, old_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
    
  calcEuclDist <- function(matr, young, old){
    tmp <- data.matrix(sqrt(matr[genes, young]))
    mean <- rowMeans(sqrt(matr[genes, young]))
    d_young <- distancevector(t(tmp), mean , d="euclid")
    names(d_young) <- young
    
    tmp <- data.matrix(sqrt(matr[genes, old]))
    mean <- rowMeans(sqrt(matr[genes, old]))
    d_old <- distancevector(t(tmp), mean , d="euclid")
    names(d_old) <- old
    
    list(young = d_young, old = d_old)
  }
  ds <- calcEuclDist(matr = ds_expr_r, old = old_r, young = young_r)
  ds
}

# Run for all celltypes ####
res <- lapply(celltypes, function(x) getEuclideanDistance(x, lowcv = T))
names(res) <- celltypes
nulls <- which(unlist(lapply(res, class)) == "NULL")
celltypes[nulls]

# Remove celltypes without enough cells ####
#res <- res[-nulls]
res_original <- res

# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")
sizes <- (-log10(adj_pvals))
sizes[which(sizes < 1)] <- 1
sizes[which(sizes > 4)] <- 4
sizes <- sizes * 0.75
farben <- rep("grey", length(adj_pvals))
farben[which(adj_pvals < 0.05)] <- "purple"

# Generate Fig 2a ####
ord <- rev(order(diffs))

# #### using ggplot2 ###
# tf_df <- do.call(c, res[ord])
# df = as.data.frame(matrix(nrow=0,ncol=3))
# colnames(df) <- c('Cell_types', 'Transcriptional_noise', 'Group')
# for (i in tf_df) {
#   d <- unlist(d)
#   r
#   
# } 
par(mar = c(5,5,2,5))
### COLOR
pal= pal_npg("nrc")(10)
show_col(pal)

pdf("C:/BGI/×ªÂ¼ÔëÉù/fig1a.pdf",width = 7, height = 5)
boxplot(do.call(c, res[ord]), las = 2, outline = F, col = c("#3C5488FF", "#DC0000FF"), 
        ylab = "Transcriptional noise", xaxt = 'n')

axis(1, at = seq(from = 1.5, to =16.5, by = 2), names(res)[ord], las = 2)
legend("topleft", c("young", "old"), col = c("#3C5488FF", "#DC0000FF"), bty = "n", pch = 19)
dev.off()

# Calculate goat means ####
showGoatLevel <- function(celltype){
  tmp <- res_original[[celltype]]
  lapply(tmp, function(x){
    nom <- names(x)
    ids <- unlist(lapply(nom, function(y) strsplit(y, ':', fixed = T)[[1]][1]))
    unlist(lapply(split(x, ids), median))
  })
}
plotGoatLevel <- function(celltype){
  noise <- showGoatLevel(celltype)
  pval <- wilcox.test(noise[[1]], noise[[2]])$p.value
  boxplot(noise, main = paste(celltype, paste('Wilcox P',signif(pval, 2)), sep = '\n'), col = c("#3C5488FF", "#DC0000FF"), names = NA)
}

# Correlate mean differences ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
diffs_goatlevel <- unlist(lapply(names(res_original), function(x){
  tmp <- showGoatLevel(x)
  log2(mean(tmp[[2]]) / mean(tmp[[1]]))
}))
names(diffs_goatlevel) <- names(diffs)

### generate dateframe ###
goatlevel <- data.frame(Goat_level = diffs_goatlevel, Cell_level = diffs)
f <- f.robftest(rlm(diffs ~ diffs_goatlevel), var = "diffs_goatlevel")
### using ggplot for figure ###
p <- ggplot(goatlevel, aes(Goat_level, Cell_level))
p <- p + geom_point(size = 5) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black")) + 
  labs(x = "Goat Level",
       y = "Cell Level",
       title = "Transcriptional noise ratio [old/young log2]") + 
  geom_text_repel(aes(Goat_level, Cell_level, label=rownames(goatlevel)),
                  arrow = arrow(length=unit(0.01, "npc")),
                  force = 1,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=unit(1.6, "lines"), )+
  #geom_abline(coefficients(rlm(diffs ~ diffs_goatlevel)), col = "red")+
  geom_abline(slope  = 1, color = 'red')+
  geom_hline(aes(yintercept=0), linetype="dashed") +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  # annotate("text", label = paste("Robust F-test P = ", signif(f$p.value, 2)), 
  #          x=0.2, y = -0.4, colour = "red")+
  theme_classic(base_size = 16)

#### save the plot ###
ggsave('goat_level.pdf', p, path = 'C:/BGI/×ªÂ¼ÔëÉù/', height = 5, width = 7)


# Generate Fig 2b ####
highlight <- c("EC","¦Â_¢ò","IC","AC","¦Â_¢ñ","¦Á","¦Ä","DC" )

pdf("C:/BGI/×ªÂ¼ÔëÉù/fig1b.pdf",width = 6, height = 5)
plot(diffs_goatlevel, diffs, main = 'Transcriptional noise ratio [old/young log2]', ylab = 'Cell level', xlab = 'Goat level', pch = 19, cex = sizes, col = rgb(0, 0, 0, 0.7))
abline(h = 0, v = 0, lty = 2)
text(diffs_goatlevel[highlight], diffs[highlight], highlight, cex = 0.8)
f <- f.robftest(rlm(diffs ~ diffs_goatlevel), var = "diffs_goatlevel")
abline(coefficients(rlm(diffs ~ diffs_goatlevel)), col = "red")
legend('bottomright', paste("Robust F test P", signif(f$p.value, 2)), bty = 'n')
dev.off()

# Define function to get Spearman correlations ####
getSpearmanCorrelations <- function(celltype){
  print(paste("Working on", celltype))
  
  tmp <- subset(seu.ica, subset = celltypes == celltype)
  expr <- GetAssayData(object = tmp, slot = "counts")
  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])
  # Down_Sample_Matrix <-function (expr_mat) {
  #   min_lib_size <- min(colSums(expr_mat))
  #   down_sample <- function(x) {
  #     prob <- min_lib_size/sum(x)
  #     return(unlist(lapply(x, function(y) {
  #       rbinom(1, y, prob)
  #     })))
  #   }
  #   down_sampled_mat <- apply(expr_mat, 2, down_sample)
  #   return(down_sampled_mat)
  # }
  # ds_expr <- Down_Sample_Matrix(expr)
  ds_expr <- SampleUMI(data = expr, max.umi =  100000)
  nsample <- min(table(tmp@meta.data$grouping)[c("10Days", "6Months")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  }
  old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "6Months")], nsample)
  young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "10Days")], nsample)
  ds_expr_r <- ds_expr[, c(young_r, old_r)]
  cor_young <- cor(method = "spearman", as.matrix(ds_expr_r[, young_r]))
  cor_young <- unlist(cor_young[upper.tri(cor_young)])
  cor_old <- cor(method = "spearman", as.matrix(ds_expr_r[, old_r]))
  cor_old <- unlist(cor_old[upper.tri(cor_old)])
  list(cor_old, cor_young)
}

# Run for all celltypes ####
res <- lapply(celltypes, getSpearmanCorrelations)
names(res) <- celltypes
nulls <- which(unlist(lapply(res, class)) == "NULL")
celltypes[nulls]

# Remove celltypes without enough cells ####
res_cor <- res

# Generate Fig 2c ####
diffs_cor <- unlist(lapply(res_cor, function(x) log2(mean(1 - x[[1]]) / mean(1 - x[[2]]))))

### make dataframe ###
correlation <- data.frame(Spearman_correlation = diffs_cor,Euclidean_distance = diffs )

### using ggplot for figures ###
p <- ggplot(correlation, aes(Spearman_correlation, Euclidean_distance))
f <- f.robftest(rlm(diffs ~ diffs_cor), var = "diffs_cor")
cor <- coefficients(rlm(diffs ~ diffs_cor))
p <- p + geom_point(size = 5) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black")) + 
  labs(x = "1 - Spearman correlation [old/young log2]",
       y = "Euclidean distance [old/young log2]",
       title = "Transcriptional noise") + 
  geom_text_repel(aes(Spearman_correlation, Euclidean_distance, 
                      label=rownames(correlation)),
                  arrow = arrow(length=unit(0.01, "npc")),
                  force = 1,
                  box.padding=unit(0.5, "lines"), 
                  point.padding=unit(1.6, "lines"), )+
  #geom_abline(coefficients(rlm(diffs ~ diffs_goatlevel)), col = "red")+
  geom_abline(slope  = cor[2], intercept = cor[1], color = 'red')+
  geom_hline(aes(yintercept=0), linetype="dashed") +
  geom_vline(aes(xintercept=0), linetype="dashed")+
  annotate("text", label = paste("Robust F-test P = ", signif(f$p.value, 2)), 
           x=-0.05, y = -0.4, colour = "red")+
  theme_classic(base_size = 16)

#### save the plot ###
ggsave('correlation.pdf', p, path = 'C:/BGI/×ªÂ¼ÔëÉù/', height = 5, width = 7)


pdf("C:/BGI/×ªÂ¼ÔëÉù/fig1c.pdf",width = 6, height = 5)
plot(diffs_cor, diffs, ylab = "Euclidean distance [old/young log2]", xlab = "1 - Spearman correlation [old/young log2]", main = "Transcriptional noise", pch = 19, cex = sizes, col = rgb(0, 0, 0, 0.7))
abline(v = 0, h = 0, lty = 2)
text(diffs_cor[highlight], diffs[highlight], highlight) #col = c("#DF536B","#61D04F","#2297E6"))
f <- f.robftest(rlm(diffs ~ diffs_cor), var = "diffs_cor")
abline(coefficients(rlm(diffs ~ diffs_cor)), col = "red")
legend("topleft", paste("Robust F test P", signif(f$p.value, 2)), bty = "n")
dev.off()

# Generate Fig 2d ###
plot(density(1- res_cor[["AC"]][[1]]), col = "red", ylab = "Density", xlab = "1 - Spearman Correlation coefficient", main = "AC", lwd = 2)
lines(density(1- res_cor[["AC"]][[2]]), col = "blue", lwd = 2)
correl <- ks.test(res_cor[["AC"]][[2]], res_cor[["AC"]][[1]])
legend("topright", paste("KS P", signif(correl$p.value, 2)), bty = "n")
legend("topleft", c("young", "old"), col = c("blue", "red"), bty = "n", pch = 19)

### using ggplot2 for Fig 2d ###
### form dataframe ###
AC_df_old <- as.data.frame(1- res_cor[["AC"]][[1]])
#AC_df_old <- as.data.frame(1- res_cor$AC[1])
AC_df_old$Group <- factor('Old')
colnames(AC_df_old)[1]<-'Correlation'
mean_old <- mean(AC_df_old$Correlation)

AC_df_young <- as.data.frame(1- res_cor[["AC"]][[2]])
AC_df_young$Group <- factor('Young')
colnames(AC_df_young)[1]<-'Correlation'
AC_df <- rbind(AC_df_old, AC_df_young)
mean_young <- mean(AC_df_young$Correlation)

## using orginal res ###
AC_df_old <- as.data.frame(1- res_original[["AC"]][[1]])
AC_df_old$Group <- factor('Old')
colnames(AC_df_old)[1]<-'Correlation'
mean_old <- mean(AC_df_old$Correlation)

AC_df_young <- as.data.frame(1- res_original[["AC"]][[2]])
AC_df_young$Group <- factor('Young')
colnames(AC_df_young)[1]<-'Correlation'
AC_df <- rbind(AC_df_old, AC_df_young)
mean_young <- mean(AC_df_young$Correlation)
### geom_density
p <- ggplot(data = AC_df, mapping = aes(x = Correlation,color = Group))
correl <- ks.test(res_cor[["AC"]][[2]], res_cor[["AC"]][[1]])
p <- p + geom_line(stat = "density")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black")) + 
  labs(x = "Transcriptional noise",
       y = "Density",
       title = "AC") +
  annotate("text", label = paste("KS test P =", signif(correl$p.value, 2)),
           x=0.9, y = 6, colour = "red")+
  #geom_abline(coefficients(rlm(diffs ~ diffs_goatlevel)), col = "red")+
  theme_classic()

#### save the plot ###
ggsave('density_orig.png', p, path = 'C:/BGI/×ªÂ¼ÔëÉù/', height = 5, width = 8)

##### obtain gene expression level for every cell ####
genes_of_interest <- c('PNLIP', 'PNLIPRP2', 'PLA2G1B', 'CTRC', 'CPB1', 'CELA1')
genes_of_interest <- c('INS', 'SST', 'PPY', 'CTRC', 'CEL',
                       'CPB1','SOX9','KRT19')
### define a function to obtain expression level of goi for each celltype ###
getExpressionLevel <- function(celltype, barcode, goi){
  print(paste("Working on", celltype))
  library(hopach)
  tmp <- subset(seu.ica, cells = barcode, subset = celltypes == celltype)
  expr <- tmp@assays$RNA@data[goi,]
  expr
}
exp_level <- c()
celltype <- celltypes
### run for all the celltypes and goi ###
for (i in genes_of_interest) {
  i <- lapply(celltype, function(x) getExpressionLevel(x,i))
  #exp_level <- append(exp_level, exp_level_new)
  #names(exp_level_new) <- celltypes 
}

### AC ###
### Obtain the transcriptional noise for AC ###
AC_TN <- res_original$AC ### 1396 cells' transcriptional noise in total
barcode_young <- names(AC_TN$young)
barcode_old <- names(AC_TN$old)
barcode_AC <- append(barcode_young, barcode_old)

### obtain the corresponding AC cells using these barcodes ###
AC <- subset(seu.ica, cells = barcode_AC)

### obtain the gene expression of goi for AC ###
goi_level_AC <- lapply(genes_of_interest, function(x) getExpressionLevel('AC',barcode_AC, x))
names(goi_level_AC) <- genes_of_interest
### select equal number of cells to transcriptional noise by sampling 

### form a dataframe-------PNLIP
ac <- data.frame(unlist(AC_TN), goi_level_AC$PNLIP)
### rename the column name of the dataframe ###
colnames(ac) <- c('Transcriptional noise', 'Gene expression level')
p1 <- ggplot(data = ac, aes(x = `Transcriptional noise`, y = `Gene expression level`))+
      geom_point(alpha=0.5,size = 0.1)+ 
      geom_smooth(method = 'loess')+ theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(),axis.line =element_line(colour = "black"))
ggsave('AC_PNLIP.png', p1, path = 'C:/BGI/×ªÂ¼ÔëÉù/', width = 10, height = 6)
### form a dataframe-------PNLIPRP2
ac <- data.frame(unlist(AC_TN), goi_level_AC$PNLIPRP2)

colnames(ac) <- c('Transcriptional noise', 'Gene expression level')
p2 <- p1 + geom_point(data = ac, aes(x = `Transcriptional noise`, y = `Gene expression level`))


#### form a general dataframe for AC for all genes ####
ac_all <- data.frame(unlist(AC_TN), goi_level_AC)
colnames(ac_all)[1] <- 'Transcriptional Noise'
ggplot(data = ac_all, aes(x = `Transcriptional Noise`, y = `Gene expression level`))+
  geom_point(alpha=0.5,size = 0.1, )+ 
  geom_smooth(method = 'loess')+ theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line =element_line(colour = "black"))


### adding a running mean to the plot
## obtain the average expression of goi in these cells ##
ave_goi_ac <- AverageExpression(AC, features = genes_of_interest)

ac_all_i<-ac_all[order(ac_all[,1],decreasing = F),] 
### reverse the dataframe ###
ac_all_i <- t(ac_all_i)
### change the colnames to the number of transcriptional noise ###
colnames(ac_all_i)<- ac_all_i[1,]
### remove the first row ###
ac_all_i <- as.data.frame(ac_all_i)
ac_all_i <- ac_all_i[-1,]
###heat map ###
g = as.ggplot(pheatmap(ac_all_i, cluster_cols = FALSE, show_colnames = F))
ggsave('heatmap_AC.png', g, path = 'C:/BGI/×ªÂ¼ÔëÉù/',width = 10,height = 5)


### generate new dataframe ###
ac_all_n <- ac_all
## PNLIP ##
colnames(ac_all_n)[2] <- 'Gene Expression'
ac_all_n$Genes <- as.factor('PNLIP')
ac_all_n <- ac_all_n[,-c(3,4,5,6,7)]
PNLIP <- ac_all_n

## PNLIPRP2 ##
PNLIPRP2 <- ac_all_n[,c(1,3)]
colnames(PNLIPRP2)[2] <- 'Gene Expression'
PNLIPRP2$Genes <- as.factor('PNLIPRP2')

## PLA2G1B ##
PLA2G1B <- ac_all_n[,c(1,4)]
colnames(PLA2G1B)[2] <- 'Gene Expression'
PLA2G1B$Genes <- as.factor('PLA2G1B')

## CTRC ##
CTRC <- ac_all_n[,c(1,5)]
colnames(CTRC)[2] <- 'Gene Expression'
CTRC$Genes <- as.factor('CTRC')

## CPB1 ##
CPB1 <- ac_all_n[,c(1,6)]
colnames(CPB1)[2] <- 'Gene Expression'
CPB1$Genes <- as.factor('CPB1')

## CELA1 ##
CELA1 <- ac_all_n[,c(1,7)]
colnames(CELA1)[2] <- 'Gene Expression'
CELA1$Genes <- as.factor('CELA1')

## INS ##
INS <- ac_all_n[,c(1,2)]
colnames(INS)[2] <- 'Gene Expression'
INS$Genes <- as.factor('INS')

## SST ##
SST <- ac_all_n[,c(1,3)]
colnames(SST)[2] <- 'Gene Expression'
SST$Genes <- as.factor('SST')

# ## GCG ##
# GCG <- ac_all_n[,c(1,4)]
# colnames(GCG)[2] <- 'Gene Expression'
# GCG$Genes <- as.factor('GCG')

## PPY ##
PPY <- ac_all_n[,c(1,5)]
colnames(PPY)[2] <- 'Gene Expression'
PPY$Genes <- as.factor('PPY')

## CTRC ##
CTRC <- ac_all_n[,c(1,6)]
colnames(CTRC)[2] <- 'Gene Expression'
CTRC$Genes <- as.factor('CTRC')

## CEL ##
CEL <- ac_all_n[,c(1,7)]
colnames(CEL)[2] <- 'Gene Expression'
CEL$Genes <- as.factor('CEL')

## CPB1 ##
CPB1 <- ac_all_n[,c(1,8)]
colnames(CPB1)[2] <- 'Gene Expression'
CPB1$Genes <- as.factor('CPB1')

## SOX9 ##
SOX9 <- ac_all_n[,c(1,9)]
colnames(SOX9)[2] <- 'Gene Expression'
SOX9$Genes <- as.factor('SOX9')

# ## KRT8 ##
# KRT8 <- ac_all_n[,c(1,10)]
# colnames(KRT8)[2] <- 'Gene Expression'
# KRT8$Genes <- as.factor('KRT8')

## KRT19 ##
KRT19 <- ac_all_n[,c(1,11)]
colnames(KRT19)[2] <- 'Gene Expression'
KRT19$Genes <- as.factor('KRT19')

### merge all the dataframes together ###
l <- data.frame()
for(i in genes_of_interest){
  df.now <- get(i)
  l <- rbind(l, df.now)
}

colnames(l) <- c('Transcriptional_noise', 'Gene_Expression', "Genes")

### method 2 for merging the dataframes ###
df_list <- list(PNLIP,PNLIPRP2,PLA2G1B,CTRC,CPB1,CELA1)

df_total <- Reduce(function(x,y) merge(x,y,all=T),df_list)

color <- c('#fc9272','#fb6a4a','#de2d26',### red
           '#6baed6','#4292c6','#084594',### blue
           '#756bb1','#54278f') ###purple
### use ggplot to put all the genes vs TN on a single plot ###
p <- ggplot(data=l, mapping = aes(
            x = Transcriptional_noise,
            y = Gene_Expression,
            color = Genes))

p <- p + geom_point(alpha=0.5,size = 0.1) +
  geom_smooth(method="loess") + theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line =element_line(colour = "black")) + 
  labs(x = "Transcriptional Noise",
       y = "Gene Expression",
       title = "AC")
p <- p + scale_color_manual(values = color)

### save the plot ###
ggsave('scatterplot_AC_allgenes.png', p, path = 'C:/BGI/×ªÂ¼ÔëÉù/',width = 8,height = 5)

### save the NEW plot ###
ggsave('scatterplot_AC_all_marker_genes.pdf', p, path = 'C:/BGI/×ªÂ¼ÔëÉù/',width = 8,height = 5)

