rm(list=ls())

library(gplots)
library(ggplot2)
library(VennDiagram)
library(viridis)
library(readxl)
library(conflicted)
library(dplyr)
library(tidyverse) # contains ggplot for volcano plots
library(RColorBrewer) # to color plots
library(ggrepel) # annotating volcano plot
library(ggfortify) # PCA plots
library(factoextra) # PCA plot ellipses
library(MSnSet.utils) # Different PCA plot option (simpler?)
library('corrr') # PCA plot 3rd option
library(ggcorrplot) # PCA plot 3rd option
library("FactoMineR") # PCA plot 3rd option
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library("xlsx")
library(pheatmap)
library('HybridMTest') # for ANOVA stats
library('BiocManager')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("HybridMTest")

conflicts_prefer(base::setdiff, base::intersect, dplyr::filter())


setwd("C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx")

# Loading full dataset

Chamberlain_proteomics_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/dataset_full.xlsx"))

Chamberlain_proteomics_data_transposed <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/dataset_full_transposed.xlsx"))

# Loading 2-group comparison DEP datasets ----

#DEPs calculated and filtered in Excel 

mdx_vs_WT_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsWTDEPs.xlsx"))

uDysmdx_vs_WT_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/uDysmdxvsWTDEPs.xlsx"))

MidiDysmdx_vs_WT_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/midiDysmdxvsWTDEPs.xlsx"))

flDysmdx_vs_WT_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/flDysmdxvsWTDEPs.xlsx"))

mdx_vs_uDysmdx_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsuDysmdxDEPs.xlsx"))

mdx_vs_MidiDysmdx_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsmidiDysmdxDEPs.xlsx"))

mdx_vs_flDysmdx_DEPs <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsflDysmdxDEPs.xlsx"))

----
  
# Loading 2-group volcano plot datasets ----
  
mdx_vs_WT_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsWT_volcano_data.xlsx"))

uDysmdx_vs_WT_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/uDysmdxvsWT_volcano_data.xlsx"))

MidiDysmdx_vs_WT_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/midiDysmdxvsWT_volcano_data.xlsx"))

flDysmdx_vs_WT_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/flDysmdxvsWT_volcano_data.xlsx"))

mdx_vs_uDysmdx_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsuDysmdx_volcano_data.xlsx"))

mdx_vs_MidiDysmdx_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsmidiDysmdx_volcano_data.xlsx"))

mdx_vs_flDysmdx_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsflDysmdx_volcano_data.xlsx"))

uDysmdx_vs_flDysmdx_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/uDysmdxvsflDysmdx_volcano_data.xlsx"))

midiDysmdx_vs_flDysmdx_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/midiDysmdxvsflDysmdx_volcano_data.xlsx"))

uDysmdx_vs_MidiDysmdx_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/uDysmdxvsMidiDysmdx_volcano_data.xlsx"))

----
  
  
##########################
# PCA Plot
##########################

PCA_no_legend <- Chamberlain_proteomics_data_transposed[2:2625]

PCA_no_missing <- PCA_no_legend[ , colSums(is.na(PCA_no_legend))==0]

pca_res_chamberlain <- prcomp(PCA_no_missing, scale. = FALSE)

PCi_chamberlain <- data.frame(pca_res_chamberlain$x,Legend=Chamberlain_proteomics_data_transposed$Group)

ggplot(PCi_chamberlain,aes(x=PC1,y=PC2,col=Legend,frame=T))+
  geom_point(size=3,alpha=0.6) + #Size and alpha just for fun
  theme_minimal() +
  theme(text = element_text(size=15)) +
  stat_ellipse(alpha=0.03,geom = "polygon", aes(fill=Legend)) +
  scale_color_manual(values = c("#46B1C9","#D66853","#FCD351","#A4BAB7","#3C0919"))

# Biplot and other PCA plot option

corr_matrix_chamberlain <- cor(PCA_no_missing)

data_normalized_chamberlain <- scale(PCA_no_missing)
head(data_normalized_chamberlain)

data.pca_chamberlain <- prcomp(data_normalized_chamberlain)
summary(data.pca_chamberlain)

# Scree plot
fviz_eig(data.pca_chamberlain, addlabels = TRUE)

# Biplot
fviz_pca_var(data.pca_chamberlain, col.var = "cos2",
             gradient.cols=c("black","turquoise","purple"),
             ggtheme = theme_minimal(),
             select.var = list(name = NULL, cos2 = NULL, contrib=10),
             repel=TRUE)

cos2_chamberlain <- fviz_cos2(data.pca_chamberlain, choice = "var", axes = 1:2, 
                            sort.val = "desc",
                            top=20) # Can change top to any numerical value
cos2_chamberlain_PC1only <- fviz_cos2(data.pca_chamberlain, choice = "var", axes = 1, 
                                    sort.val = "desc",
                                    top=20) # Can change top to any numerical value

topfeatures_chamberlain <- as.data.frame(cos2_chamberlain$data)
topfeatures_chamberlain_PC1 <- as.data.frame(cos2_chamberlain_PC1only$data)

#write.xlsx(topfeatures_chamberlain, 
#           "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/PC1_PC2_top_features.xlsx")
#write.xlsx(topfeatures_chamberlain_PC1, 
#          "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/PC1_top_features.xlsx")

# Most useful and visually appealing biplot
# Can change contrib number to add or subtract vectors
fviz_pca_biplot(data.pca_chamberlain, axes = c(1,2), label = "var", repel = TRUE,
                col.var = "black",
                geom = "point",
                pointsize = 2,
                habillage=Chamberlain_proteomics_data_transposed$Group,
                palette = c("#212D40","#D66853","#364156","#7D4E57","#440154"), 
                addEllipses=TRUE, ellipse.level=0.95,
                ellipse.type = "t",
                select.var = list(contrib = 1),
                ggtheme = theme_minimal())

-----

  
##########################
# Volcano Plots
##########################

# Defining Volcano plot function 
volcano_fun <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`L2FC` >= 1 & `log10pval` >= 1.3, "A",
                     if_else(`L2FC` <= -1 & `log10pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`L2FC`, y=`log10pval`, colour = threshold)) +
    geom_point(alpha = 1, size = 2.5, shape = 16) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 6), xlim = c(-4, 5)) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#7D4E57", "B"= "#5E7297", 
                                   "C" = "#A2A7A5"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    #xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    theme_minimal() + #Set the theme
    theme(legend.position = c(.98, .98),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",linewidth =1)) +
    theme(text = element_text(size=14)) 
}

# Defining Volcano plot function 2
volcano_fun_2 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`L2FC` >= 1 & `log10pval` >= 1.3, "A",
                     if_else(`L2FC` <= -1 & `log10pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`L2FC`, y=`log10pval`, colour = threshold)) +
    geom_point(alpha = 1, size = 2.5, shape = 16) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = c(-1,1), linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 6), xlim = c(-4, 5)) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#7D4E57", "B"= "#5E7297", 
                                   "C" = "#A2A7A5"),
                        labels = c("Not Significant", "Upregulated", "Downregulated")) +
    #xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    theme_minimal() + #Set the theme
    theme(legend.position = c(.98, .98),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",linewidth =1)) +
    theme(text = element_text(size=14)) 
}

# mdx vs WT volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdx_vs_WT_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
mdx_vs_WT_volcano_data$diffexpressed[mdx_vs_WT_volcano_data$L2FC > 1 & 
                                       mdx_vs_WT_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
mdx_vs_WT_volcano_data$diffexpressed[mdx_vs_WT_volcano_data$L2FC < -1 & 
                                       mdx_vs_WT_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(mdx_vs_WT_volcano_data)

# uDysmdx vs WT volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
uDysmdx_vs_WT_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
uDysmdx_vs_WT_volcano_data$diffexpressed[uDysmdx_vs_WT_volcano_data$L2FC > 1 & 
                                       uDysmdx_vs_WT_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
uDysmdx_vs_WT_volcano_data$diffexpressed[uDysmdx_vs_WT_volcano_data$L2FC < -1 & 
                                       uDysmdx_vs_WT_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(uDysmdx_vs_WT_volcano_data)

# midiDysmdx vs WT volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
MidiDysmdx_vs_WT_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
MidiDysmdx_vs_WT_volcano_data$diffexpressed[MidiDysmdx_vs_WT_volcano_data$L2FC > 1 & 
                                           MidiDysmdx_vs_WT_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
MidiDysmdx_vs_WT_volcano_data$diffexpressed[MidiDysmdx_vs_WT_volcano_data$L2FC < -1 & 
                                           MidiDysmdx_vs_WT_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(MidiDysmdx_vs_WT_volcano_data)

# flDysmdx vs WT volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
flDysmdx_vs_WT_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
flDysmdx_vs_WT_volcano_data$diffexpressed[flDysmdx_vs_WT_volcano_data$L2FC > 1 & 
                                              flDysmdx_vs_WT_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
flDysmdx_vs_WT_volcano_data$diffexpressed[flDysmdx_vs_WT_volcano_data$L2FC < -1 & 
                                              flDysmdx_vs_WT_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(flDysmdx_vs_WT_volcano_data)

# mdx vs uDysmdx volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdx_vs_uDysmdx_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
mdx_vs_uDysmdx_volcano_data$diffexpressed[mdx_vs_uDysmdx_volcano_data$L2FC > 1 & 
                                            mdx_vs_uDysmdx_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
mdx_vs_uDysmdx_volcano_data$diffexpressed[mdx_vs_uDysmdx_volcano_data$L2FC < -1 & 
                                            mdx_vs_uDysmdx_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(mdx_vs_uDysmdx_volcano_data)

# mdx vs MidiDysmdx volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdx_vs_MidiDysmdx_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
mdx_vs_MidiDysmdx_volcano_data$diffexpressed[mdx_vs_MidiDysmdx_volcano_data$L2FC > 1 & 
                                            mdx_vs_MidiDysmdx_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
mdx_vs_MidiDysmdx_volcano_data$diffexpressed[mdx_vs_MidiDysmdx_volcano_data$L2FC < -1 & 
                                            mdx_vs_MidiDysmdx_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(mdx_vs_MidiDysmdx_volcano_data)

# mdx vs flDysmdx volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdx_vs_flDysmdx_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
mdx_vs_flDysmdx_volcano_data$diffexpressed[mdx_vs_flDysmdx_volcano_data$L2FC > 1 & 
                                            mdx_vs_flDysmdx_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
mdx_vs_flDysmdx_volcano_data$diffexpressed[mdx_vs_flDysmdx_volcano_data$L2FC < -1 & 
                                            mdx_vs_flDysmdx_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun(mdx_vs_flDysmdx_volcano_data)

# uDysmdx vs flDysmdx volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
uDysmdx_vs_flDysmdx_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
uDysmdx_vs_flDysmdx_volcano_data$diffexpressed[uDysmdx_vs_flDysmdx_volcano_data$L2FC > 1 & 
                                                uDysmdx_vs_flDysmdx_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
uDysmdx_vs_flDysmdx_volcano_data$diffexpressed[uDysmdx_vs_flDysmdx_volcano_data$L2FC < -1 & 
                                                uDysmdx_vs_flDysmdx_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun_2(uDysmdx_vs_flDysmdx_volcano_data)

# flDysmdx vs MidiDysmdx volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
midiDysmdx_vs_flDysmdx_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
midiDysmdx_vs_flDysmdx_volcano_data$diffexpressed[midiDysmdx_vs_flDysmdx_volcano_data$L2FC > 1 & 
                                                 midiDysmdx_vs_flDysmdx_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
midiDysmdx_vs_flDysmdx_volcano_data$diffexpressed[midiDysmdx_vs_flDysmdx_volcano_data$L2FC < -1 & 
                                                 midiDysmdx_vs_flDysmdx_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun_2(midiDysmdx_vs_flDysmdx_volcano_data)

# uDysmdx vs MidiDysmdx volcano plot

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
uDysmdx_vs_MidiDysmdx_volcano_data$diffexpressed <- NA
# if log2Foldchange > 1 and -log10pval > 1.3, set as "UP"
uDysmdx_vs_MidiDysmdx_volcano_data$diffexpressed[uDysmdx_vs_MidiDysmdx_volcano_data$L2FC > 1 & 
                                                    uDysmdx_vs_MidiDysmdx_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < -1 and -log10pval > 1.3, set as "DOWN"
uDysmdx_vs_MidiDysmdx_volcano_data$diffexpressed[uDysmdx_vs_MidiDysmdx_volcano_data$L2FC < -1 & 
                                                    uDysmdx_vs_MidiDysmdx_volcano_data$log10pval > 1.3] <- "DOWN"

volcano_fun_2(uDysmdx_vs_MidiDysmdx_volcano_data)


##########################
# Venn Diagrams
##########################

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
# Suppresses data file generation each time a Venn diagram is generated
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Venn diagram all groups

All_groups_4D_venn <- display_venn(All_groups_4D_venn_list <- list(
  'mdx4cv vs WT' = mdx_vs_WT_ANOVA_filtered$Gene.Symbol,
  'μDys-mdx4cv vs WT' = uDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol,
  'MidiDys-mdx4cv vs WT' = MidiDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol,
  'flDys-mdx4cv vs WT' = flDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx4cv vs WT DEPs",
    "μDys-mdx4cv vs WT DEPs",
    "MidiDys-mdx4cv vs WT DEPs", 
    "flDys-mdx4cv vs WT DEPs"),
  fill = c("#C44917","#59594A","#F2E86D","#51A3A3"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-12, -3, -20, 20),
  cat.dist = c(0.12, .13, 0.088, 0.084),
  main = "DEP comparison between dystrophin transgenic constructs",
  main.fontface = "bold"
) 
------

##########################################
# Row-wise one-way ANOVA for full dataset
##########################################

#fit_aov <- function(col) {
#  aov(col ~ Group, data = Chamberlain_proteomics_data_transposed)
#}
  
#anovas <- map(Chamberlain_proteomics_data_transposed[, 2:ncol(Chamberlain_proteomics_data_transposed)], fit_aov)

#summary(anovas$Zyx)

#TukeyHSD(anovas)

# Define row.oneway.anova function from https://rdrr.io/bioc/HybridMTest/src/R/row.oneway.anova.R
# Also pre-set function in HybridMTest package that is part of BioConductoR

row.oneway.anova <-
  function(Y,grplbl)
  {
    ugrps<-unique(grplbl)
    ngrps<-length(ugrps) #number of groups
    ngenes<-dim(Y)[1]
    GrandM<-rowMeans(Y)   # overall mean
    
    SST<-rowSums((Y-GrandM)^2) # total sum of squares for each gene
    
    grp.mean<-matrix(NA,ngenes,ngrps)  # group mean matrix, rows for genes, each column for a different group
    grp.SSW<-matrix(NA,ngenes,ngrps)  # within-group sums of squares for each gene
    n<-rep(NA,ngrps)  # vector with group-specific sample sizes
    for (i in 1:ngrps)
    {
      grp.mtch<-(grplbl==ugrps[i])
      n[i]<-sum(grp.mtch)
      grp.mean[,i]<-rowMeans(Y[,grp.mtch])
      grp.SSW[,i]<-rowSums((Y[,grp.mtch]-grp.mean[,i])^2)
    }
    
    df1<-(ngrps-1)
    df2<-sum(n)-df1-1
    
    SSW<-rowSums(grp.SSW)
    SSB<-SST-SSW
    MSE<-SSW/df2
    MSB<-SSB/df1
    
    F.stat<-MSB/MSE
    pval<-1-pf(F.stat,df1,df2)
    
    ebp<-grenander.ebp(unlist(pval))
    res<-cbind.data.frame(stat=F.stat,pval=pval,ebp=ebp$ebp)
    return(res)
  }

data_subjects_columns <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/Chamberlain_protein_abundance_data_subject_columns.xlsx"))

# Removing row names (gene symbols for each protein)
data_subjects_columns_exc_names <- data_subjects_columns[2:31]

one_way_anova_result <- row.oneway.anova(data_subjects_columns_exc_names,Chamberlain_proteomics_data_transposed[,1])

# Adding gene names back
one_way_anova_result$gene_name <- data_subjects_columns[,1]

# Adjusting p-value for multiple comparisons using Benjamini-Hochberg correction method
one_way_anova_result$adjpval <- p.adjust(one_way_anova_result$pval,method = "BH")

# Adding column to display significant and non-significant ANOVA results

one_way_anova_result$sig_result <- NA

one_way_anova_result$sig_result[one_way_anova_result$adjpval < 0.05] <- "Significant"
                                
one_way_anova_result$sig_result[one_way_anova_result$adjpval >= 0.05] <- "Not Significant"

Chamberlain_proteomics_data_with_anova <- Chamberlain_proteomics_data[,c("Description", "Accession","Gene.Symbol")]

Chamberlain_proteomics_data_with_anova$adjpval <-one_way_anova_result[,"adjpval"]

Chamberlain_proteomics_data_with_anova$sig_result <-one_way_anova_result[,"sig_result"]

# Remove rows that are not signficant

ANOVA_sig_only_proteinlist <- Chamberlain_proteomics_data_with_anova %>%
  dplyr::filter(.,!grepl("Not Significant",sig_result))

----
  
############################################################################

# Filtering 2-group DEP datasets by proteins significant by one-way ANOVA

mdx_vs_WT_ANOVA_filtered = mdx_vs_WT_DEPs[mdx_vs_WT_DEPs$Gene.Symbol %in% 
                                            ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(mdx_vs_WT_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/mdxvsWT_sigANOVA.xlsx",
#                      sheetName="DEPs sig by ANOVA and BH corr")

full_mdx_vs_WT_ANOVA_filtered <- Chamberlain_proteomics_data %>%
  filter(Gene.Symbol %in% mdx_vs_WT_ANOVA_filtered$Gene.Symbol)

#write.xlsx(full_mdx_vs_WT_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/full_mdxvsWT_sigANOVA.xlsx",
#                      sheetName="DEPs sig by ANOVA and BH corr")

uDysmdx_vs_WT_ANOVA_filtered = 
  uDysmdx_vs_WT_DEPs[uDysmdx_vs_WT_DEPs$Gene.Symbol %in% 
                       ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(uDysmdx_vs_WT_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/uDysmdxvsWT_sigANOVA.xlsx",
#           sheetName="DEPs sig by ANOVA and BH corr")

MidiDysmdx_vs_WT_ANOVA_filtered = 
  MidiDysmdx_vs_WT_DEPs[MidiDysmdx_vs_WT_DEPs$Gene.Symbol %in% 
                          ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(MidiDysmdx_vs_WT_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/MidiDysmdxvsWT_sigANOVA.xlsx",
#           sheetName="DEPs sig by ANOVA and BH corr")

flDysmdx_vs_WT_ANOVA_filtered = 
  flDysmdx_vs_WT_DEPs[flDysmdx_vs_WT_DEPs$Gene.Symbol %in% 
                        ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(flDysmdx_vs_WT_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/flDysmdxvsWT_sigANOVA.xlsx",
#           sheetName="DEPs sig by ANOVA and BH corr")

mdx_vs_uDysmdx_ANOVA_filtered = 
  mdx_vs_uDysmdx_DEPs[mdx_vs_uDysmdx_DEPs$Gene.Symbol %in% 
                        ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(mdx_vs_uDysmdx_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/mdxvsuDysmdx_sigANOVA.xlsx",
#           sheetName="DEPs sig by ANOVA and BH corr")

mdx_vs_MidiDysmdx_ANOVA_filtered = 
  mdx_vs_MidiDysmdx_DEPs[mdx_vs_MidiDysmdx_DEPs$Gene.Symbol %in% 
                           ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(mdx_vs_MidiDysmdx_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/mdxvsMidiDysmdx_sigANOVA.xlsx",
#           sheetName="DEPs sig by ANOVA and BH corr")

mdx_vs_flDysmdx_ANOVA_filtered = 
  mdx_vs_flDysmdx_DEPs[mdx_vs_flDysmdx_DEPs$Gene.Symbol %in% 
                           ANOVA_sig_only_proteinlist$Gene.Symbol, ]

#write.xlsx(mdx_vs_flDysmdx_ANOVA_filtered, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/mdxvsflDysmdx_sigANOVA.xlsx",
#          sheetName="DEPs sig by ANOVA and BH corr")

----
  
##########################################################################
# Filtering DEP lists by proteins overlapping between specific mdx groups
##########################################################################

overlap_mdxvsWT_and_uDysmdxvsWT <- 
  data.frame('overlap_mdxvsWT_and_uDysmdxvsWT' = intersect
             (mdx_vs_WT_DEPs$Gene.Symbol,
               uDysmdx_vs_WT_DEPs$Gene.Symbol)) %>% 
  'colnames<-' ("Gene symbol")

overlap_mdxvsWT_and_uDysmdxvsWT_and_midiDysmdx_vs_WT <-
  data.frame('overlap_mdxvsWT_and_uDysmdxvsWT_and_midiDysmdx_vs_WT' = intersect
             (overlap_mdxvsWT_and_uDysmdxvsWT$`Gene symbol`,
               MidiDysmdx_vs_WT_DEPs$Gene.Symbol)) %>% 
  'colnames<-' ("Gene symbol")

overlap_uDysmdxvsWT_and_flDysmdxvsWT <- 
  data.frame('overlap_uDysmdxvsWT_and_flDysmdxvsWT' = intersect
             (uDysmdx_vs_WT_DEPs$Gene.Symbol,
               flDysmdx_vs_WT_DEPs$Gene.Symbol)) %>% 
  'colnames<-' ("Gene symbol")

unique_overlap_uDysmdxvsWT_and_flDysmdxvsWT <-
  overlap_uDysmdxvsWT_and_flDysmdxvsWT %>% 
  base::subset(!(.$`Gene symbol` %in% 
                   overlap_mdxvsWT_and_uDysmdxvsWT_and_midiDysmdx_vs_WT$`Gene symbol`)) %>%
  as.data.frame(.) %>%
  base::subset(!(.$`Gene.Symbol` %in% mdx_vs_WT_DEPs$`Gene.Symbol`)) %>%
  as.data.frame(.)

overlap_MidiDysmdxvsWT_and_flDysmdxvsWT <- 
  data.frame('overlap_MidiDysmdxvsWT_and_flDysmdxvsWT' = intersect
             (MidiDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol,
               flDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol)) %>% 
  'colnames<-' ("Gene symbol")

overlap_mdxvsWT_and_MidiDysmdxvsWT_and_flDysmdx_vs_WT <-
  data.frame('overlap_mdxvsWT_and_MidiDysmdxvsWT_and_flDysmdx_vs_WT' = intersect
             (overlap_MidiDysmdxvsWT_and_flDysmdxvsWT$`Gene symbol`,
               mdx_vs_WT_ANOVA_filtered$Gene.Symbol)) %>% 
  'colnames<-' ("Gene symbol")

overlap_split_intein_constructs_and_mdx_no_uDysmdxvsWT <-
  overlap_mdxvsWT_and_MidiDysmdxvsWT_and_flDysmdx_vs_WT %>% 
  base::subset(!(.$`Gene symbol` %in% uDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol)) %>%
  as.data.frame(.) #%>%
  #base::subset(!(.$`Gene.Symbol` %in% mdx_vs_WT_DEPs$`Gene.Symbol`)) %>%
  #as.data.frame(.)

####################################################################

# Lists of unrestored proteins in mdx mice treated with AAV vectors

# uDysmdx AAV

unchanged_mdx_vs_uDysmdx_proteins <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsuDysmdx_unchanged.xlsx"))

unrestored_uDysmdx_proteins <- 
  data.frame('unrestored_uDysmdx_proteins' = intersect
             (uDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol,
               unchanged_mdx_vs_uDysmdx_proteins$Gene.Symbol)) %>%
  'colnames<-' ("Gene symbol")
#write.xlsx(unrestored_uDysmdx_proteins, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/unrestored_uDysmdx_proteins.xlsx",
#           sheetName="Protein list")

# MidiDysmdx AAV

unchanged_mdx_vs_MidiDysmdx_proteins <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsMidiDysmdx_unchanged.xlsx"))

unrestored_MidiDysmdx_proteins <- 
  data.frame('unrestored_MidiDysmdx_proteins' = intersect
             (MidiDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol,
               unchanged_mdx_vs_MidiDysmdx_proteins$Gene.Symbol)) %>%
  'colnames<-' ("Gene symbol")
#write.xlsx(unrestored_MidiDysmdx_proteins, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/unrestored_MidiDysmdx_proteins.xlsx",
#           sheetName="Protein list")

# flDysmdx AAV

unchanged_mdx_vs_flDysmdx_proteins <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/Data files for R/mdxvsflDysmdx_unchanged.xlsx"))

unrestored_flDysmdx_proteins <- 
  data.frame('unrestored_flDysmdx_proteins' = intersect
             (flDysmdx_vs_WT_ANOVA_filtered$Gene.Symbol,
               unchanged_mdx_vs_flDysmdx_proteins$Gene.Symbol)) %>%
  'colnames<-' ("Gene symbol")
#write.xlsx(unrestored_flDysmdx_proteins, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Proteomics Chamberlain TG mdx/unrestored_flDysmdx_proteins.xlsx",
#           sheetName="Protein list")
