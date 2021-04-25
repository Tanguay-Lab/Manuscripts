###########################################################################
##  Script and Data Analysis by Cheryl L Dunham                          ##
##                                                                       ##
##  Budding Bioinformatician                                             ##
##  Tanguay Laboratory | Sinnhuber Aquatic Research Laboratory (SARL)    ##
##  Environmental and Molecular Toxicology                               ##
##  Oregon State University                                              ##
##                                                                       ##
##  28645 E Hwy 34                                                       ##
##  Corvallis, OR 97333                                                  ##
##                                                                       ##
##  (541) 737-6500                                                       ##
##  dunhamcg@gmail.com                                                   ##
###########################################################################


##Use to analyse Zebrafish miRNA-seq data


##----1) Install any missing packages----
#Use to install multiple CRAN packages
#install.packages(c('xlsx','dplyr')) 
#To install a single CRAN package
#install.packages('xlsx') 

#To install Bioconductor packages
#install.packages('BiocManager')
#To install a single Bioconductor package
#BiocManager::install('biomaRt') 
#To install multiple Bioconductor packages
#BiocManager::install(c('biomaRt, edgeR')) 
#To upgrade current packages
#BiocManager::install() 
#To upgrade current packages to specific version of Bioconductor
#BiocManager::install(version=3.9) 

#To install GitHub packages
#devtools::install_github("AliciaSchep/iheatmapr") 


##----2) Load packages and specify file paths----
#List all required packages to load in loop
all_packages <- c('dplyr',
                  'xlsx',   #requires updated version of Java (most likely 64 bit)
                  'rJava',  #required to help control Java memory problems that arise from using xlsx package
                  'dichromat',  #colorblind friendly color palettes
                  'scales',  #enables you to see what color codes look like (makes a color strip image)
                  'ggplot2',
                  'ggrepel',
                  'dendextend',
                  ##Bioconductor packages:
                  'biomaRt',
                  'DESeq2',
                  'edgeR',
                  'Glimma',
                  'ComplexHeatmap')
#Load all packages at once
sapply(all_packages, require, character.only = TRUE) #require() similar to the function library() except it returns FALSE instead of a warning if a package fails to load

#I append this to the beginning of all files to ensure all file names go into whatever folder you specify
output_filepath <- 'C:/Users/dunhamc/Google Drive/miRNAs/FRCs/mirDeep2/Output_Mean_28_06_2019_t_15_44_06'
#File path of data to be analyzed, this also helps with documentation later on
input_filepath <- 'C:/Users/dunhamc/Google Drive/miRNAs/FRCs/miRDeep2/Mapper/Prediction_2012/miRNAs_expressed_all_samples_28_06_2019_t_15_44_06.tab'
#Label to append to file names
identifier <- '_FRC_miRNA'
#Label to add to images
appender <- ' (miRNA-seq)'
#Server file path indicating the location of trimmed reads
trimmer_filepath <- '/nfs1/EMT/Tanguay_Lab/Bioinfo/dunhamc/FRCs/miRNA/Skewer3/FA_Files/No_whitespace'
#Server file path indicating the location of aligned reads
mapper_filepath <- '/nfs1/EMT/Tanguay_Lab/Bioinfo/dunhamc/FRCs/miRNA/miRDeep2/Mapper'
#Server file path indicating the location of miRNA quantification files
quantifier_filepath <- '/nfs1/EMT/Tanguay_Lab/Bioinfo/dunhamc/FRCs/miRNA/miRDeep2/Mapper/Prediction_2012'


##----3) Import and format sequence data----
#Read in miRDeep2 expression table
##NOTE: The hashtag (#) must be removed from the front of the miRNA column of the spreadsheet for the table to import properly
experiment_counts <- read.table(input_filepath, sep ='\t', stringsAsFactors = FALSE, header = TRUE)

#Average counts from duplicated mature miRNAs
mean_counts <- aggregate(x = experiment_counts, by = list(experiment_counts$miRNA), FUN = mean)
#Use warnings() to check that warnings are only due to "returning NA" in columns 2 and 4
warnings()
#Remove unnecessary miRDeep2 text columns
mean_counts <- mean_counts[,-c(2,4)]
#Rename column 1
names(mean_counts)[1] <- 'miRNA'
#Write the table to file
write.table(mean_counts, file = paste0(output_filepath,'/Mean_count_table',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Label rows of the mean count table with miRBase IDs
row.names(mean_counts) <- mean_counts$miRNA
#Create save-compatible count matrix excluding normalized values (edgeR only works on raw count data)
mean_counts <- mean_counts[,c(1,4:(ncol(mean_counts)-((ncol(mean_counts)-3)/2)))]
#Write count matrix to file
write.table(mean_counts, file = paste0(output_filepath,'/Mean_count_matrix',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
#Create edgeR count matrix by removing text from the matrix
mean_counts <- mean_counts[,-1]


##----4) Add sample information and link sample info to count table----
sampleInfo <- read.delim(file = 'C:/Users/dunhamc/Google Drive/RNA_Seq/FRCs/sampleInfo.txt', header = TRUE, stringsAsFactors = FALSE)
#Add experimental groups/factors for differential testing
sampleInfo$Group <- factor(sampleInfo$Group)
#Relevel the variable so that the control group "Con" is the reference against which all other samples will be tested
sampleInfo$Group <- relevel(sampleInfo$Group, ref = 'Con')
sampleInfo$Lane <- factor(sampleInfo$Lane)
sampleInfo$Chemical <- factor(sampleInfo$Chemical)
#Relevel the variable so that the control group "DMSO" is the reference against which all other samples will be tested
sampleInfo$Chemical <- relevel(sampleInfo$Chemical, ref = 'DMSO')
#Label the data in sampleInfo with the abbreviated ID provided in the sample information sheet
row.names(sampleInfo) <- sampleInfo$subID
#Verify IDs match in the sampleInfo and mean_counts tables
sum(!(sampleInfo$subID %in% names(mean_counts)))
sum(!(names(mean_counts) %in% sampleInfo$subID))

#Ensure the column info from count data matches row info from sampleInfo sheet
sampleInfo <- sampleInfo[names(mean_counts),]
#Write sample information sheet to file for documentation
write.xlsx(x = sampleInfo, file = paste0(output_filepath,'/sampleInfo',identifier,'.xlsx'), row.names = FALSE)
#Write count matrix to file
write.table(experiment_counts, file = paste0(output_filepath,'/Count_matrix',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


##----5) Make the edgeR digital gene expression (DGE) data object----
DGE <- DGEList(counts = mean_counts, 
               genes = row.names(mean_counts), 
               group = sampleInfo$Group) #add more sample info
#Add additional sample information to the DGE object
DGE$samples$Chemical_ID <- sampleInfo$Chemical
DGE$samples$Chemical_Name <- sampleInfo$Chemical_Name
DGE$samples$subID <- sampleInfo$subID
DGE$samples$Lane <- sampleInfo$Lane


##----6) Additional filtering based on edgeR vignette and (Chen 2016)----
#Determine counts per million corresponding to a minimum of 10-20 read counts 
mean_cpm_gene_count_of_15 <- round(as.numeric(cpm(15, mean(DGE$samples$lib.size))),2)
#Determine the number of samples in the smallest group
num_samples_smallest_group <- min(table(sampleInfo$Group))

#Only keep miRNAs that have a minimum of 10-20 counts in at least 4 samples to be included in analysis
keep.exprs <- rowSums(cpm(DGE)> mean_cpm_gene_count_of_15) >= num_samples_smallest_group
DGE <- DGE[keep.exprs, , keep.lib.sizes=FALSE]

#Write background miRNAs to file
write.table(DGE$genes$genes, file = paste0(output_filepath,'/Zfish_background_miRNAs',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


##----7) Normalize libraries----
#Normalize libraries using the TMM method
DGE <- edgeR::calcNormFactors(DGE) #use cpm(DGE) to get normalized count table
#Write normalization documentation file
writeLines(capture.output({
    
    cat('Script filepath:', rstudioapi::getSourceEditorContext()$path, '\n\n')
    cat('Skewer filepath:', trimmer_filepath, '\n\n')
    cat('Mapper filepath', mapper_filepath, '\n\n')
    cat('miRDeep2 quantifier filepath: ', quantifier_filepath, '\n\n')
    cat('Input filepath:', input_filepath, '\n\n')
    cat('Output filepath:', output_filepath, '\n\n')
    cat('_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\n')
    cat('Count per million reads per sample for mean gene count of 15:', mean_cpm_gene_count_of_15,'\n\n')
    cat('Minimum number of biological replicates per group: ', num_samples_smallest_group,'\n\n')
    cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n')
    print(DGE$samples)
    
}), con = paste0(output_filepath,'/DGE_libsize_normfactors',identifier,'.txt'))


##----Create normalized count data for use with heatmaps and PCA plots----
#Calculate cpm to factor in TMM normalization
rltdata <- edgeR::cpm(y = DGE, normalized.lib.sizes = TRUE) #get normalized count data in counts per million units
#Multiply normalized cpm numbers by library size per million to get count values normalized by library size
for (bio_sample in colnames(rltdata)){
    rltdata[,bio_sample] <- (rltdata[,bio_sample])*((DGE$samples[bio_sample,'lib.size']/1000000))
}
#Assign TMM-normalized counts to the 'normalized_counts' slot of the DGE object
DGE$normalized_counts <- rltdata
#Regular log transformation of count data (corrects for heteroskedasticy)
rltdata <- DESeq2::rlogTransformation(round(rltdata)) 
#Add miRNA labels to TMM-normalized, rlog transformed data
row.names(rltdata) <- row.names(DGE$counts)
#Add rlog normalized data to DGE data object
DGE$rlog_transformed <- rltdata
#Write TMM-normalized, rlog transformed count table to file
write.table(cbind(row.names(DGE$rlog_transformed), DGE$rlog_transformed), file = paste0(output_filepath,'/DGE_rlog_TMM_normalized_counts',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----8) PCA and interactive MDS plot----
#Create directory for storing interactive plots
if(!dir.exists(paste0(output_filepath,'/Interactive_Plots'))){dir.create(paste0(output_filepath,'/Interactive_Plots'))}

#Create interactive MDS plot of all miRNAs
glMDSPlot(cpm(DGE, log = TRUE), top = Inf, 
          col = colors[DGE$samples$group], 
          pch = pch[DGE$samples$group], dim = c(1,2), 
          groups=DGE$samples[,c('group','Chemical_ID','Chemical_Name','Lane')],
          main = paste0('MDS plot of all genes', appender), 
          html = paste0('MDS_plot_allGenes',identifier),
          labels = DGE$samples$subID,
          launch = FALSE,
          path = output_filepath,
          folder = 'Interactive_Plots')

#Create interactive MDS plot of the top 100 most variable miRNAs using pairwise selection between each gene
glMDSPlot(cpm(DGE, log = TRUE), top = 100, 
          gene.selection = 'pairwise',
          col = colors[DGE$samples$group], 
          pch = pch[DGE$samples$group], dim = c(1,2), 
          groups=DGE$samples[,c('group','Chemical_ID','Chemical_Name','Lane')],
          main = paste0('MDS plot of 100 most variable genes (pairwise selection)', appender), 
          html = paste0('MDS_plot_100_pairwise',identifier),
          labels = DGE$samples$subID,
          launch = FALSE,
          path = output_filepath,
          folder = 'Interactive_Plots')

#Create interactive MDS plot of the top 100 most variable miRNAs using common gene selection among all genes
glMDSPlot(cpm(DGE, log = TRUE), top = 100, 
          gene.selection = 'common',
          col = colors[DGE$samples$group], 
          pch = pch[DGE$samples$group], dim = c(1,2), 
          groups=DGE$samples[,c('group','Chemical_ID','Chemical_Name','Lane')],
          main = paste0('MDS plot of 100 most variable genes (common selection) ', appender), 
          html = paste0('MDS_plot_100_common',identifier),
          labels = DGE$samples$subID,
          launch = FALSE,
          path = output_filepath,
          folder = 'Interactive_Plots')


##----Create PCA plots----
#Perform principle components analysis on the Tmm-normalized, rlog-transformed data
pca <- prcomp(t(rltdata)) #all data
#Extract principle component proportion of variance
sumpca<-data.frame(summary(pca)$importance)

#Make colorblind-friendly color palette for each chemical
cat_12 <- colorschemes$Categorical.12
color_palette <- cat_12[c(7,1:5,8:12)]
#To visualize the chosen colors
show_col(color_palette)

#Create PCA plot of all miRNAs
PCA_plot <- ggplot(data = data.frame(unlist(pca$x[,1:2])), 
                   aes(x = PC1, y = PC2)) + #graph first two PCs
    geom_point(aes(fill = factor(DGE$samples$group)),
               size = 3.5,
               shape = 21,
               stroke = 0.75) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    scale_fill_manual(values = color_palette) +
    xlab(paste('PC1 (', 
               round(sumpca['Proportion of Variance','PC1']*100),
               '%)', sep = '')) +
    ylab(paste('PC2 (',
               round(sumpca['Proportion of Variance','PC2']*100),
               '%)', sep = '')) +
    ggtitle(paste0('PCA plot of TMM-normalized, \nrlog-transformed genes',appender)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme(text = element_text(family = 'sans'),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key = element_rect(fill = NA, color = NA),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          axis.ticks = element_line(color = 'black'),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.border = element_blank(),
          panel.grid = element_blank())

PCA_plot

#Write the PCA plot to PDF
pdf(paste0(output_filepath,'/PCA_plot_all_genes',identifier,'.pdf'), width = 6, height = 4)
PCA_plot
dev.off()

#Create PCA plot with abbreviated IDs added to the plot
PCA_plot_ids <- ggplot(data = data.frame(unlist(pca$x[,1:2])), 
                   aes(x = PC1, y = PC2)) + #graph first two PCs
    geom_point(aes(fill = factor(DGE$samples$group)),
               size = 4.75,
               shape = 21,
               stroke = 1) +
    geom_text_repel(aes(label = DGE$samples$subID), 
                    size = 3, 
                    box.padding = unit(0.5,'line'),
                    point.padding = unit(0.5,'line')) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    scale_fill_manual(values = color_palette) +
    xlab(paste('PC1 (', 
               round(sumpca['Proportion of Variance','PC1']*100),
               '%)', sep = '')) +
    ylab(paste('PC2 (',
               round(sumpca['Proportion of Variance','PC2']*100),
               '%)', sep = '')) +
    ggtitle(paste0('PCA plot of TMM-normalized, \nrlog-transformed genes',appender)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme(text = element_text(family = 'sans'),
          plot.title = element_text(size = 19, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = NA, color = NA),
          legend.key.height = unit(25, units = 'pt'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.ticks = element_line(color = 'black'),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.border = element_blank(),
          panel.grid = element_blank())

PCA_plot_ids

#Write PCA plot of all genes to file
pdf(paste0(output_filepath,'/PCA_plot_all_genes_ids',identifier,'.pdf'), width = 8, height = 6)
PCA_plot_ids
dev.off()
    

##----Create PCA plots of most variable miRNAs----
#Identify the most variable genes as determined by row variance
topVarGenes <- order(rowVars(rltdata), decreasing = TRUE) 
#Create a table of TMM-normalized, rlog-transformed genes
gen_mat <- as.data.frame(rltdata)
#Reorder the gene table by the most variable genes
gen_mat <- gen_mat[topVarGenes,]
#Scale the transformed count data
gen_mat <- t(scale(t(gen_mat)))

#Perform PCA analysis on the top 100 miRNAs based on variance
pca_100 <- prcomp(t(gen_mat[1:100,]))
#Extract principle component proportion of variance for the top 100 genes
sumpca_100<-data.frame(summary(pca_100)$importance)

#Create PCA plot of top 100 miRNAs sorted by row variance
PCA_plot_100 <- ggplot(data = data.frame(unlist(pca_100$x[,1:2])), 
                   aes(x = PC1, y = PC2)) + #graph first two PCs
    geom_point(aes(fill = factor(DGE$samples$group)),
               size = 3.5,
               shape = 21,
               stroke = 0.75) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    scale_fill_manual(values = color_palette) +
    xlab(paste('PC1 (', 
               round(sumpca_100['Proportion of Variance','PC1']*100),
               '%)', sep = '')) +
    ylab(paste('PC2 (',
               round(sumpca_100['Proportion of Variance','PC2']*100),
               '%)', sep = '')) +
    ggtitle(paste0('PCA plot of TMM-normalized, rlog-transformed,\ntop 100 genes based on gene variance', appender)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme(text = element_text(family = 'sans'),
          plot.title = element_text(size = 14, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key = element_rect(fill = NA, color = NA),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          axis.ticks = element_line(color = 'black'),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.border = element_blank(),
          panel.grid = element_blank())
PCA_plot_100

#Write PCA plot of the top 500 most variable genes to PDF
pdf(paste0(output_filepath,'/PCA_plot_100mv',identifier,'.pdf'), width = 6, height = 4)
PCA_plot_100
dev.off()

#Create PCA plot of the top 100 most variable miRNAs with abbreviated sample IDs added to the plot
PCA_plot_100_ids <- ggplot(data = data.frame(unlist(pca_100$x[,1:2])), 
                       aes(x = PC1, y = PC2)) + #graph first two PCs
    geom_point(aes(fill = factor(DGE$samples$group)),
               size = 4.75,
               shape = 21,
               stroke = 1) +
    geom_text_repel(aes(label = DGE$samples$subID), 
                    size = 3, 
                    box.padding = unit(0.5,'line'),
                    point.padding = unit(0.5,'line')) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    scale_fill_manual(values = color_palette) +
    xlab(paste('PC1 (', 
               round(sumpca_100['Proportion of Variance','PC1']*100),
               '%)', sep = '')) +
    ylab(paste('PC2 (',
               round(sumpca_100['Proportion of Variance','PC2']*100),
               '%)', sep = '')) +
    ggtitle(paste0('PCA plot of TMM-normalized, rlog-transformed,\ntop 100 genes based on gene variance', appender)) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme(text = element_text(family = 'sans'),
          plot.title = element_text(size = 19, hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key = element_rect(fill = NA, color = NA),
          legend.key.height = unit(25, units = 'pt'),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          axis.ticks = element_line(color = 'black'),
          panel.background = element_rect(fill = 'white', color = 'black'),
          panel.border = element_blank(),
          panel.grid = element_blank())
PCA_plot_100_ids

#Write PCA plot of the top 500 most variable genes to PDF
pdf(paste0(output_filepath,'/PCA_plot_100mv_ids',identifier,'.pdf'), width = 8, height = 6)
PCA_plot_100_ids
dev.off()


##----9) Make the design matrix----
design <- model.matrix(~0 + DGE$samples$group) #the zero means that there is no intercept in the model, this helps with contrasts in later steps
#Assign sample names to the rows of the design matrix
rownames(design) <- colnames(DGE)
#Assign group names to the columns of the design matrix
colnames(design) <- levels(DGE$samples$group)
#Write the design matrix table to documentation file
write.table(cbind('Sample' = row.names(design), as.data.frame(design)), file = paste0(output_filepath,'/Design_matrix',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----10) Estimate dispersion----
#Estimate common, trended, and tagwise negative binomial dispersions by weighted likelihood empirical Bayes
DGE <- estimateDisp(DGE, design, robust = TRUE)
plotBCV(DGE, main = 'Dispersion estimation')


##----11) Fit data to the design model----
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data
fit <- glmQLFit(DGE, design, robust = TRUE)
plotQLDisp(glmfit = fit, main = 'Quasi-liklihood Dispersion')


##----12) Design contrast matrix to extract comparisons of interest----
my_contrasts <- makeContrasts(`BDE47-Con` = BDE47-Con,
                              `IPP-Con` = IPP-Con,
                              `TBBPA-Con` = TBBPA-Con,
                              `TBBPA_DBPE-Con` = TBBPA_DBPE-Con,
                              `TBPH-Con` = TBPH-Con,
                              `TCEP-Con` = TCEP-Con,
                              `TCPP-Con` = TCPP-Con,
                              `TDBPP-Con` = TDBPP-Con,
                              `TiBP-Con` = TiBP-Con,
                              `TPP-Con` = TPP-Con,
                              levels = design)
#Write contrast table to documentation file
write.table(cbind(data.frame('Group' = row.names(my_contrasts)), as.data.frame(my_contrasts)), file = paste0(output_filepath,'/Contrasts',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----13) Set up ComplexHeatmap parameters----
##----Set heatmap global parameters
ht_opt('heatmap_row_names_gp' = gpar(fontfamily = 'sans'),
       'heatmap_row_title_gp' = gpar(fontfamily = 'sans'),
       'legend_title_gp' = gpar(fontfamily = 'sans', fontsize = 11),
       'legend_labels_gp' = gpar(fontfamily = 'sans', fontsize = 9),
       'heatmap_column_names_gp' = gpar(fontfamily = 'sans'),
       'heatmap_column_title_gp' = gpar(fontfamily = 'sans'))#END ht_opt(RESET = TRUE)

#Set up column annotation information for treatment group
show_col(colorschemes$Categorical.12) #colorblind-friendly color palette
chemical_colors <- c('Con' = I(colorschemes$Categorical.12[7]),
                     'BDE47' = I(colorschemes$Categorical.12[1]),
                     'IPP' = I(colorschemes$Categorical.12[2]),
                     'TBBPA' = I(colorschemes$Categorical.12[3]),
                     'TBBPA_DBPE' = I(colorschemes$Categorical.12[4]),
                     'TBPH' = I(colorschemes$Categorical.12[5]),
                     'TCEP' = I(colorschemes$Categorical.12[8]),
                     'TCPP' = I(colorschemes$Categorical.12[9]),
                     'TDBPP' = I(colorschemes$Categorical.12[10]),
                     'TiBP' = I(colorschemes$Categorical.12[11]),
                     'TPP' = I(colorschemes$Categorical.12[12]))

#Create column annotation information object for 100 miRNAs
col_anno <- HeatmapAnnotation(df = data.frame(Chemical = DGE$samples$group),
                              name = 'Chemical',
                              col = list(Chemical = chemical_colors),
                              annotation_legend_param = list(
                                  Chemical = list(nrow = 1, 
                                                  title = "Chemical", 
                                                  title_gp = gpar(fontsize = 11, 
                                                                  fontface = 'plain'),
                                                  title_position = "topleft",
                                                  labels_gp = gpar(fontsize = 9, 
                                                                   fontfamily = 'sans'),
                                                  by_row = TRUE,
                                                  legend_direction = 'horizontal')),
                              show_legend = TRUE,
                              which = 'column',
                              # gp = gpar(col = 'black'),
                              show_annotation_name = TRUE,
                              annotation_name_gp = gpar(fontsize = 11),
                              annotation_name_side = 'left',
                              height = unit(3, "mm"),
                              simple_anno_size_adjust = TRUE
)#end column annotation

#Create column annotation information object for 30 miRNAs
col_anno30 <- HeatmapAnnotation(df = data.frame(Chemical = DGE$samples$group),
                                name = 'Chemical',
                                col = list(Chemical = chemical_colors),
                                annotation_legend_param = list(
                                    Chemical = list(nrow = 1, 
                                                    title = "Chemical", 
                                                    title_gp = gpar(fontsize = 8, 
                                                                    fontface = 'plain'),
                                                    title_position = "topleft",
                                                    labels_gp = gpar(fontsize = 6, 
                                                                     fontfamily = 'sans'),
                                                    by_row = TRUE,
                                                    legend_direction = 'horizontal')),
                                show_legend = TRUE,
                                which = 'column',
                                # gp = gpar(col = 'black'),
                                show_annotation_name = TRUE,
                                annotation_name_gp = gpar(fontsize = 8),
                                annotation_name_side = 'left',
                                height = unit(2.5, "mm"),
                                simple_anno_size_adjust = TRUE
)#end column annotation

#Create column annotation information object for 30 genes to be used in PowerPoint
col_anno30PPT <- HeatmapAnnotation(df = data.frame(Chemical = DGE$samples$group),
                                   name = 'Chemical',
                                   col = list(Chemical = chemical_colors),
                                   annotation_legend_param = list(
                                       Chemical = list(nrow = 1, 
                                                       title = "Chemical", 
                                                       title_gp = gpar(fontsize = 12.5, 
                                                                       fontface = 'plain'),
                                                       title_position = "topleft",
                                                       labels_gp = gpar(fontsize = 10.5, 
                                                                        fontfamily = 'sans'),
                                                       by_row = TRUE,
                                                       legend_direction = 'horizontal')),
                                   show_legend = TRUE,
                                   which = 'column',
                                   # gp = gpar(col = 'black'),
                                   show_annotation_name = TRUE,
                                   annotation_name_gp = gpar(fontsize = 12.5),
                                   annotation_name_side = 'left',
                                   height = unit(3.25, "mm"),
                                   simple_anno_size_adjust = TRUE
)#end column annotation


##----14) Find significant genes for each comparison with a fold change threshold of log2(1), or no fold change cutoff----
fold_change = 1
#Create table for recording all significantly expressed genes
all_sig_genes <- data.frame()

#Test for differential expression and create images and tables for each experimental group
for(comparison in colnames(my_contrasts)){
    Sys.sleep(.01)
    print(comparison)
    #Make a label for each contrast
    cutoff_contrast <- paste0(comparison, '_log2(',fold_change,')')
    
    ##Necessary to add FDR values of all genes to MDS and Volcano plots##
    #Test for differential expression
    comparisons <- edgeR::glmTreat(glmfit = fit, 
                                   contrast = my_contrasts[,comparison],
                                   lfc = 0)
    #Extract all differential expression values regardless of p-value or fold-chagne cutoffs
    sig_table <- data.frame(topTags(comparisons, n = Inf, 
                                    adjust.method = 'BH',
                                    p.value = 1))
    #Create negative log10(p-value) column for use in Volcano plots
    sig_table$negLogPval <- -log10(sig_table$FDR)
    #Order the table by gene ID
    sig_table <- sig_table[order(sig_table$genes),]
    FDR <- sig_table #necessary to add FDR values to MDS plots
    ##Above is necessary to add FDR values of all genes to MD and Volcano plots##
    
    ##----create output folder for each contrast----
    if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast))){dir.create(paste0(output_filepath,'/', cutoff_contrast))}
    if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))){dir.create(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))}
    #Write table of all differentially expressed genes without a fold-change cutoff to file
    write.table(x = FDR, file = paste0(output_filepath,'/',  cutoff_contrast,'/', cutoff_contrast, '_all_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
    
    #Extract the significant differentially expressed genes that pass the desired Benjamini-Hochberg FDR cutoff of 0.05
    sig_table <- data.frame(topTags(comparisons,
                                    n = Inf, adjust.method = "BH",
                                    p.value = 0.05))
    
    if(ncol(sig_table) > 0){
        #Create data objects for each DGELRT object created in the loop
        assign(cutoff_contrast, comparisons)
        #Add negative log10(FDR) values to the significant miRNA table
        sig_table$negLogPval <- -log10(sig_table$FDR)
        #Create data table for each table of significant differentially expressed miRNAs
        assign(paste(cutoff_contrast, 'sigGenes', sep = '_'), sig_table) 
        #Create data table for each table of significant differentially expressed miRNAs
        #Add significant differentially expressed miRNAs to the table of all differentially expressed miRNAs
        all_sig_genes <- rbind(all_sig_genes, sig_table) 
        
        #Write table of significant differentially expressed miRNAs to text file
        write.table(x = sig_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write table of significant differentially expressed miRNAs to Excel file
        xlsx::write.xlsx2(x = sig_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sig_log_ratios',identifier,'.xlsx'), row.names = FALSE)
        .jcall("java/lang/System", method = "gc") #helps with java memory problems
        #Write table of significant differentially expressed upregulated miRNAs to text file
        write.table(x = sig_table[sig_table$logFC > 0,], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'UR_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write table of significant differentially expressed downregulated miRNAs to text file
        write.table(x = sig_table[sig_table$logFC < 0,], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'DR_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        
        #Write list of significant differentially expressed miRNAs to text file
        write.table(x = sig_table$genes, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sigGene_list',identifier,'.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed upregulated miRNAs to text file
        write.table(x = sig_table$genes[sig_table$logFC > 0], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'UR_sigGene_list',identifier,'.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed downregulated miRNAs to text file
        write.table(x = sig_table$genes[sig_table$logFC < 0], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'DR_sigGene_list',identifier,'.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    
##----Create heatmaps of top fold changes and FDR----
    if(nrow(sig_table) >= 2){
        #Order miRNAs based on their absolute log2(fold change)
        heatmap_selector_FC <- sig_table[order(abs(sig_table$logFC), decreasing = TRUE),]
        #Order miRNAs based on their FDR-adjusted p-values
        heatmap_selector_FDR <- sig_table[order(sig_table$FDR),] 
        
        ##----Create heatmap of top 30 genes sorted by fold change----
        #Select the top 30 miRNAs
        heatmap_selector_FC30 <- heatmap_selector_FC[1:30,'genes']
        #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 miRNAs
        heatmap_selector_FC30 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FC30
        #Select data from the TMM-normalized, rlog transformed count table matching the top 30 miRNAs
        rlogdataFC30 <- data.frame(DGE$rlog_transformed[heatmap_selector_FC30,])
        #Scale the transformed count data
        rlogdataFC30 <- t(scale(t(rlogdataFC30)))
        
        #Set row dendrogram line thickness
        row_dend_FC = as.dendrogram(hclust(dist(rlogdataFC30))) %>% set('branches_lwd', 1.5)
        #Set column dendrogram line thickness
        col_dend_FC = as.dendrogram(hclust(dist(t(rlogdataFC30)))) %>% set('branches_lwd', 1.5)
        
        ##If there are 30 or fewer total significant differentially expressed genes, the following heatmaps are created
        if(nrow(heatmap_selector_FC) <= 30) {
            #Identify the number of significant differentially expressed miRNAs
            number_of_genes <- nrow(rlogdataFC30)
            ##Create main ComplexHeatmap layer of all significant differentially expressed miRNAs
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30_FC',
                            column_title = paste0(cutoff_contrast, ':\nAll ',number_of_genes,' differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 10),
                            row_dend_side = 'left',
                            row_dend_width = unit(10, 'mm'),
                            show_row_dend = TRUE,
                            column_dend_side = 'top',
                            column_dend_height = unit(7.5,'mm'),
                            show_column_dend = TRUE,
                            row_labels = rownames(rlogdataFC30),
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 6,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 8),
                            top_annotation = col_anno30,
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = 
                                list(color_bar = 'continuous',
                                     title = 'Z-score',
                                     title_gp = gpar(fontsize = 8, 
                                                     fontfamily = 'sans'),
                                     title_position = 'topleft', 
                                     labels_gp = gpar(fontsize = 6, 
                                                      fontfamily = 'sans'),
                                     legend_direction = 'horizontal')
            )#END hmap
            
            ##----Draw the  heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_all',number_of_genes,identifier,'.pdf'), width = 6, height = 5)
            draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##Create main ComplexHeatmap layer with dimensions for PowerPoint
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30_FC',
                            column_title = paste0(cutoff_contrast, ':\nAll ',number_of_genes,' differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14.5),
                            cluster_rows = row_dend_FC,
                            row_dend_side = 'left',
                            row_dend_width = unit(15, 'mm'),
                            show_row_dend = TRUE,
                            cluster_columns = col_dend_FC,
                            column_dend_side = 'top',
                            column_dend_height = unit(15,'mm'),
                            show_column_dend = TRUE,
                            row_labels = rownames(rlogdataFC30),
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 10.5,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 12.5),
                            top_annotation = col_anno30PPT,
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = 
                                list(color_bar = 'continuous',
                                     title = 'Z-score',
                                     title_gp = gpar(fontsize = 12.5, 
                                                     fontfamily = 'sans'),
                                     title_position = 'topleft', 
                                     labels_gp = gpar(fontsize = 10.5, 
                                                      fontfamily = 'sans'),
                                     legend_direction = 'horizontal')
            )#END hmap
            
            ##----Draw the heatmap and save as a PDF---
            #Draw the combined heatmap and write to file
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_all',number_of_genes,'_PPT',identifier,'.pdf'), width = 9, height = 7.5)
            draw(hmap, padding = unit(c(2,2,2,2), unit = 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
        }else{
            ##Create main ComplexHeatmap layer of the top 30 miRNAs sorted by absolute fold change
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30_FC',
                            column_title = paste0(cutoff_contrast, ': Top 30 fold change \nof all differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 10),
                            row_dend_side = 'left',
                            row_dend_width = unit(10, 'mm'),
                            show_row_dend = TRUE,
                            column_dend_side = 'top',
                            column_dend_height = unit(7.5,'mm'),
                            show_column_dend = TRUE,
                            row_labels = rownames(rlogdataFC30),
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 6,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 8),
                            top_annotation = col_anno30,
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = 
                                list(color_bar = 'continuous',
                                     title = 'Z-score',
                                     title_gp = gpar(fontsize = 8, 
                                                     fontfamily = 'sans'),
                                     title_position = 'topleft', 
                                     labels_gp = gpar(fontsize = 6, 
                                                      fontfamily = 'sans'),
                                     legend_direction = 'horizontal')
            )#END hmap
            
            ##----Draw the heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FC',identifier,'.pdf'), width = 6, height = 5)
            draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##Create main ComplexHeatmap layer with dimensions for PowerPoint for the top 30 miRNAs sorted by absolute fold change 
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30_FC',
                            column_title = paste0(cutoff_contrast, ': Top 30 fold change \nof all differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14.5),
                            cluster_rows = row_dend_FC,
                            row_dend_side = 'left',
                            row_dend_width = unit(15, 'mm'),
                            show_row_dend = TRUE,
                            cluster_columns = col_dend_FC,
                            column_dend_side = 'top',
                            column_dend_height = unit(15,'mm'),
                            show_column_dend = TRUE,
                            row_labels = rownames(rlogdataFC30),
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 10.5,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 12.5),
                            top_annotation = col_anno30PPT,
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = 
                                list(color_bar = 'continuous',
                                     title = 'Z-score',
                                     title_gp = gpar(fontsize = 12.5, 
                                                     fontfamily = 'sans'),
                                     title_position = 'topleft', 
                                     labels_gp = gpar(fontsize = 10.5, 
                                                      fontfamily = 'sans'),
                                     legend_direction = 'horizontal')
            )#END hmap
            
            ##----Draw the heatmap and save as a PDF---
            #Draw the combined heatmap and write to file
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FC_PPT',identifier,'.pdf'), width = 9, height = 7.5)
            draw(hmap, padding = unit(c(2,2,2,2), unit = 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##----Create heatmap of top 30 genes sorted by FDR-adjusted p-value----
            #Select the top 30 miRNAs
            heatmap_selector_FDR30 <- heatmap_selector_FDR[1:30,'genes']
            #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 miRNAs
            heatmap_selector_FDR30 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FDR30
            #Select data from the TMM-normalized, rlog transformed count table matching the top 30 miRNAs
            rlogdataFDR30 <- data.frame(DGE$rlog_transformed[heatmap_selector_FDR30,])
            #Scale the transformed count data
            rlogdataFDR30 <- t(scale(t(rlogdataFDR30)))
            
            #Set row dendrogram line thickness
            row_dend_FDR = as.dendrogram(hclust(dist(rlogdataFDR30))) %>% set('branches_lwd', 1.5)
            #Set column dendrogram line thickness
            col_dend_FDR = as.dendrogram(hclust(dist(t(rlogdataFDR30)))) %>% set('branches_lwd', 1.5)
            
            ##Create main ComplexHeatmap layer of the top 30 miRNAs sorted by FDR-adjusted p-value
            hmap <- Heatmap(matrix = rlogdataFDR30, 
                            name = 'top_30_FDR',
                            column_title = paste0(cutoff_contrast, ': Top 30 FDR-adjusted p-value \nof all differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 10),
                            row_dend_side = 'left',
                            row_dend_width = unit(10, 'mm'),
                            show_row_dend = TRUE,
                            column_dend_side = 'top',
                            column_dend_height = unit(7.5,'mm'),
                            show_column_dend = TRUE,
                            row_labels = rownames(rlogdataFDR30),
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 6,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 8),
                            top_annotation = col_anno30,
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = 
                                list(color_bar = 'continuous',
                                     title = 'Z-score',
                                     title_gp = gpar(fontsize = 8, 
                                                     fontfamily = 'sans'),
                                     title_position = 'topleft', 
                                     labels_gp = gpar(fontsize = 6, 
                                                      fontfamily = 'sans'),
                                     legend_direction = 'horizontal')
            )#END hmap
            
            ##----Draw the heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FDR',identifier,'.pdf'), width = 6, height = 5)
            draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##Create main ComplexHeatmap layer with dimensions for PowerPoint for the top 30 genes sorted by FDR-adjusted p-value
            hmap <- Heatmap(matrix = rlogdataFDR30, 
                            name = 'top_30_FDR',
                            column_title = paste0(cutoff_contrast, ': Top 30 FDR-adjusted p-value \nof all differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14.5),
                            cluster_rows = row_dend_FDR,
                            row_dend_side = 'left',
                            row_dend_width = unit(15, 'mm'),
                            show_row_dend = TRUE,
                            cluster_columns = col_dend_FDR,
                            column_dend_side = 'top',
                            column_dend_height = unit(15,'mm'),
                            show_column_dend = TRUE,
                            row_labels = rownames(rlogdataFDR30),
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 10.5,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 12.5),
                            top_annotation = col_anno30PPT,
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = 
                                list(color_bar = 'continuous',
                                     title = 'Z-score',
                                     title_gp = gpar(fontsize = 12.5, 
                                                     fontfamily = 'sans'),
                                     title_position = 'topleft', 
                                     labels_gp = gpar(fontsize = 10.5, 
                                                      fontfamily = 'sans'),
                                     legend_direction = 'horizontal')
            )#END hmap
            
            ##----Draw the heatmap and save as a PDF---
            #Draw the combined heatmap and write to file
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FDR_PPT',identifier,'.pdf'), width = 9, height = 7.5)
            draw(hmap, padding = unit(c(2,2,2,2), unit = 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##----Create heatmap of top 100 genes sorted by absolute fold change----
            #Select the top 100 miRNAs
            heatmap_selector_FC100 <- heatmap_selector_FC[1:100,'genes']
            #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 100 miRNAs
            heatmap_selector_FC100 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FC100
            #Select data from the TMM-normalized, rlog transformed count table matching the top 100 miRNAs
            rlogdataFC100 <- data.frame(DGE$rlog_transformed[heatmap_selector_FC100,])
            #Scale the transformed count data
            rlogdataFC100 <- t(scale(t(rlogdataFC100)))
            
            ##If there are 100 or fewer total significant differentially expressed miRNAs, the following heatmaps are created
            if(nrow(heatmap_selector_FC) <= 100){
                #Identify the number of significant differentially expressed genes
                number_of_genes <- nrow(rlogdataFC100)
                ##Create main ComplexHeatmap layer of all significant differentially expressed miRNAs
                hmap <- Heatmap(matrix = as.matrix(rlogdataFC100), 
                                name = 'all_de',
                                column_title = paste0(cutoff_contrast, ':\nAll ',number_of_genes,' differentially expressed genes',appender), 
                                column_title_side = 'top',
                                column_title_gp = gpar(fontsize = 14),
                                row_dend_side = 'left',
                                row_dend_width = unit(13, "mm"),
                                show_row_dend = TRUE,
                                column_dend_side = 'top',
                                column_dend_height = unit(8,'mm'),
                                show_column_dend = TRUE,
                                row_labels = rownames(rlogdataFC100),
                                row_names_side = 'right',
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 5,
                                                    fontface = 'italic'),
                                column_names_gp = gpar(fontsize = 11),
                                top_annotation = col_anno,
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = 
                                    list(color_bar = 'continuous',
                                         title = 'Z-score',
                                         title_gp = gpar(fontsize = 11, 
                                                         fontfamily = 'sans'),
                                         title_position = 'topleft', 
                                         labels_gp = gpar(fontsize = 9, 
                                                          fontfamily = 'sans'),
                                         legend_direction = 'horizontal')
                )#END hmap
                
                ##----Draw the heatmap and save as a PDF---
                pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_all',number_of_genes,identifier,'.pdf'), width = 8.5, height = 11)
                draw(hmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
                dev.off()
                
            }else{
                ##Create main ComplexHeatmap layer of the top 100 genes sorted by absolute fold change
                hmap <- Heatmap(matrix = as.matrix(rlogdataFC100), 
                                name = 'top_100_fc',
                                column_title = paste0(cutoff_contrast, ': Top 100 fold change \nof all differentially expressed genes',appender), 
                                column_title_side = 'top',
                                column_title_gp = gpar(fontsize = 14),
                                row_dend_side = 'left',
                                row_dend_width = unit(13, "mm"),
                                show_row_dend = TRUE,
                                column_dend_side = 'top',
                                column_dend_height = unit(8,'mm'),
                                show_column_dend = TRUE,
                                row_labels = rownames(rlogdataFC100),
                                row_names_side = 'right',
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 5,
                                                    fontface = 'italic'),
                                column_names_gp = gpar(fontsize = 11),
                                top_annotation = col_anno,
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = 
                                    list(color_bar = 'continuous',
                                         title = 'Z-score',
                                         title_gp = gpar(fontsize = 11, 
                                                         fontfamily = 'sans'),
                                         title_position = 'topleft', 
                                         labels_gp = gpar(fontsize = 9, 
                                                          fontfamily = 'sans'),
                                         legend_direction = 'horizontal')
                )#END hmap
                
                ##----Draw the heatmap and save as a PDF---
                pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top100FC',identifier,'.pdf'), width = 8.5, height = 11)
                draw(hmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
                dev.off()
                
                ##----Create heatmap of top 100 miRNAs sorted by FDR-adjusted p-value----
                #Select the top 100 miRNAs
                heatmap_selector_FDR100 <- heatmap_selector_FDR[1:100,'genes']
                #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 100 miRNAs
                heatmap_selector_FDR100 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FDR100
                #Select data from the TMM-normalized, rlog transformed count table matching the top 100 miRNAs
                rlogdataFDR100 <- data.frame(DGE$rlog_transformed[heatmap_selector_FDR100,])
                #Scale the transformed count data
                rlogdataFDR100 <- t(scale(t(rlogdataFDR100)))
                
                ##Create main ComplexHeatmap layer of the top 100 genes sorted by FDR-adjusted p-value
                hmap <- Heatmap(matrix = as.matrix(rlogdataFDR100), 
                                name = 'top_100_fc',
                                column_title = paste0(cutoff_contrast, ': Top 100 FDR-adjusted p-value \nof all differentially expressed genes',appender), 
                                column_title_side = 'top',
                                column_title_gp = gpar(fontsize = 14),
                                row_dend_side = 'left',
                                row_dend_width = unit(13, "mm"),
                                show_row_dend = TRUE,
                                column_dend_side = 'top',
                                column_dend_height = unit(8,'mm'),
                                show_column_dend = TRUE,
                                row_labels = rownames(rlogdataFDR100),
                                row_names_side = 'right',
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 5,
                                                    fontface = 'italic'),
                                column_names_gp = gpar(fontsize = 11),
                                top_annotation = col_anno,
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = 
                                    list(color_bar = 'continuous',
                                         title = 'Z-score',
                                         title_gp = gpar(fontsize = 11, 
                                                         fontfamily = 'sans'),
                                         title_position = 'topleft', 
                                         labels_gp = gpar(fontsize = 9, 
                                                          fontfamily = 'sans'),
                                         legend_direction = 'horizontal')
                )#END hmap
                
                ##----Draw the heatmap and save as a PDF---
                pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top100FDR',identifier,'.pdf'), width = 8.5, height = 11)
                draw(hmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
                dev.off()
            }
        }
        
        ##----Generate MD plots of differentially expressed genes----
        #Create directory for Interactive Glimma plots
        if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))){dir.create(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))}
        #Determine which genes are considered significant based on an FDR-adjusted p-value cutoff of 0.05
        dt <- decideTestsDGE(comparisons)
        #Add GeneID values to the differential expression table 
        comparisons$table$GeneID <- row.names(comparisons$table)
        #Add an miRNA miRBase identifier column to the comparisons differential expression table
        comparisons$table$miRNA <- comparisons$genes$miRNA
        #Add an miRNA miRBase precursor identifier column to the comparisons differential expression table
        comparisons$table$precursor <- comparisons$genes$precursor
        #Add the FDR-adjusted p-values to the differential expression table
        comparisons$table$FDR <- FDR$FDR
        #Add the -log10(FDR) to the differential expression table
        comparisons$table$negLogPval <- FDR$negLogPval
        #Add the p-values to the differential expression table
        comparisons$table$pValue <- FDR$PValue
        

        #Create interactive mean-difference (MD) plot comparing the mean log2(CPM) to the log2(fold change)
        glXYPlot(x = comparisons$table$logCPM,
                 y = comparisons$table$logFC,
                 xlab = 'log2(counts per million)',
                 ylab = 'log2(fold change)',
                 counts = DGE$normalized_counts,
                 groups = DGE$samples$group,
                 samples = DGE$samples$subID,
                 status = dt,
                 transform = TRUE,
                 launch = FALSE,
                 anno = comparisons$table,
                 display.columns = c('GeneID','pValue', 'FDR', 'logFC','logCPM','negLogPval'),
                 path = paste0(output_filepath,'/', cutoff_contrast),
                 folder = 'Interactive_Plots',
                 html = paste0('MD_Plot_', cutoff_contrast,identifier),
                 main = paste0(comparison, ': (', comparisons$comparison, ')<br> fold change cutoff = log<sub>2</sub>(',fold_change,')',appender))
        .jcall("java/lang/System", method = "gc")
        
        #Create interactive Volcano plot comparing the log2(fold change) and the -log10(FDR-adjusted p-value)
        glXYPlot(x = comparisons$table$logFC,
                 y = comparisons$table$negLogPval,
                 xlab = 'log2(fold change)',
                 ylab = '-log10(p-value)',
                 counts = DGE$normalized_counts,
                 groups = DGE$samples$group,
                 samples = DGE$samples$subID,
                 status = dt,
                 transform = TRUE,
                 launch = FALSE,
                 anno = comparisons$table,
                 display.columns = c('GeneID','pValue', 'FDR', 'logFC','logCPM','negLogPval'),
                 path = paste0(output_filepath,'/', cutoff_contrast),
                 folder = 'Interactive_Plots',
                 html = paste0('Volcano-Plot_', cutoff_contrast, identifier),
                 main = paste0(comparison, ': (', comparisons$comparison, ')<br> fold change cutoff = log<sub>2</sub>(',fold_change,')',appender))
        .jcall("java/lang/System", method = "gc")
        
    }
}


##----15) Heatmaps of significant miRNAs----
#Order all significantly differentially expressed genes by their FDR-adjusted p-value
all_sig_genesFDR <- all_sig_genes[order(all_sig_genes$FDR, decreasing = FALSE),]
#Remove duplicated values from the significant gene table
all_sig_genesFDR <- all_sig_genesFDR[!duplicated(all_sig_genesFDR$genes),'genes']

#Select the combined significantly differentially expressed miRNAs
heatmap_selector <- row.names(DGE$rlog_transformed) %in% all_sig_genesFDR
#Select data from the TMM-normalized, rlog transformed count table matching the combined significantly differentially expressed miRNAs
rlogdata <- data.frame(DGE$rlog_transformed[heatmap_selector,])
#Scale the transformed count data
rlogdata <- t(scale(t(rlogdata)))
#Identify the number of combined significant differentially expressed miRNAs
number_of_genes <- nrow(rlogdata)

##Create main ComplexHeatmap layer of all combined significant differentially expressed genes
hmap <- Heatmap(matrix = as.matrix(rlogdata), 
                name = 'all_miRNAs',
                column_title = paste0('All ',number_of_genes,' combined differentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 14),
                row_dend_side = 'left',
                row_dend_width = unit(13, "mm"),
                show_row_dend = TRUE,
                column_dend_side = 'top',
                column_dend_height = unit(8,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdata),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 1.5,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 11),
                top_annotation = col_anno,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 11, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 9),
                         legend_direction = 'horizontal')
                
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath, '/Heatmap_allDE_genes',identifier,'.pdf'), width = 8.5, height = 11)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()

##----Create heatmap of top 100 combined significantly differentially expressed miRNAs sorted by FDR-adjusted p-value----
#Determine the data from the TMM-normalized, rlog transformed count table that matches the top 100 combined significantly differentially expressed miRNAs
heatmap_selector_FDR <- row.names(DGE$rlog_transformed) %in% all_sig_genesFDR[1:100]
#Select data from the TMM-normalized, rlog transformed count table matching the top 100 combined significantly differentially expressed miRNAs
rlogdataFDR <- data.frame(DGE$rlog_transformed[heatmap_selector_FDR,])
#Scale the transformed count data
rlogdataFDR <- t(scale(t(rlogdataFDR)))

##Create main ComplexHeatmap layer of the top 100 combined significantly differentially expressed miRNAs sorted by FDR-adjusted p-value
hmap <- Heatmap(matrix = rlogdataFDR, 
                name = 'top_100_FDR',
                column_title = paste0('Top 100 FDR-adjusted p-value of combined\ndifferentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 14),
                row_dend_side = 'left',
                row_dend_width = unit(13, 'mm'),
                show_row_dend = TRUE,
                column_dend_side = 'top',
                column_dend_height = unit(8,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdataFDR),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 5,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 11),
                top_annotation = col_anno,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 11, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 9),
                         legend_direction = 'horizontal')
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath, '/Heatmap_top100FDR_DE_genes',identifier,'.pdf'), width = 8.5, height = 11)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()


##----Create heatmap of top 30 combined significantly differentially expressed miRNAs sorted by FDR-adjusted p-value----
#Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 combined significantly differentially expressed miRNAs
heatmap_selector_FDR <- row.names(DGE$rlog_transformed) %in% all_sig_genesFDR[1:30]
#Select data from the TMM-normalized, rlog transformed count table matching the top 30 combined significantly differentially expressed miRNAs
rlogdataFDR <- data.frame(DGE$rlog_transformed[heatmap_selector_FDR,])
#Scale the transformed count data
rlogdataFDR <- t(scale(t(rlogdataFDR)))
#Set row dendrogram line thickness
row_dend_FDR = as.dendrogram(hclust(dist(rlogdataFDR))) %>% set('branches_lwd', 1.5)
#Set column dendrogram line thickness
col_dend_FDR = as.dendrogram(hclust(dist(t(rlogdataFDR)))) %>% set('branches_lwd', 1.5)

##Create main ComplexHeatmap layer of the top 30 combined significantly differentially expressed miRNAs sorted by FDR-adjusted p-value
hmap <- Heatmap(matrix = rlogdataFDR, 
                name = 'top_30_FDR',
                column_title = paste0('Top 30 FDR-adjusted p-value of combined\ndifferentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 10),
                row_dend_side = 'left',
                row_dend_width = unit(10, 'mm'),
                show_row_dend = TRUE,
                column_dend_side = 'top',
                column_dend_height = unit(7.5,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdataFDR),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 6,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 8),
                top_annotation = col_anno30,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 8, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 8),
                         legend_direction = 'horizontal')
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath,'/Heatmap_top30FDR_DE_genes',identifier,'.pdf'), width = 6, height = 5)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()

##Create main ComplexHeatmap layer with dimensions for PowerPoint for the top 30 combined significantly differentially expressed miRNAs sorted by FDR-adjusted p-value
hmap <- Heatmap(matrix = rlogdataFDR, 
                name = 'top_30_FDR',
                column_title = paste0('Top 30 FDR-adjusted p-value of combined\ndifferentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 14.5),
                cluster_rows = row_dend_FDR,
                row_dend_side = 'left',
                row_dend_width = unit(15, 'mm'),
                show_row_dend = TRUE,
                cluster_columns = col_dend_FDR,
                column_dend_side = 'top',
                column_dend_height = unit(15,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdataFDR),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 10.5,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 12.5),
                top_annotation = col_anno30PPT,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 12.5, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 10.5),
                         legend_direction = 'horizontal')
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath,'/Heatmap_top30FDR_PPT_DE_genes',identifier,'.pdf'), width = 9, height = 7.5)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()


##----16) Heatmaps of all significant miRNAs: FC----
#Order miRNAs based on their absolute log2(fold change)
all_sig_genesFC <- all_sig_genes[order(abs(all_sig_genes$logFC), decreasing = TRUE),]
#Remove duplicated values from the significant gene table
all_sig_genesFC <- all_sig_genesFC[!duplicated(all_sig_genesFC$genes),'genes']

##----Create heatmap of top 100 combined significantly differentially expressed miRNAs sorted by absolute fold change----
#Determine the data from the TMM-normalized, rlog transformed count table that matches the top 100 combined significantly differentially expressed miRNAs
heatmap_selector_FC <- row.names(DGE$rlog_transformed) %in% all_sig_genesFC[1:100]
#Select data from the TMM-normalized, rlog transformed count table matching the top 100 combined significantly differentially expressed miRNAs
rlogdataFC <- data.frame(DGE$rlog_transformed[heatmap_selector_FC,])
#Scale the transformed count data
rlogdataFC <- t(scale(t(rlogdataFC)))

##Create main ComplexHeatmap layer of the top 100 combined significantly differentially expressed miRNAs sorted by absolute fold change
hmap <- Heatmap(matrix = rlogdataFC, 
                name = 'top_100_FC',
                column_title = paste0('Top 100 fold change of all combined\ndifferentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 14),
                row_dend_side = 'left',
                row_dend_width = unit(13, 'mm'),
                show_row_dend = TRUE,
                column_dend_side = 'top',
                column_dend_height = unit(8,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdataFC),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 5,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 11),
                top_annotation = col_anno,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 11, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 9),
                         legend_direction = 'horizontal')
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath, '/Heatmap_top100FC_DE_genes',identifier,'.pdf'), width = 8.5, height = 11)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()


##----Create heatmap of top 30 combined significantly differentially expressed miRNAs sorted by absolute fold change----
#Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 combined significantly differentially expressed miRNAs
heatmap_selector_FC <- row.names(DGE$rlog_transformed) %in% all_sig_genesFC[1:30]
#Select data from the TMM-normalized, rlog transformed count table matching the top 30 combined significantly differentially expressed miRNAs
rlogdataFC <- data.frame(DGE$rlog_transformed[heatmap_selector_FC,])
#Scale the transformed count data
rlogdataFC <- t(scale(t(rlogdataFC)))
#Set row dendrogram line thickness
row_dend_FC = as.dendrogram(hclust(dist(rlogdataFC))) %>% set('branches_lwd', 1.5)
#Set column dendrogram line thickness
col_dend_FC = as.dendrogram(hclust(dist(t(rlogdataFC)))) %>% set('branches_lwd', 1.5)

##Create main ComplexHeatmap layer of the top 30 combined significantly differentially expressed miRNAs sorted by absolute fold change
hmap <- Heatmap(matrix = rlogdataFC, 
                name = 'top_30_FC',
                column_title = paste0('Top 30 fold change of all combined\ndifferentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 10),
                row_dend_side = 'left',
                row_dend_width = unit(10, 'mm'),
                show_row_dend = TRUE,
                column_dend_side = 'top',
                column_dend_height = unit(7.5,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdataFC),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 6,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 8),
                top_annotation = col_anno30,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 8, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 6),
                         legend_direction = 'horizontal')
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath,'/Heatmap_top30FC_DE_genes',identifier,'.pdf'), width = 6, height = 5)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()

##Create main ComplexHeatmap layer with dimensions for PowerPoint for the top 30 combined significantly differentially expressed miRNAs sorted by absolute fold change 
hmap <- Heatmap(matrix = rlogdataFC, 
                name = 'top_30_FC',
                column_title = paste0('Top 30 fold change of all combined\ndifferentially expressed genes',appender), 
                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 14.5),
                cluster_rows = row_dend_FDR,
                row_dend_side = 'left',
                row_dend_width = unit(15, 'mm'),
                show_row_dend = TRUE,
                cluster_columns = col_dend_FDR,
                column_dend_side = 'top',
                column_dend_height = unit(15,'mm'),
                show_column_dend = TRUE,
                row_labels = rownames(rlogdataFC),
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 10.5,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 12.5),
                top_annotation = col_anno30PPT,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = 
                    list(color_bar = 'continuous',
                         title = 'Z-score',
                         title_gp = gpar(fontsize = 12.5, 
                                         fontface = 'plain'),
                         title_position = 'topleft', 
                         labels_gp = gpar(fontsize = 10.5),
                         legend_direction = 'horizontal')
)#END hmap

##----Draw the heatmap and save as a PDF---
pdf(file = paste0(output_filepath,'/Heatmap_top30FC_PPT_DE_genes',identifier,'.pdf'), width = 9, height = 7.5)
draw(hmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()


##----17) Write session information----
name <- 'miRNAseq'
writeLines(capture.output({
    cat('###########################################################################\n##  Script and Data Analysis by Cheryl L Dunham                          ##\n##                                                                       ##\n##  Budding Bioinformatician                                             ##\n##  Tanguay Laboratory | Sinnhuber Aquatic Research Laboratory (SARL)    ##\n##  Environmental and Molecular Toxicology                                   ##\n##  Oregon State University                                              ##\n##                                                                       ##\n##  28645 E Hwy 34                                                       ##\n##  Corvallis, OR 97333                                                  ##\n##                                                                       ##\n##  (541) 737-6500                                                       ##\n##  dunhamcg@gmail.com                                                   ##\n###########################################################################\n\n\n')
    cat('Script filepath:', rstudioapi::getSourceEditorContext()$path, '\n\n')
    cat('Input filepath:', input_filepath, '\n\n')
    cat('Output filepath:', output_filepath, '\n\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\n')
    print(sessionInfo())
    cat('\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nRStudio Version Information\n\n')
    print(RStudio.Version())
    cat('_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nAttached Package Citation Information\n\n')
    print(lapply(all_packages, citation))
    
}), paste0(output_filepath,'/', name,'_sessionInfo_',Sys.Date(),'.txt'))


#Save a copy of the R Global Environment to work with later
save.image(paste0(output_filepath,'/',name,'_',Sys.Date(),'.RData'))


#Remove objects from your global environment (can reload them from the .RData file created above)
rm(list = ls())