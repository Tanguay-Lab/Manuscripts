###########################################################################
##  Script and Data Analysis by Cheryl L Dunham                          ##
##                                                                       ##
##  Budding Bioinformatician                                             ##
##  Tanguay Laboratory | Sinnhuber Aquatic Research Laboratory (SARL)    ##
##  Environmental and Molecular Toxicology                               ##
##  Oregon State University                                              ##
##                                                                       ##
##  28645 E Hwy 34                                                       ##
##  Corvallis, OR 97333                                                  ##                                                                     ##
##  (541) 737-6500                                                       ##
##  dunhamcg@gmail.com                                                   ##
###########################################################################



##----1) Install any missing packages----
#Use to install multiple CRAN packages
#install.packages(c('xlsx','dplyr')) 
#Use to install a single CRAN package
#install.packages('xlsx') 

#To install Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite('edgeR') 
#Use to upgrade current Bioconductor packages
#biocLite()


##----2) Load packages and specify file paths----
#List all required packages to load in loop
all_packages <- c('dplyr', 
                  'xlsx',   #requires updated version of Java (most likely 64 bit)
                  'rJava',  #required to help control Java memory problems that arise from using xlsx package
                  'dichromat',  #color blind friendly color palettes
                  'scales',  #enables you to see what color codes look like (makes a color strip image)
                  'ggplot2',  
                  #Bioconductor packages:
                  'biomaRt',
                  'DESeq2',
                  'edgeR',
                  'ComplexHeatmap')
#Load all packages at once
sapply(all_packages, require, character.only = TRUE)

#I append this to the beginning of all files to ensure all file names go into whatever folder you specify
output_filepath <- 'C:/Users/dunhamc/Dropbox/SlincR (2nd paper)/Reviewer Comments/Resubmission/Data'
#File path of data to be analyzed
input_filepath <- 'C:/Users/dunhamc/Google Drive/RNA_Seq/Sox9blncRNA_RNA-Seq/Data/Zv10_89/HTseq_count'


##----3) Import and combine sequence data----
#Set the working directory to the folder containing HTSeq-count output .tab files
setwd(input_filepath)

HTseq_combo <- function(){
    #Combines data from multiple HTseq-count output files (sans the --additional-attr option) into a single dataframe
    #be sure to be in the correct working directory before running the function 
    
    #List all file names from the current working folder that end in .tab file extension
    filenames <- list.files(pattern = '.tab')
    #Read in the count table from the first file in the folder
    initial_count <- read.table(filenames[1], sep ='\t', stringsAsFactors = FALSE, header = FALSE)
    #Creates the dataframe to hold information from all .tab files in the current working folder
    counts_combo <- as.data.frame(matrix(ncol = length(filenames) + 1, nrow = length(initial_count[,1])))
    #Sets the first column of the combined dataframe to the Ensembl gene ID (the first column of HTSeq-count files)
    counts_combo[,1] <- initial_count[,1]
    #Sets the second column of the combined dataframe to the corresponding gene count values for the first .tab file in the current folder
    counts_combo[,2] <- initial_count[,2]
    #Sets the column names of the first 2 columns in the combined dataframe to 'Ensembl_ID' and the file name of the first .tab file in the current folder
    names(counts_combo)[1:2] <- c('Ensembl_ID',filenames[1])
    
    #Create loop to add gene count information to the combined dataframe from the rest of the .tab files in the current folder
    for (i in 2:length(filenames)){
        #Read in HTSeq-count data table
        counts <- read.table(filenames[i], sep = '\t', stringsAsFactors = FALSE, header = FALSE)
        #Add gene count information to the combined dataframe
        counts_combo[,i+1] <- counts[,2]
        #Link the HTSeq-count file name to its corresponding gene counts 
        names(counts_combo)[i+1] <- filenames[i]
    }
    
    #Remove non-count-data information from the combined data table (last 5 rows of the dataframe)
    counts_combo <- counts_combo[1:(length(counts_combo[,1])-5),]
    #Label rows of the combined data table with the corresponding Ensembl gene ID
    row.names(counts_combo) <- counts_combo[,1]
    #Remove the first column of the combined dataframe, leaving only count information in the data table
    counts_combo <- counts_combo[,-1]
    return(counts_combo)
}

#Merge count data into one data table
merged_counts <- HTseq_combo()

#Set the working directory
setwd(output_filepath)


##----4) Add sample information and link sample info to count table----
sampleInfo <- read.xlsx(file = 'C:/Users/dunhamc/Google Drive/RNA_Seq/Sox9blncRNA_RNA-Seq/Data/Zv10_89/Output/sampleInfo.xlsx', sheetIndex = 'Sheet1', header = TRUE, stringsAsFactors = FALSE)

#Add experimental groups/factors for differential testing
sampleInfo$Morphant <- factor(sampleInfo$Morphant, levels = c('CM','SM'))
sampleInfo$Treatment <- factor(sampleInfo$Treatment, levels = c('DMSO','TCDD'))
sampleInfo$Group <- factor(sampleInfo$Group)
#Ensure the 'Filename' column is of the character class and not a factor
sampleInfo$Filename <- as.character(sampleInfo$Filename)

#Verify that the file names associated with the merged count table matches with the filenames provided in the sampleInfo spreadsheet
colnames(merged_counts) %in% sampleInfo$Filename
sampleInfo$Filename %in% colnames(merged_counts)

#Create gene count matrix arranged by experimental sample order from sample information spreadsheet
experiment_counts <- merged_counts[,sampleInfo$Filename]
#Verify the column info from count data matches row info from sampleInfo sheet
colnames(experiment_counts) == sampleInfo$Filename
sum(colnames(experiment_counts) != sampleInfo$Filename)
#Rename experimental count table with corresponding sample ID abbreviation
colnames(experiment_counts) <- sampleInfo$subID
#Write count matrix to file
write.table(experiment_counts, file = paste0(output_filepath,'/Count_matrix.txt'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


##----5) Make the edgeR digital gene expression (DGE) data object----
#See the Bioconductor workflow, 'From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline' <https://www.bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html>, for an up-to-date RNAseq differential expression workflow
DGE <- DGEList(counts = experiment_counts, 
               genes = row.names(experiment_counts), 
               group = sampleInfo$Group) 
#Add additional sample information to the DGE object
DGE$samples$Treatment <- sampleInfo$Treatment
DGE$samples$Morphant <- sampleInfo$Morphant
DGE$samples$ID <- sampleInfo$ID
DGE$samples$subID <- sampleInfo$subID


##----6) Add gene annotation information----
#Zv10 (GRCz10) Ensembl release 89 annotation information
ensemblDrer89 <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl', host = 'may2017.archive.ensembl.org') 

##--Biomart options--
#sets the data columns to be added to the new dataframe (attributes) and filters the data based on user specifications and associated values (from the main Ensembl Biomart site)
#Filters
#listFilters(ensemblMmus89)  #to see all possible filter names
#Filters for use in all annotation retrieval
filts <- 'ensembl_gene_id' 
vals <- row.names(DGE) #values corresponding to filters (filts)

#Attributes  
#listAttributes(ensemblMmus89)  #to see all possible attribute names
#Attributes for zebrafish annotation information
zebrafish_attributes <- c('ensembl_gene_id','external_gene_name','gene_biotype')
#Attributes for zebrafish to human ortholog information
human_attributes <- c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name')
#Attributes for zebrafish to mouse ortholog information
mouse_attributes <- c('ensembl_gene_id','mmusculus_homolog_ensembl_gene','mmusculus_homolog_associated_gene_name')

##----Add zebrafish gene annotation information----
#Retrieve zebrafish Ensembl gene ID, gene name/symbol, and gene type for those genes included in the DGE object
martTable1 <- getBM(attributes = zebrafish_attributes, filters = filts, values = vals, mart = ensemblDrer89)

#Ensure the table is ordered by Ensembl gene ID
martTable1 <- martTable1[order(martTable1$ensembl_gene_id),]
#Check for duplicated Ensembl gene IDs
sum(duplicated(martTable1$ensembl_gene_id))

#Verify that the DGE gene ID annotation information corresponds with the zebrafish BioMart annotation table before creating a master annotation table
sum(DGE$genes$genes != martTable1$ensembl_gene_id)

#Create master BioMart annotation table
martTable <- data.frame(genes = DGE$genes$genes)
#Add gene symbol information to the master BioMart annotation table
martTable$Symbol <- martTable1$external_gene_name

#Ensure that the Symbol column does not have NAs or blanks
sum(is.na(martTable$Symbol))
sum(martTable$Symbol == '')
if(sum(martTable$Symbol == '') > 0){
    for (i in seq_along((martTable$Symbol))){
        if (martTable$Symbol[i] == ''){
            martTable$Symbol[i] <- martTable$genes[i]
        }
    }
}

#Ensures that each value in the Symbol column is unique
sum(duplicated(martTable$Symbol))
#Add each Symbol to gene_symbols
gene_symbols <- martTable$Symbol[1] 
for (i in 2:length(martTable$Symbol)){
    print(i)
    if (martTable$Symbol[i] %in% gene_symbols){
        #If the Symbol is not unique (i.e. already included in gene_symbols), then add a blank space to the end of its name, so that it appears the same but registers as unique for row names to be included in heatmaps
        while(martTable$Symbol[i] %in% gene_symbols){
            martTable$Symbol[i] <- paste0(martTable$Symbol[i], ' ')
        }
    }
    #Add each symbol to the cumulative gene_symgols list
    gene_symbols <- c(gene_symbols, martTable$Symbol[i])
}
##Verify that each value in the Symbol column is unique
sum(duplicated(martTable$Symbol))

#Add gene type to the master BioMart annotation table
martTable$Biotype <- martTable1$gene_biotype


##----Add human gene annotation information----
martTable2 <- getBM(attributes = human_attributes, filters = filts, values = vals, mart = ensemblDrer89)
#Order results by zebrafish Ensembl gene ID
martTable2 <- martTable2[order(martTable2$ensembl_gene_id),]

#Check for duplicated zebrafish Ensembl ID values
sum(duplicated(martTable2$ensembl_gene_id))
#Create column in zebrafish to human BioMart annotation table for combining multiple human orthologs (gene ID) into one line
martTable2$hsapiens_homolog_ensembl_gene_all <- martTable2$hsapiens_homolog_ensembl_gene
#Create column in zebrafish to human BioMart annotation table for combining multiple human orthologs (gene name) into one line
martTable2$hsapiens_homolog_associated_gene_name_all <- martTable2$hsapiens_homolog_associated_gene_name


##--Combine multiple human gene ID and symbol information into one line
for (i in 2:length(DGE$genes$genes)){
    print(i)
    if (martTable2$ensembl_gene_id[(i-1)] == martTable2$ensembl_gene_id[i]){
        while(martTable2$ensembl_gene_id[(i-1)] == martTable2$ensembl_gene_id[i]){
            martTable2$hsapiens_homolog_ensembl_gene[i] <- martTable2$hsapiens_homolog_ensembl_gene[i-1]
            martTable2$hsapiens_homolog_ensembl_gene_all[i] <- paste(martTable2$hsapiens_homolog_ensembl_gene_all[i-1], martTable2$hsapiens_homolog_ensembl_gene_all[i], sep = ', ')
            
            martTable2$hsapiens_homolog_associated_gene_name[i] <- martTable2$hsapiens_homolog_associated_gene_name[i-1]
            martTable2$hsapiens_homolog_associated_gene_name_all[i] <- paste(martTable2$hsapiens_homolog_associated_gene_name_all[i-1], martTable2$hsapiens_homolog_associated_gene_name_all[i], sep = ', ')
            martTable2 <- martTable2[-(i-1),]
        }
    }
}

#Check for duplicated zebrafish Ensembl ID values
sum(duplicated(martTable2$ensembl_gene_id))
#Verify that the master BioMart annotation table information corresponds with the zebrafish to human ortholog BioMart annotation table before adding its columns to the master table
sum(martTable$genes != martTable2$ensembl_gene_id)

#Add human Ensembl gene ID ortholog information to the master BioMart annotation table
martTable$Hsap_ID <- martTable2$hsapiens_homolog_ensembl_gene
#Add human Ensembl gene ID multiple ortholog information to the master BioMart annotation table
martTable$Hsap_IDs <- martTable2$hsapiens_homolog_ensembl_gene_all
#Add human gene symbol ortholog information to the master BioMart annotation table
martTable$Hsap_symbol <- martTable2$hsapiens_homolog_associated_gene_name
#Add human gene symbol multiple ortholog information to the master BioMart annotation table
martTable$Hsap_symbols <- martTable2$hsapiens_homolog_associated_gene_name_all

##----Add mouse gene annotation information----
martTable3 <- getBM(attributes = mouse_attributes, filters = filts, values = vals, mart = ensemblDrer89)
#Order results by zebrafish Ensembl gene ID
martTable3 <- martTable3[order(martTable3$ensembl_gene_id),]

#Check for duplicated zebrafish Ensembl ID values
sum(duplicated(martTable3$ensembl_gene_id))
#Create column in zebrafish to mouse BioMart annotation table for combining multiple mouse orthologs (gene ID) into one line
martTable3$mmusculus_homolog_ensembl_gene_all <- martTable3$mmusculus_homolog_ensembl_gene
#Create column in zebrafish to mouse BioMart annotation table for combining multiple mouse orthologs (gene name) into one line
martTable3$mmusculus_homolog_associated_gene_name_all <- martTable3$mmusculus_homolog_associated_gene_name

##--Combine multiple mouse gene ID and symbol information into one line
for (i in 2:length(DGE$genes$genes)){
    print(i)
    if (martTable3$ensembl_gene_id[(i-1)] == martTable3$ensembl_gene_id[i]){
        while(martTable3$ensembl_gene_id[(i-1)] == martTable3$ensembl_gene_id[i]){
            martTable3$mmusculus_homolog_ensembl_gene[i] <- martTable3$mmusculus_homolog_ensembl_gene[i-1]
            martTable3$mmusculus_homolog_ensembl_gene_all[i] <- paste(martTable3$mmusculus_homolog_ensembl_gene_all[i-1], martTable3$mmusculus_homolog_ensembl_gene_all[i], sep = ', ')
            martTable3$mmusculus_homolog_associated_gene_name[i] <- martTable3$mmusculus_homolog_associated_gene_name[i-1]
            martTable3$mmusculus_homolog_associated_gene_name_all[i] <- paste(martTable3$mmusculus_homolog_associated_gene_name_all[i-1], martTable3$mmusculus_homolog_associated_gene_name_all[i], sep = ', ')
            martTable3 <- martTable3[-(i-1),]
        }
    }
}

#Check for duplicated zebrafish Ensembl ID values
sum(duplicated(martTable3$ensembl_gene_id))
#Verify that the master BioMart annotation table information corresponds with the zebrafish to mouse ortholog BioMart annotation table before adding its columns to the master table
sum(martTable$genes != martTable3$ensembl_gene_id)

#Add mouse Ensembl gene ID ortholog information to the master BioMart annotation table
martTable$Mmus_ID <- martTable3$mmusculus_homolog_ensembl_gene
#Add mouse Ensembl gene ID multiple ortholog information to the master BioMart annotation table
martTable$Mmus_IDs <- martTable3$mmusculus_homolog_ensembl_gene_all
#Add mouse gene symbol ortholog information to the master BioMart annotation table
martTable$Mmus_symbol <- martTable3$mmusculus_homolog_associated_gene_name
#Add mouse gene symbol multiple ortholog information to the master BioMart annotation table
martTable$Mmus_symbols <- martTable3$mmusculus_homolog_associated_gene_name_all

#Ensure that the BioMart annotation table does not contain NA values
martTable[is.na(martTable)] <- ''

##----OR: Add annotation information from file----
#martTable <- read.delim(file = file.choose(), header = TRUE, stringsAsFactors = FALSE)

#Verify that the master BioMart annotation table corresponds with the order of genes in the DGE object
sum(DGE$genes$genes != martTable$genes)
#Add BioMart annotation information to the DGE object
DGE$genes$Symbol <- martTable$Symbol
DGE$genes$Biotype <- martTable$Biotype
DGE$genes$Hsap_ID <- martTable$Hsap_ID
DGE$genes$Hsap_IDs <- martTable$Hsap_IDs
DGE$genes$Hsap_symbol <- martTable$Hsap_symbol
DGE$genes$Hsap_symbols <- martTable$Hsap_symbols
DGE$genes$Mmus_ID <- martTable$Mmus_ID
DGE$genes$Mmus_IDs <- martTable$Mmus_IDs
DGE$genes$Mmus_symbol <- martTable$Mmus_symbol
DGE$genes$Mmus_symbols <- martTable$Mmus_symbols

#Write DGE annotation table to file
write.table(DGE$genes, file = paste0(output_filepath,'/Zfish_martTable_Zv10_89.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----7) Additional filtering based on edgeR vignette and (Chen 2016)----
#Determine counts per million corresponding to a minimum of 10-20 read counts 
mean_cpm_gene_count_of_15 <- as.numeric(cpm(15, mean(DGE$samples$lib.size)))
#Determine the number of samples in the smallest group
num_samples_smallest_group <- min(table(sampleInfo$Group))  #4 is the minimum of samples in each group

#Only keep genes that have a minimum of 10-20 counts in at least 4 samples to be included in analysis
keep.exprs <- rowSums(cpm(DGE)> mean_cpm_gene_count_of_15) >= num_samples_smallest_group
DGE <- DGE[keep.exprs, , keep.lib.sizes=FALSE]

#Write background genes to file
write.table(DGE$genes$genes, file = paste0(output_filepath,'/Zfish_background_genes.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#Write human ortholog background genes to file
write.table(DGE$genes$Hsap_ID[DGE$genes$Hsap_ID != ''], file = paste0(output_filepath,'/Metacore_background_genes.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#Write annotation table to file
write.table(DGE$genes, file = paste0(output_filepath,'/DGE_annotation.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----8) Normalize libraries----
DGE <- calcNormFactors(DGE)
DGE$samples
write.table(DGE$samples, file = paste0(output_filepath,'/DGE_libsize_normfactors.txt'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

##----Create normalized count data for use with heatmaps and PCA plots----
#Calculate cpm to factor in TMM normalization
rltdata <- cpm(y = DGE, normalized.lib.sizes = TRUE) #get normalized count data in counts per million units
#Multiply normalized cpm numbers by library size per million to get count values normalized by library size
for (bio_sample in colnames(rltdata)){
    rltdata[,bio_sample] <- (rltdata[,bio_sample])*((DGE$samples[bio_sample,'lib.size']/1000000))
}
#Assign TMM-normalized counts to the 'normalized_counts' slot of the DGE object
DGE$normalized_counts <- rltdata
#Regular log normalization of count data (corrects for heteroskedasticy)
rltdata <- rlogTransformation(round(rltdata))
#Add gene labels to TMM-normalized, rlog transformed data
row.names(rltdata) <- row.names(DGE$counts)
#Assign TMM-normalized, rlog transformed counts to the 'rlog_normalized' slot of the DGE object
DGE$rlog_normalized <- rltdata
#Write TMM-normalized, rlog transformed count table to file
write.table(DGE$rlog_normalized, file = paste0(output_filepath,'/DGE_rlog_normalized_counts.txt'), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


##----9) Make the design matrix----
design <- model.matrix(~0 + DGE$samples$group) #the zero means that there is no intercept in the model, this helps with contrasts in later steps
#Assign sample names to the rows of the design matrix
rownames(design) <- colnames(DGE)
#Assign group names to the columns of the design matrix
colnames(design) <- levels(DGE$samples$group)
design


##----10) Estimate dispersions----
DGE <- estimateDisp(DGE, design, robust = TRUE)
DGE$common.dispersion
sqrt(DGE$common.dispersion)
plotBCV(DGE, main = 'Dispersion estimation')


##----11) Fit data to the design model----
fit <- glmQLFit(DGE, design, robust = TRUE)
plotQLDisp(glmfit = fit, main = 'Quasi-liklihood Dispersion')


##----12) Design contrast matrix to extract comparisons of interest----
my_contrasts <- makeContrasts(`CM_TCDD-CM_DMSO` = CM_TCDD - CM_DMSO,
                              `SM_TCDD-SM_DMSO` = SM_TCDD - SM_DMSO, 
                              `SM_DMSO-CM_DMSO` = SM_DMSO - CM_DMSO,
                              `SM_TCDD-CM_TCDD` = SM_TCDD - CM_TCDD, 
                              levels = design)


##----13) Set up ComplexHeatmap parameters----
#Set heatmap global text options
ht_global_opt(heatmap_row_names_gp = gpar(fontfamily = 'sans'),
              heatmap_row_title_gp = gpar(fontfamily = 'sans'),
              heatmap_legend_title_gp = gpar(fontfamily = 'sans', fontsize = 8),
              heatmap_legend_labels_gp = gpar(fontfamily = 'sans', fontsize = 8),
              heatmap_column_names_gp = gpar(fontfamily = 'sans'),
              heatmap_column_title_gp = gpar(fontfamily = 'sans'),
              annotation_legend_title_gp = gpar(fontfamily = 'sans', fontsize = 8),
              annotation_legend_labels_gp = gpar(fontfamily = 'sans', fontsize = 8))

#Set up color schemes for row annotation information
biotypes <- data.frame(biotypes = c('antisense', 'lincRNA', 'misc_RNA', 'Mt_rRNA', 'Mt_tRNA', 'processed_transcript', 'protein_coding', 'rRNA','sense_overlapping'))
Biotype_colScale <- data.frame(Biotype = biotypes$biotypes,
                               Color = I(colorschemes$Categorical.12[c(12,5,3,7,8,9,10,11,6)]))
row.names(Biotype_colScale) <- Biotype_colScale$Biotype
show_col(Biotype_colScale$Color)

#Set up column annotation information
col_anno_info <- data.frame(Treatment = sampleInfo$Treatment, Morphant = sampleInfo$Morphant)

col_anno <- HeatmapAnnotation(
    df = col_anno_info,
    col = list(
        Treatment = c('DMSO' = I(colorschemes$Categorical.12[7]),
                      'TCDD' = I(colorschemes$Categorical.12[12])),
        Morphant = c('CM' = I(colorschemes$Categorical.12[4]),
                     'SM'= I(colorschemes$Categorical.12[2]))),
    show_legend = TRUE,
    annotation_legend_param = list(
        Treatment = list(nrow = 1, 
                         title = "Treatment", 
                         title_position = "topcenter",
                         legend_direction = 'horizontal', 
                         fontfamily = 'sans'),
        Morphant = list(nrow = 1, 
                        title = "Morphant", 
                        title_position = "topcenter",
                        legend_direction = 'horizontal', 
                        fontfamily = 'sans')),
    show_annotation_name = TRUE,
    annotation_name_side = 'left',
    annotation_name_gp = gpar(fontsize = 8),
    height = unit(5, "mm")
)#END col_anno

#Creates 'padding' object so column annotation names are not cut off
padding <- unit.c(unit(1, "mm"), 
                  grobWidth(textGrob("  Treatment  ")) - unit(22.5, "mm"),
                  unit(1, "mm"),
                  unit(1, "mm"))


##----14) Find significant genes for each comparison without a fold change cutoff----
fold_change = 1 

for(comparison in colnames(my_contrasts)){
    Sys.sleep(.01)
    print(comparison)
    #Make a label for each contrast
    cutoff_contrast <- paste0(comparison, '_log2(',fold_change,')')
    
    #Test for differential expression using quasi-likelihood F-tests
    comparisons <- glmQLFTest(glmfit = fit, 
                              contrast = my_contrasts[,comparison])
    #Extract the significant differentially expressed genes that pass the Benjamini-Hochberg FDR cutoff of 0.05
    sig_table <- data.frame(topTags(comparisons,
                                    n = Inf, adjust.method = "BH",
                                    p.value = 0.05)) 
    
    ##----create output folder for each contrast----
    if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast))){dir.create(paste0(output_filepath,'/', cutoff_contrast))}
    
    if(ncol(sig_table) > 0){
        #Create data objects for each DGELRT object created in the loop
        assign(cutoff_contrast, comparisons)
        #Create data table for each table of significant differentially expressed genes
        assign(paste(cutoff_contrast, 'sigGenes', sep = '_'), sig_table)
        #Write table of significant differentially expressed genes to text file
        write.table(x = sig_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sig_log_ratios.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write table of significant differentially expressed genes to Excel file
        xlsx::write.xlsx2(x = sig_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sig_log_ratios.xlsx'), row.names = FALSE)
        .jcall("java/lang/System", method = "gc") #helps with java memory problems
        
        #Create table for use with MetaCore
        metacore_table <- sig_table[sig_table$Hsap_ID != '', c('Hsap_ID','logFC','FDR')]
        #Write significanly differentially expressed upregulated human orthologs to text file
        write.table(x = metacore_table[(metacore_table$logFC > 0),],file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'UR_Metacore_sigGenes.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write significanly differentially expressed downregulated human orthologs to text file
        write.table(x = metacore_table[(metacore_table$logFC < 0),],file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'DR_Metacore_sigGenes.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write significanly differentially expressed human orthologs to text file
        write.table(x = metacore_table,file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_Metacore_sigGenes.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed genes to text file
        write.table(x = sig_table$genes, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sigGene_list.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed upregulated genes to text file
        write.table(x = sig_table$genes[sig_table$logFC > 0], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'UR_sigGene_list.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed downregulated genes to text file
        write.table(x = sig_table$genes[sig_table$logFC < 0], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'DR_sigGene_list.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        ##----Heatmap of top genes sorted by fold change----
        #Order genes based on their absolute log2(fold change)
        heatmap_selector_FC <- sig_table[order(abs(sig_table$logFC), decreasing = TRUE),] 
        #Select the top 30 genes
        heatmap_selector_FC <- heatmap_selector_FC[1:30,'genes']
        #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 genes
        heatmap_selector_FC <- row.names(rltdata) %in% heatmap_selector_FC
        #Select data from the TMM-normalized, rlog transformed count table matching the top 30 genes
        rlogdataFC <- data.frame(rltdata[heatmap_selector_FC,])
        #Select data from the annotation table matching the top 30 genes
        phenodataFC <- DGE$genes[heatmap_selector_FC, 'Biotype',drop = FALSE]
        #Label count data with gene symbols
        row.names(rlogdataFC) <- DGE$genes$Symbol[heatmap_selector_FC]
        #Label phenotype data with gene symbols
        row.names(phenodataFC) <- DGE$genes$Symbol[heatmap_selector_FC]
        #Scale the transformed count data
        rlogdataFC <- t(scale(t(rlogdataFC)))
        
        ##Set up row annotation information
        row_anno_info <- phenodataFC
        row_key <- names(table(as.character(row_anno_info$Biotype)))
        row_cols <- character()
        for (i in seq.int(1:length(row_key))){
            row_col <- Biotype_colScale[row_key[i],'Color']
            row_cols <- c(row_cols,row_col)
        }
        row_col_list <- list(Biotype = row_cols)
        names(row_col_list[[1]]) <- row_key
        
        ##Create main heatmap layer
        hmap <- Heatmap(matrix = rlogdataFC, 
                        name = 'main_mat',
                        heatmap_legend_param = 
                            list(title = 'Z-score',
                                 color_bar = 'continuous', 
                                 legend_direction = 'horizontal', 
                                 title_position = 'topcenter', 
                                 legend_side = 'bottom', 
                                 fontfamily = 'sans'),
                        column_title = cutoff_contrast,
                        column_title_side = 'top',
                        column_title_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(fontsize = 8),
                        column_dend_height = unit(8,'mm'),
                        top_annotation = col_anno,
                        show_row_names = FALSE,
                        show_row_dend = TRUE,
                        row_dend_width = unit(13, "mm"),
                        show_heatmap_legend = TRUE
        )#END hmap
        
        ##Create row annotations layer
        rmap <- Heatmap(matrix = row_anno_info, 
                        name = 'Biotype', 
                        col = row_col_list$Biotype,
                        column_title_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = 8),
                        show_row_names = TRUE,
                        row_title_gp = gpar(fontsize = 8),
                        row_names_side = 'right',
                        row_names_gp = gpar(fontsize = 8),
                        heatmap_legend_param = (
                            Biotype = list(nrow = if(length(row_key) > 5){2}else{1}, 
                                           title = "Biotype", 
                                           title_position = "topcenter",
                                           legend_direction = 'horizontal', 
                                           fontfamily = 'sans')),
                        width = unit(2.5,'mm')
        )#END rmap
        
        #Draw the combined heatmap and write to file
        pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FC.pdf'), width = 3.4, height = 6)            
        draw(hmap + rmap, padding = padding, gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
        dev.off()
    }
}    


##----15) Write session information----
save.image(paste0(output_filepath,'/RNAseq_RData_',Sys.Date(),'.RData'))
writeLines(capture.output({
    cat('###########################################################################\n##  Script and Data Analysis by Cheryl Dunham                            ##\n##                                                                       ##\n##  Budding Bioinformatician                                             ##\n##  Tanguay Laboratory | Sinnhuber Aquatic Research Laboratory (SARL)    ##\n##  28645 E Hwy 34                                                       ##\n##  Corvallis, OR 97333                                                  ##\n##  Oregon State University                                              ##\n##                                                                       ##\n##  (541) 737-6500                                                       ##\n##  dunhamcg@gmail.com                                                   ##\n###########################################################################\n\n')
    cat('Input filepath:', input_filepath, '\n\n')
    cat('Output filepath:', output_filepath, '\n\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\n')
    print(sessionInfo())
    cat('\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nRStudio Version Information\n\n')
    print(RStudio.Version())
    cat('_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nAttached Package Citation Information\n\n')
    print(lapply(all_packages, citation))
    
}), paste0(output_filepath,'/RNAseq_sessionInfo_',Sys.Date(),'.txt'))
