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


##Use to analyse Zebrafish RNA-seq data


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
output_filepath <- 'C:/Users/dunhamc/Google Drive/RNA_Seq/FRCs/Output_20190703'
#File path of data to be analyzed, this also helps with documentation later on
input_filepath <- 'C:/Users/dunhamc/Google Drive/RNA_Seq/FRCs/HTSeq_2.1.0'
#Label to append to file names
identifier <- '_FRC_RNA'
#Label to add to images
appender <- ' (FRC RNA-seq)'
#Server file path indicating the location of trimmed reads
skewer_filepath <- '/nfs1/EMT/Tanguay_Lab/Bioinfo/dunhamc/FRCs/RNAseq/Skewer'
#Server file path indicating the location of aligned reads
hisat2_filepath <- '/nfs1/EMT/Tanguay_Lab/Bioinfo/dunhamc/FRCs/RNAseq/HISAT2_2.1.0'
#Server file path indicating the location of gene quantification files
htseq_filepath <- '/nfs1/EMT/Tanguay_Lab/Bioinfo/dunhamc/FRCs/RNAseq/HTSeq_2.1.0'


##----3) Import and combine sequence data----
filenames <- list.files(path = input_filepath, pattern = '.tab')
#Set the working directory
setwd(input_filepath)

HTseq_combo <- function(){
    #Combines data from multiple HTseq-count output files (sans the --additional-attr option) into a single dataframe
    #Be sure to be in the correct working directory before running the function
    
    #List all file names from the current working folder that end in .tab file extension
    filenames <- list.files(pattern = '.tab')
    #Read in the count table from the first file in the folder
    initial_count <- read.table(filenames[1], sep ='\t', stringsAsFactors = FALSE, header = FALSE)
    #Create the dataframe to hold information from all .tab files in the current working folder
    counts_combo <- as.data.frame(matrix(ncol = length(filenames) + 2, nrow = length(initial_count[,1])+1))
    #Set the first column of the combined dataframe to the Ensembl gene ID (the first column of HTSeq-count files)
    counts_combo[,1] <- c(initial_count[,1],'')
    #Set the second column of the combined dataframe to the gene symbol additional attribute
    counts_combo[,2] <- c(initial_count[,2],'')
    #Set the third column of the combined dataframe to the corresponding gene count values for the first .tab file in the current folder
    counts_combo[,3] <- as.numeric(c(initial_count[,3],''))
    #Create a statistic for total annotated reads
    counts_combo[nrow(counts_combo),1] <- 'Total_annotated_reads'
    #Calculate the total number of annotated reads detected by HTseq-count
    counts_combo[nrow(counts_combo),3] <- sum(initial_count[1:(nrow(initial_count)-5),3])
    #Set the name of the first two columns of the combined dataframe to the Ensembl_ID and gene_name (the first two columns of HTSeq-count files) and the file name of the first .tab file in the current folder 
    names(counts_combo)[1:3] <- c('Ensembl_ID','gene_name',filenames[1])
    
    #Create loop to add gene count information to the combined dataframe from the rest of the .tab files in the current folder
    for (i in 2:length(filenames)){
        #Read in HTSeq-count data table
        counts <- read.table(filenames[i], sep = '\t', stringsAsFactors = FALSE, header = FALSE)
        #Add gene count information to the combined dataframe
        counts_combo[,i+2] <- as.numeric(c(counts[,3],''))
        #Add total annotated reads statistic to the combined dataframe
        counts_combo[nrow(counts_combo),i+2] <- sum(counts[1:(nrow(counts)-5),3])
        #Link the HTSeq-count file name to its corresponding gene counts 
        names(counts_combo)[i+2] <- filenames[i]
    }
    
    #Create count statistics table
    counts_stats <- counts_combo[(length(counts_combo[,1])-5):(length(counts_combo[,1])),]
    #Save complete HTseq-count table 
    counts_table <- counts_combo
    #Remove non-count data information from the combined data table (last 6 rows of the dataframe)
    counts_combo <- counts_combo[1:(length(counts_combo[,1])-6),]
    #Label rows of the combined data table with the corresponding Ensembl gene ID
    row.names(counts_combo) <- counts_combo[,1]
    #Combine the three data tables into a list
    count_list <- list(counts_combo,counts_stats,counts_table)
    #Return the list of three data tables
    return(count_list)
}

#Merge count data and retrieve count data statistics
master_counts <- HTseq_combo()
#Retrieve merged count data table
merged_counts <- as.data.frame(master_counts[[1]], stringsAsFactors = FALSE)
#Retrieve merged count data statistics
count_stats <- as.data.frame(master_counts[[2]], stringsAsFactors = FALSE)
#Retrieve merged count data table and statistics
combined_count_table <- as.data.frame(master_counts[[3]], stringsAsFactors = FALSE)

#Set the working directory to your output file path
setwd(output_filepath)


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

#Ensure the 'Filename' column is of the character class and not a factor
sampleInfo$Filename <- as.character(sampleInfo$Filename)

#Verify that the file names associated with the combined count table matches with the file names provided in the sampleInfo spreadsheet
colnames(combined_count_table) %in% sampleInfo$Filename

#Reorder the merged count and statistics table by the desired filename order from your sample information sheet
combined_count_table <- combined_count_table[,c('Ensembl_ID','gene_name',sampleInfo$Filename)]
#Write merged count and statistics table to file
write.table(combined_count_table, file = paste0(output_filepath,'/Combined_HTSeq-count_table',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Reorder the count statistics table by the desired filename order from your sample information sheet
count_stats <- count_stats[,c('Ensembl_ID','gene_name',sampleInfo$Filename)]
#Verify the column info from count data matches row info from sampleInfo sheet
colnames(count_stats[,-c(1:2)]) == sampleInfo$Filename
sum(colnames(count_stats[,-c(1:2)]) != sampleInfo$Filename)
#Change the column names of the count stat table to the abbreviated ID provided in the sample info sheat
colnames(count_stats) <- c('Ensembl_ID','gene_name',sampleInfo$subID)
#Write count statistics table to file
write.table(count_stats, file = paste0(output_filepath,'/Count_stats',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Verify that the file names associated with the merged count table matches with the filenames provided in the sampleInfo spreadsheet
colnames(merged_counts[-c(1,2)]) %in% sampleInfo$Filename
sampleInfo$Filename %in% colnames(merged_counts)

#Create gene count matrix arranged by experimental sample order from sample information spreadsheet
experiment_counts <- merged_counts[,sampleInfo$Filename]
#Verify the column info from count data matches row info from sampleInfo sheet
colnames(experiment_counts) == sampleInfo$Filename
sum(colnames(experiment_counts) != sampleInfo$Filename)

#Rename experimental count table with corresponding sample ID abbreviation
colnames(experiment_counts) <- sampleInfo$subID
#Write sample information sheet to file for documentation
write.xlsx(x = sampleInfo, file = paste0(output_filepath,'/sampleInfo',identifier,'.xlsx'), row.names = FALSE)
#Write count matrix to file
write.table(cbind(data.frame('Ensembl_ID' = row.names(experiment_counts), experiment_counts)), file = paste0(output_filepath,'/Count_matrix',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----5) Make the edgeR digital gene expression (DGE) data object----
#See the Bioconductor workflow, 'From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline' <https://www.bioconductor.org/packages/release/workflows/vignettes/RnaSeqGeneEdgeRQL/inst/doc/edgeRQL.html>, for an up-to-date RNAseq differential expression workflow
DGE <- DGEList(counts = experiment_counts, 
               genes = row.names(experiment_counts), 
               group = sampleInfo$Group) 
#Add additional sample information to the DGE object
DGE$samples$Chemical_ID <- sampleInfo$Chemical
DGE$samples$Chemical_Name <- sampleInfo$Chemical_Name
DGE$samples$subID <- sampleInfo$subID
DGE$samples$Lane <- sampleInfo$Lane


##----6) Add gene annotation information----
#Zv10 (GRCz10) Ensembl release 91 annotation information
ensemblDrer91 <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'drerio_gene_ensembl', host = 'dec2017.archive.ensembl.org') 

##--Biomart options--
#Set the data columns to be added to the new dataframe (attributes) and filter the data based on user specifications and associated values (from the main Ensembl Biomart site: dec2017.archive.ensembl.org/biomart/martview/)

#Filters
#listFilters(ensemblDrer91)  #to see all possible filter names
#Filters for use in all annotation retrieval
filts <- 'ensembl_gene_id' 
vals <- row.names(DGE) #values corresponding to filters (filts)

#Attributes  
#listAttributes(ensemblDrer91)  #to see all possible attribute names
#Attributes for zebrafish annotation information
zebrafish_ID_attributes <- c('ensembl_gene_id', 'external_gene_name')
zebrafish_phenotype_attributes <- c('ensembl_gene_id', 'description', 'gene_biotype', 'phenotype_description')
#Attributes for zebrafish to human ortholog information
human_attributes <- c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_associated_gene_name')
#Attributes for zebrafish Gene Ontology (GO) information
GO_attributes <-  c('ensembl_gene_id','go_id','name_1006','definition_1006','namespace_1003') 

##----Add zebrafish gene annotation information----
#Retrieve zebrafish Ensembl gene ID, gene name/symbol, and gene type for those genes included in the DGE object
martTable1a <- getBM(attributes = zebrafish_ID_attributes, filters = filts, values = vals, mart = ensemblDrer91)

#Ensure the table is ordered by Ensembl gene ID
martTable1a <- martTable1a[order(martTable1a$ensembl_gene_id),]
#Check for duplicated Ensembl gene IDs
sum(duplicated(martTable1a$ensembl_gene_id))

#Verify that the DGE gene ID annotation information corresponds with the zebrafish BioMart annotation table before creating a master annotation table
sum(DGE$genes$genes != martTable1a$ensembl_gene_id)

#Create master BioMart annotation table
martTable <- data.frame(Genes = DGE$genes$genes)
#Reassign the Genes column as a character class
martTable$Genes <- as.character(martTable$Genes)
#Add gene symbol information to the master BioMart annotation table
martTable$Symbol <- martTable1a$external_gene_name

#Ensure that the Symbol column does not have NAs or blanks
sum(is.na(martTable$Symbol))
sum(martTable$Symbol == '')
#Assign the Ensembl gene ID as the symbol for those genes lacking corresponding symbols in BioMart
if(sum(martTable$Symbol == '') > 0){
    for (i in seq_along((martTable$Symbol))){
        if (martTable$Symbol[i] == ''){
            martTable$Symbol[i] <- martTable$Genes[i]
        }
    }
}

#Ensure that each value in the Symbol column is unique
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
#Verify that each value in the Symbol column is unique
sum(duplicated(martTable$Symbol))

#Add zebrafish phenotype information
martTable1b <- getBM(attributes = zebrafish_phenotype_attributes, filters = filts, values = vals, mart = ensemblDrer91)
#Ensure the table is ordered by Ensembl gene ID
martTable1b <- martTable1b[order(martTable1b$ensembl_gene_id),]
#Check for duplicated Ensembl gene IDs
sum(duplicated(martTable1b$ensembl_gene_id))

#Combine phenotpye information into one line
for (i in 2:length(DGE$genes$genes)){
    print(i)
    if (martTable1b$ensembl_gene_id[(i-1)] == martTable1b$ensembl_gene_id[i]){
        while(martTable1b$ensembl_gene_id[(i-1)] == martTable1b$ensembl_gene_id[i]){
            martTable1b$phenotype_description[i] <- paste(martTable1b$phenotype_description[i], martTable1b$phenotype_description[i-1], sep = ', ')
            martTable1b <- martTable1b[-(i-1),]
        }
    }
}
#Check for duplicated Ensembl gene IDs
sum(duplicated(martTable1b$ensembl_gene_id))

#Add annotation information to the main BioMart table
martTable$Description <- martTable1b$description
martTable$Biotype <- martTable1b$gene_biotype
martTable$Phenotype <- martTable1b$phenotype_description

##----Add human gene annotation information----
martTable2 <- getBM(attributes = human_attributes, filters = filts, values = vals, mart = ensemblDrer91)
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
            martTable2$hsapiens_homolog_ensembl_gene_all[i] <- paste(martTable2$hsapiens_homolog_ensembl_gene_all[i], martTable2$hsapiens_homolog_ensembl_gene_all[i-1], sep = ', ')
            martTable2$hsapiens_homolog_associated_gene_name_all[i] <- paste(martTable2$hsapiens_homolog_associated_gene_name_all[i], martTable2$hsapiens_homolog_associated_gene_name_all[i-1], sep = ', ')
            martTable2 <- martTable2[-(i-1),]
        }
    }
}

#Check for duplicated zebrafish Ensembl ID values
sum(duplicated(martTable2$ensembl_gene_id))
#Verify that the gene order of the master BioMart annotation table  corresponds with the zebrafish to human ortholog BioMart annotation table before adding its columns to the master table
sum(martTable$Genes != martTable2$ensembl_gene_id)

#Add human Ensembl gene ID ortholog information to the master BioMart annotation table
martTable$Hsap_ID <- martTable2$hsapiens_homolog_ensembl_gene
#Add human Ensembl gene ID multiple ortholog information to the master BioMart annotation table
martTable$Hsap_IDs <- martTable2$hsapiens_homolog_ensembl_gene_all
#Add human gene symbol ortholog information to the master BioMart annotation table
martTable$Hsap_symbol <- martTable2$hsapiens_homolog_associated_gene_name
#Add human gene symbol multiple ortholog information to the master BioMart annotation table
martTable$Hsap_symbols <- martTable2$hsapiens_homolog_associated_gene_name_all

##----Add GO annotation information----
GOmartTable <- getBM(attributes = GO_attributes, filters = filts, values = vals, mart = ensemblDrer91)
names(GOmartTable) <- c('ensembl_gene_id','go_id','go_name','go_definition','go_domain')
#Remove NA values from GOmartTable
GOmartTable[is.na(GOmartTable)] <- ''
#Write Gene Ontology table to file
write.table(GOmartTable, file = paste0(output_filepath,'/Zfish_GOTable_GRCz10_91.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Order GOmartTable by GO_ID, removing missing values
GOmartTable <- GOmartTable[GOmartTable$go_id != '',]
#Reorder GOmartTable by Ensembl Gene ID
GOmartTable <- GOmartTable[order(GOmartTable$ensembl_gene_id),]
sum(duplicated(GOmartTable$ensembl_gene_id))

##--Combine GO ID infoormation into one line
for (i in 2:(length(unique(GOmartTable$ensembl_gene_id)))){
    print(i)
    if (GOmartTable$ensembl_gene_id[(i-1)] == GOmartTable$ensembl_gene_id[i]){
        while(GOmartTable$ensembl_gene_id[(i-1)] == GOmartTable$ensembl_gene_id[i]){
            GOmartTable$go_id[i] <- paste0(GOmartTable$go_id[i-1], '; ', GOmartTable$go_id[i])
            GOmartTable$go_name[i] <- paste0(GOmartTable$go_name[i-1], '; ', GOmartTable$go_name[i])
            GOmartTable$go_definition[i] <- paste0(GOmartTable$go_definition[i-1], '; ', GOmartTable$go_definition[i])
            GOmartTable$go_domain[i] <- paste0(GOmartTable$go_domain[i-1], '; ', GOmartTable$go_domain[i])
            GOmartTable <- GOmartTable[-(i-1),]
        }
    }
}
sum(duplicated(GOmartTable$ensembl_gene_id))
#Label rows with Ensembl gene ID
row.names(GOmartTable) <- GOmartTable$ensembl_gene_id
#Write organized Gene Ontology table to file
write.table(GOmartTable, file = paste0(output_filepath,'/Zfish_GOTableCombo_GRCz10_91.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#Add GO information to the master BioMart table
martTable$GO_IDs <- NA
martTable$GO_Terms <- NA
martTable$GO_Definitions <- NA
martTable$GO_Domains <- NA
row.names(martTable) <- martTable$Genes
n = 1
for(id in GOmartTable$ensembl_gene_id){
    print(n)    
    martTable[id, "GO_IDs"] <- GOmartTable[id, "go_id"]
    martTable[id, "GO_Terms"] <- GOmartTable[id, "go_name"]
    martTable[id, "GO_Definitions"] <- GOmartTable[id, "go_definition"]
    martTable[id, "GO_Domains"] <- GOmartTable[id, "go_domain"]
    n = n+1
}


##----OR: Add annotation information from file----
#martTable <- read.delim(file = file.choose(), header = TRUE, stringsAsFactors = FALSE)
martTable <- read.delim(file = 'C:/Users/dunhamc/Google Drive/RNA_Seq/Zfish_martTable_Zv10_91.txt', header = TRUE, stringsAsFactors = FALSE)

#Ensure that the BioMart annotation table does not contain NA values
martTable[is.na(martTable)] <- ''
unique(martTable$Biotype)
#Replace biotypes with '_' in the names
cat(paste(shQuote(unique(martTable$Biotype)), collapse=", "),'\n')
martTable <- mutate_all(martTable, 
                        ~case_when(. == 'protein_coding' ~ 'protein coding',
                                   . == 'polymorphic_pseudogene' ~ 'polymorphic pseudogene',
                                   . == 'processed_transcript' ~ 'processed transcript',
                                   . == 'unprocessed_pseudogene' ~ 'unprocessed pseudogene',
                                   . == 'antisense_RNA' ~ 'antisense RNA',
                                   . == 'processed_pseudogene' ~ 'processed pseudogene',
                                   . == 'transcribed_unprocessed_pseudogene' ~ 'transcribed unprocessed pseudogene',
                                   . == 'sense_intronic' ~ 'sense intronic',
                                   . == 'IG_V_pseudogene' ~ 'IG V pseudogene',
                                   . == 'TR_J_gene' ~ 'TR J gene',
                                   . == 'IG_J_pseudogene' ~ 'IG J pseudogene',
                                   . == 'TR_V_gene' ~ 'TR V gene',
                                   . == 'Mt_tRNA' ~ 'Mt tRNA',
                                   . == 'Mt_rRNA' ~ 'Mt rRNA',
                                   . == 'misc_RNA' ~ 'misc RNA',
                                   . == 'IG_C_pseudogene' ~ 'IG C pseudogene',
                                   . == 'IG_pseudogene' ~ 'IG pseudogene',
                                   . == 'sense_overlapping' ~ 'sense overlapping',
                                   . == 'TR_V_pseudogene' ~ 'TR V pseudogene',
                                   . == 'TR_D_gene' ~ 'TR D gene',
                                   . == 'IG_C_gene' ~ 'IG C gene',
                                   TRUE ~.))
unique(martTable$Biotype)
row.names(martTable) <- martTable[,1]

#Verify that the master BioMart annotation table corresponds with the order of genes in the DGE object
sum(DGE$genes$genes != row.names(martTable))
#Add BioMart annotation information to the DGE object
DGE$genes$Symbol <- martTable$Symbol
DGE$genes$Description <- martTable$Description
DGE$genes$Biotype <- martTable$Biotype
DGE$genes$Phenotype <- martTable$Phenotype
DGE$genes$Hsap_ID <- martTable$Hsap_ID
DGE$genes$Hsap_IDs <- martTable$Hsap_IDs
DGE$genes$Hsap_symbol <- martTable$Hsap_symbol
DGE$genes$Hsap_symbols <- martTable$Hsap_symbols
DGE$genes$GO_IDs <- martTable$GO_IDs
DGE$genes$GO_Terms <- martTable$GO_Terms
DGE$genes$GO_Definitions <- martTable$GO_Definitions
DGE$genes$GO_Domains <- martTable$GO_Domains

#Write DGE annotation table to file
write.table(DGE$genes, file = paste0(output_filepath,'/Zfish_martTable_Zv10_91.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----7) Additional filtering based on edgeR vignette and (Chen 2016)----
#Determine counts per million corresponding to a minimum of 10-20 read counts 
mean_cpm_gene_count_of_15 <- round(as.numeric(cpm(15, mean(DGE$samples$lib.size))),2)
#Determine the number of samples in the smallest group
num_samples_smallest_group <- min(table(DGE$samples$group))

#Only keep genes that have a minimum of 10-20 counts in at least 4 samples to be included in analysis
keep.exprs <- rowSums(cpm(DGE)> mean_cpm_gene_count_of_15) >= num_samples_smallest_group
DGE <- DGE[keep.exprs, , keep.lib.sizes=FALSE]

#Write background genes to file
write.table(DGE$genes$genes, file = paste0(output_filepath,'/Zfish_background_genes',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#Write human ortholog background genes to file
write.table(DGE$genes$Hsap_ID[DGE$genes$Hsap_ID != ''], file = paste0(output_filepath,'/Metacore_background_genes',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#Write annotation table to file
write.table(DGE$genes, file = paste0(output_filepath,'/DGE_genes_pheno',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----8) Normalize libraries----
#Normalize libraries using the TMM method
DGE <- calcNormFactors(DGE) #use cpm(DGE) to get normalized count table
#Write normalization documentation file
writeLines(capture.output({

    cat('Script filepath:', rstudioapi::getSourceEditorContext()$path, '\n\n')
    cat('Skewer filepath:', skewer_filepath, '\n\n')
    cat('HISAT2 server filepath:', hisat2_filepath, '\n\n')
    cat('HTSeq-count server filepath: ', htseq_filepath, '\n\n')
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
rltdata <- cpm(y = DGE, normalized.lib.sizes = TRUE) 
#Multiply normalized cpm numbers by library size per million to get count values normalized by library size
for (bio_sample in colnames(rltdata)){
    rltdata[,bio_sample] <- (rltdata[,bio_sample])*((DGE$samples[bio_sample,'lib.size']/1000000))
}
#Assign TMM-normalized counts to the 'normalized_counts' slot of the DGE object
DGE$normalized_counts <- rltdata
#Regular log normalization of count data (corrects for heteroskedasticy)
rltdata <- rlogTransformation(round(DGE$normalized_counts)) 
#Add gene labels to TMM-normalized, rlog transformed data
row.names(rltdata) <- row.names(DGE$counts)
#Assign TMM-normalized, rlog transformed counts to the 'rlog_normalized' slot of the DGE object
DGE$rlog_transformed <- rltdata
#Write TMM-normalized, rlog transformed count table to file
write.table(cbind(data.frame('Ensembl_ID' = row.names(DGE$rlog_transformed)), DGE$rlog_transformed), file = paste0(output_filepath,'/DGE_rlog_TMM_normalized_counts',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----9) PCA and interactive MDS plots----
#Create interactive MDS plot of all genes
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

#Create interactive MDS plot of the top 500 most variable genes using pairwise selection between each gene
glMDSPlot(cpm(DGE, log = TRUE), top = 500, 
          gene.selection = 'pairwise',
          col = colors[DGE$samples$group], 
          pch = pch[DGE$samples$group], dim = c(1,2), 
          groups=DGE$samples[,c('group','Chemical_ID','Chemical_Name','Lane')],
          main = paste0('MDS plot of 500 genes with the most variance (pairwise selection)', appender), 
          html = paste0('MDS_plot_500_pairwise',identifier),
          labels = DGE$samples$subID,
          launch = FALSE,
          path = output_filepath,
          folder = 'Interactive_Plots')

#Create interactive MDS plot of the top 500 most variable genes using common gene selection among all genes
glMDSPlot(cpm(DGE, log = TRUE), top = 500, 
          gene.selection = 'common',
          col = colors[DGE$samples$group], 
          pch = pch[DGE$samples$group], dim = c(1,2), 
          groups=DGE$samples[,c('group','Chemical_ID','Chemical_Name','Lane')],
          main = paste0('MDS plot of 500 genes with the most variance (common selection) ', appender), 
          html = paste0('MDS_plot_500_common',identifier),
          labels = DGE$samples$subID,
          launch = FALSE,
          path = output_filepath,
          folder = 'Interactive_Plots')
.jcall("java/lang/System", method = "gc") #helps with java memory problems


##----Create PCA plots----
#Perform principle components analysis on the Tmm-normalized, rlog-transformed data
pca <- prcomp(t(rltdata))
#Extract principle component proportion of variance
sumpca<-data.frame(summary(pca)$importance)

#Make colorblind-friendly color palette for each chemical
cat_12 <- colorschemes$Categorical.12
color_palette <- cat_12[c(7,1:5,8:12)]
show_col(color_palette)

#Create PCA plot of all genes
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


##----Create PCA plots of most variable genes----
#Identify the most variable genes as determined by row variance
topVarGenes <- order(rowVars(rltdata), decreasing = TRUE) 
#Create a table of TMM-normalized, rlog-transformed genes
gen_mat <- as.data.frame(rltdata)
#Reorder the gene table by the most variable genes
gen_mat <- gen_mat[topVarGenes,]
#Scale the transformed count data
gen_mat <- t(scale(t(gen_mat)))

#Perform PCA analysis on the top 500 genes based on variance
pca_500 <- prcomp(t(gen_mat[1:500,])) 
#Extract principle component proportion of variance for the top 500 genes
sumpca_500<-data.frame(summary(pca_500)$importance)

#Create PCA plot of top 500 genes sorted by row variance
PCA_plot_500 <- ggplot(data = data.frame(unlist(pca_500$x[,1:2])), 
                       aes(x = PC1, y = PC2)) + #graph first two PCs
    geom_point(aes(fill = factor(DGE$samples$group)),
               size = 3.5,
               shape = 21,
               stroke = 0.75) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black", size = 0.2) +
    scale_fill_manual(values = color_palette) +
    xlab(paste('PC1 (', 
               round(sumpca_500['Proportion of Variance','PC1']*100),
               '%)', sep = '')) +
    ylab(paste('PC2 (',
               round(sumpca_500['Proportion of Variance','PC2']*100),
               '%)', sep = '')) +
    ggtitle(paste0('PCA plot of TMM-normalized, rlog-transformed,\ntop 500 genes based on gene variance', appender)) +
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
PCA_plot_500

#Write PCA plot of the top 500 most variable genes to PDF
pdf(paste0(output_filepath,'/PCA_plot_500mv',identifier,'.pdf'), width = 6, height = 4)
PCA_plot_500
dev.off()

#Create PCA plot of the top 500 most variable genes with abbreviated sample IDs added to the plot
PCA_plot_500_ids <- ggplot(data = data.frame(unlist(pca_500$x[,1:2])), 
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
               round(sumpca_500['Proportion of Variance','PC1']*100),
               '%)', sep = '')) +
    ylab(paste('PC2 (',
               round(sumpca_500['Proportion of Variance','PC2']*100),
               '%)', sep = '')) +
    ggtitle(paste0('PCA plot of TMM-normalized, rlog-transformed,\ntop 500 genes based on gene variance', appender)) +
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
PCA_plot_500_ids

#Write PCA plot of the top 500 most variable genes to PDF
pdf(paste0(output_filepath,'/PCA_plot_500mv_ids',identifier,'.pdf'), width = 8, height = 6)
PCA_plot_500_ids
dev.off()


##----10) Make the design matrix----
design <- model.matrix(~0 + DGE$samples$group) #the zero means that there is no intercept in the model, which helps with contrasts in later steps
#Assign sample names to the rows of the design matrix
rownames(design) <- colnames(DGE)
#Assign group names to the columns of the design matrix
colnames(design) <- levels(DGE$samples$group)
#Write the design matrix table to documentation file
write.table(cbind(data.frame('Sample' = row.names(design)), as.data.frame(design)), file = paste0(output_filepath,'/Design_matrix',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----11) Estimate dispersion----
#Estimate common, trended, and tagwise negative binomial dispersions by weighted likelihood empirical Bayes
DGE <- estimateDisp(DGE, design, robust = TRUE)
plotBCV(DGE, main = 'Dispersion estimation')


##----12) Fit data to the design model----
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data
fit <- glmQLFit(DGE, design, robust = TRUE)
plotQLDisp(glmfit = fit, main = 'Quasi-liklihood Dispersion')


##----13) Design contrast matrix to extract comparisons of interest----
my_contrasts <- makeContrasts(`BDE47-Con` = BDE47 - Con,
                              `IPP-Con` = IPP - Con,
                              `TBBPA-Con` = TBBPA - Con,
                              `TBBPA_DBPE-Con` = TBBPA_DBPE - Con,
                              `TBPH-Con` = TBPH - Con,
                              `TCEP-Con` = TCEP - Con,
                              `TCPP-Con` = TCPP - Con,
                              `TDBPP-Con` = TDBPP - Con,
                              `TiBP-Con` = TiBP - Con,
                              `TPP-Con` = TPP - Con,
                              levels = design)
#Write contrast table to documentation file
write.table(cbind(data.frame('Group' = row.names(my_contrasts)), as.data.frame(my_contrasts)), file = paste0(output_filepath,'/Contrasts',identifier,'.txt'), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


##----14) Set up ComplexHeatmap parameters----
#Set heatmap global text parameters
ht_opt('heatmap_row_names_gp' = gpar(fontfamily = 'sans'),
       'heatmap_row_title_gp' = gpar(fontfamily = 'sans'),
       'legend_title_gp' = gpar(fontfamily = 'sans', fontsize = 11),
       'legend_labels_gp' = gpar(fontfamily = 'sans', fontsize = 9),
       'heatmap_column_names_gp' = gpar(fontfamily = 'sans'),
       'heatmap_column_title_gp' = gpar(fontfamily = 'sans'))#END ht_opt(RESET = TRUE)

#Set up color schemes for row annotation information
show_col(colorschemes$Categorical.12) #colorblind-friendly color palette
Biotype_colors_de <- c('antisense RNA' = I(colorschemes$Categorical.12[12]),
                        'lincRNA' = I(colorschemes$Categorical.12[5]),
                        'miRNA' = I(colorschemes$Categorical.12[1]),
                        'misc RNA' = I(colorschemes$Categorical.12[3]),
                        'Mt rRNA' = I(colorschemes$Categorical.12[7]),
                        'Mt tRNA' = I(colorschemes$Categorical.12[7]),
                        'polymorphic pseudogene' = I(colorschemes$Categorical.12[2]),
                        'processed transcript' = I(colorschemes$Categorical.12[10]),
                        'protein coding' = I(colorschemes$Categorical.12[9]),
                        'pseudogene' = I(colorschemes$Categorical.12[8]),
                        'rRNA' = I(colorschemes$Categorical.12[6]),
                        'snRNA' = I(colorschemes$Categorical.12[4]),
                        'TEC' = I(colorschemes$Categorical.12[11]),
                        'unprocessed pseudogene' = I(colorschemes$Categorical.12[1]))

#Set up column annotation information for treatment group
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
show_col(chemical_colors)

#Create column annotation information object for 100 genes
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

#Create column annotation information object for 30 genes
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
                                   height = unit(3, "mm"),
                                   simple_anno_size_adjust = TRUE
)#end column annotation


##----15) Find significant genes for each comparison with a fold change threshold of log2(1.5)----
#Specify the desired fold change
fold_change = 1.5
#Specify the desired FDR-adjusted p-value cutoff
FDR_cutoff <- 0.05

#Create table for recording all significantly expressed genes
all_sig_genes <- data.frame()
#Set up table to record each significantly expressed gene per group
sig_gene_table <- data.frame(matrix(NA, nrow=6158, ncol=4*max(table(DGE$samples$group))))
#Assign the initial column number outside the for loop
sig_gene_table_column_number <- 1

#Test for differential expression and create images and tables for each experimental group
for(comparison in colnames(my_contrasts)){
    Sys.sleep(.01)
    print(comparison)
    #Create ID for use in sig_gene_table
    id <- strsplit(x = comparison, split = '-')[[1]][1]
    #Make a label for each contrast
    cutoff_contrast <- paste0(comparison, '_log2(',fold_change,')')
    
    ##Necessary to add FDR values of all genes to MDS and Volcano plots##
    #Test for differential expression with a log(fold change) of 0 to include points in MDS and Volcano plots that do not make the fold-change cutoff
    comparisons <- glmTreat(glmfit = fit, 
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
    #Save the table with original FDR without a fold-change cutoff
    FDR <- sig_table #necessary to add all FDR values to MD plots
    ##Above is necessary to add FDR values of all genes to MD and Volcano plots##
    
    ##----Create output folders for each contrast----
    if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast))){dir.create(paste0(output_filepath,'/', cutoff_contrast))}
    if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))){dir.create(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))}
    #Write table of all differentially expressed genes without a fold-change cutoff to file
    write.table(x = FDR, file = paste0(output_filepath,'/',  cutoff_contrast,'/', comparison, '_all_log_ratios_no_cutoff',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
    #Create table for use with MetaCore
    metacore_table <- FDR[FDR$Hsap_ID != '', c('Hsap_ID','logFC','FDR')]
    #Write all differentially expressed human orthologs to text file
    write.table(x = metacore_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_Metacore_all_genes_no_cutoff',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
    
    
    if(fold_change != 1) {
        ##--Test for differential expression relative to a threshold----
        comparisons <- glmTreat(glmfit = fit, 
                                contrast = my_contrasts[,comparison],
                                lfc = log2(fold_change))
        
        #Extract the all differentially expressed genes that pass a fold change threshold of log2(fold_change)
        sig_table <- data.frame(topTags(comparisons,
                                        n = Inf, adjust.method = "BH",
                                        p.value = 1))
        
        #----Add new FDR values to original differential expression table to take into account a fold-change cutoff
        #Arrange significant gene table by gene id
        sig_table <- sig_table[order(sig_table$genes),]
        FDR <- FDR[order(FDR$genes),]
        
        #----Add new FDR values to original differential expression table to take into account a fold-change cutoff
        for (gene_id in row.names(FDR)){
            if((sig_table[gene_id,'logFC'] >= log2(fold_change)) | (sig_table[gene_id,'logFC'] <= -log2(fold_change))){
                FDR[gene_id, 'FDR'] <- sig_table[gene_id,'FDR']
                FDR[gene_id,'PValue'] <- sig_table[gene_id,'PValue']
                FDR[gene_id, 'logFC'] <- sig_table[gene_id, 'logFC']
                FDR[gene_id, 'logCPM'] <- sig_table[gene_id, 'logCPM']
            }
        }
    }
    #Calculate negative log10(FDR-adjusted p-value) 
    FDR$negLogPval <- -log10(FDR$FDR)
    
    #Extract the significant differentially expressed genes that pass the desired Benjamini-Hochberg FDR cutoff of 0.05
    sig_table <- data.frame(topTags(comparisons,
                                    n = Inf, adjust.method = "BH",
                                    p.value = FDR_cutoff))
    
    if(ncol(sig_table) > 0){
        #Create data objects for each DGELRT object created in the loop
        assign(cutoff_contrast, comparisons)
        
        #Add negative log10(FDR) values to the significant gene table
        sig_table$negLogPval <- -log10(sig_table$FDR)
        #Reorder table by FDR
        sig_table <- sig_table[order(sig_table$FDR, decreasing = FALSE),] 
        #Create data table for each table of significant differentially expressed genes
        assign(paste(cutoff_contrast, 'sigGenes', sep = '_'), sig_table)
        #Add significant differentially expressed genes to the table of all differentially expressed genes
        all_sig_genes <- rbind(all_sig_genes, sig_table) 
        
        #Add the significantly differentially expressed genes of the comparison to the sig_gene_table
        sig_gene_table[1:nrow(sig_table), sig_gene_table_column_number] <- sig_table$genes
        #Name the column by id created from comparison name, i.e. 'IPP genes'
        names(sig_gene_table)[sig_gene_table_column_number] <- paste0(id, ' genes')
        #Advance to the next column number
        sig_gene_table_column_number <- sig_gene_table_column_number + 1
        #Add significant gene symbol information to the table
        sig_gene_table[1:nrow(sig_table), sig_gene_table_column_number] <- sig_table$Symbol
        #Name the column by id created from comparison name, i.e. 'IPP Symbol'
        names(sig_gene_table)[sig_gene_table_column_number] <- paste0(id, ' Symbol')
        #Advance to the next column number
        sig_gene_table_column_number <- sig_gene_table_column_number + 1
        #Add significant gene log2(fold change) values to the sig_gene_table
        sig_gene_table[1:nrow(sig_table), sig_gene_table_column_number] <- sig_table$logFC
        #Name the column by id created from comparison name, i.e. 'IPP logFC'
        names(sig_gene_table)[sig_gene_table_column_number] <- paste0(id, ' logFC')
        #Advance to the next column number
        sig_gene_table_column_number <- sig_gene_table_column_number + 1        
        #Add significant gene FDR values to the sig_gene_table
        sig_gene_table[1:nrow(sig_table), sig_gene_table_column_number] <- sig_table$FDR
        #Name the column by id created from comparison name, i.e. 'IPP FDR'
        names(sig_gene_table)[sig_gene_table_column_number] <- paste0(id, ' FDR')
        #Advance to the next column number
        sig_gene_table_column_number <- sig_gene_table_column_number + 1 
        
        #Write table of significant differentially expressed genes to text file
        write.table(x = sig_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write table of significant differentially expressed genes to Excel file
        xlsx::write.xlsx(x = sig_table, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sig_log_ratios',identifier,'.xlsx'), row.names = FALSE)
        .jcall("java/lang/System", method = "gc") #helps with java memory problems
        
        #Create table for use with MetaCore
        metacore_table <- sig_table[sig_table$Hsap_ID != '', c('Hsap_ID','logFC','FDR')]
        #Write significanly differentially expressed upregulated human orthologs to text file
        write.table(x = metacore_table[(metacore_table$logFC > 0),],file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'UR_Metacore_sigGenes',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write significanly differentially expressed downregulated human orthologs to text file
        write.table(x = metacore_table[(metacore_table$logFC < 0),],file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'DR_Metacore_sigGenes',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write significanly differentially expressed human orthologs to text file
        write.table(x = metacore_table,file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'all_Metacore_sigGenes',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed genes to text file
        write.table(x = sig_table$genes, file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_sigGene_list',identifier,'.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed upregulated genes to text file
        write.table(x = sig_table$genes[sig_table$logFC > 0], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'UR_sigGene_list',identifier,'.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        #Write list of significant differentially expressed downregulated genes to text file
        write.table(x = sig_table$genes[sig_table$logFC < 0], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, 'DR_sigGene_list',identifier,'.txt'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        #Create tables for use with the Bioinformatics Resource Manager
        if(fold_change == 1){
            #Write significanly differentially expressed downregulated genes to text file for transfer to BRM
            write.table(x = sig_table[(sig_table$logFC < 0),c('genes','Symbol','logFC','logCPM','F','PValue','FDR')], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_BRM_DR_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
            #Write significanly differentially expressed upregulated genes to text file for transfer to BRM
            write.table(x = sig_table[(sig_table$logFC > 0),c('genes','Symbol','logFC','logCPM','F','PValue','FDR')], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_BRM_UR_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
            #Write all significanly differentially expressed  genes to text file for transfer to BRM
            write.table(x = sig_table[,c('genes','Symbol','logFC','logCPM','F','PValue','FDR')], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_BRM_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        }else{
            #Write significanly differentially expressed downregulated genes to text file for transfer to BRM
            write.table(x = sig_table[(sig_table$logFC < 0),c('genes','Symbol','logFC','logCPM','PValue','FDR')], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_BRM_DR_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
            #Write significanly differentially expressed upregulated genes to text file for transfer to BRM
            write.table(x = sig_table[(sig_table$logFC > 0),c('genes','Symbol','logFC','logCPM','PValue','FDR')], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_BRM_UR_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
            #Write all significanly differentially expressed genes to text file for transfer to BRM
            write.table(x = sig_table[,c('genes','Symbol','logFC','logCPM','PValue','FDR')], file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast, '_BRM_sig_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
            
        }
    }
    
    ##----Create heatmaps and interactive plots of top significantly differentially expressed genes----
    if(nrow(sig_table) >= 2){
        #Order genes based on their absolute log2(fold change)
        heatmap_selector_FC <- sig_table[order(abs(sig_table$logFC), decreasing = TRUE),]
        #Order genes based on their FDR-adjusted p-values
        heatmap_selector_FDR <- sig_table[order(sig_table$FDR),] 
        
        ##----Create heatmap of top 30 genes sorted by fold change----
        #Select the top 30 genes
        heatmap_selector_FC30 <- heatmap_selector_FC[1:30,'genes']
        #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 genes
        heatmap_selector_FC30 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FC30
        #Select data from the TMM-normalized, rlog transformed count table matching the top 30 genes
        rlogdataFC30 <- data.frame(DGE$rlog_transformed[heatmap_selector_FC30,])
        #Scale the transformed count data
        rlogdataFC30 <- t(scale(t(rlogdataFC30)))
        
        #Select data from the annotation table matching the top 30 genes
        phenodataFC30 <- DGE$genes[heatmap_selector_FC30, 'Biotype',drop = FALSE]
        #Label phenotype data with gene symbols
        row.names(phenodataFC30) <- DGE$genes$Symbol[heatmap_selector_FC30]
        #Set row dendrogram line thickness
        row_dend_FC = as.dendrogram(hclust(dist(rlogdataFC30))) %>% set('branches_lwd', 1.5)
        #Set column dendrogram line thickness
        col_dend_FC = as.dendrogram(hclust(dist(t(rlogdataFC30)))) %>% set('branches_lwd', 1.5)
        
        #Set up row annotation legend information
        row_key <- names(table(as.character(phenodataFC30$Biotype)))
        ##If there are 30 or fewer total significant differentially expressed genes, the following heatmaps are created
        if(nrow(heatmap_selector_FC) <= 30) {
            #Identify the number of significant differentially expressed genes
            number_of_genes <- length(phenodataFC30$Biotype)
            ##Create main ComplexHeatmap layer of all significant differentially expressed genes
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30',
                            column_title = paste0(cutoff_contrast, ':\nAll ',number_of_genes,' differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 10),
                            row_dend_side = 'left',
                            row_dend_width = unit(10, 'mm'),
                            show_row_dend = TRUE,
                            column_dend_side = 'top',
                            column_dend_height = unit(7.5,'mm'),
                            show_column_dend = TRUE,
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
            
            ##Create row annotation ComplexHeatmap layer
            rmap <- Heatmap(matrix = as.matrix(phenodataFC30),
                            col = Biotype_colors_de,
                            name = 'Biotype', 
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 6,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 8),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = (
                                Biotype = list(nrow = if(length(row_key) > 6){2}else{1}, 
                                               title = "Biotype", 
                                               title_gp = gpar(fontsize = 8, 
                                                               fontface = 'plain'),
                                               title_position = "topleft",
                                               labels_gp = gpar(fontsize = 6),
                                               legend_direction = 'horizontal', 
                                               fontfamily = 'sans')),
                            width = unit(2.5,'mm')
            )#END rmap
            
            ##----Draw the combined heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_all',number_of_genes,identifier,'.pdf'), width = 6, height = 5)
            draw(hmap + rmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##Create main ComplexHeatmap layer with dimensions for PowerPoint
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30',
                            column_title = paste0(cutoff_contrast, ':\nAll ',number_of_genes,' differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14.5),
                            cluster_rows = row_dend_FC,
                            row_dend_side = 'left',
                            row_dend_width = unit(15, "mm"),
                            show_row_dend = TRUE,
                            cluster_columns = col_dend_FC,
                            column_dend_side = 'top',
                            column_dend_height = unit(15,'mm'),
                            show_column_dend = TRUE,
                            show_row_names = FALSE,
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
            
            
            ##Create row annotation ComplexHeatmap layer for PowerPoint
            rmap <- Heatmap(matrix = as.matrix(phenodataFC30),
                            col = Biotype_colors_de,
                            name = 'Biotype', 
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 10.5,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 12.5),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = (
                                Biotype = list(nrow = if(length(row_key) > 6){2}else{1}, 
                                               title = "Biotype", 
                                               title_gp = gpar(fontsize = 12.5, 
                                                               fontface = 'plain'),
                                               title_position = "topleft",
                                               labels_gp = gpar(fontsize = 10.5, 
                                                                fontfamily = 'sans'),
                                               by_row = TRUE,
                                               legend_direction = 'horizontal')),
                            width = unit(3,'mm')
            )#END rmap
            
            ##----Draw the combined heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_all',number_of_genes,'_PPT',identifier,'.pdf'), width = 9, height = 7.5)
            draw(hmap + rmap, padding = unit(c(2,2,2,2), unit = 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
        }else{
            ##Create main ComplexHeatmap layer of the top 30 genes sorted by absolute fold change
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
            
            ##Create row annotation ComplexHeatmap layer of the top 30 genes sorted by absolute fold change
            rmap <- Heatmap(matrix = as.matrix(phenodataFC30),
                            col = Biotype_colors_de,
                            name = 'Biotype', 
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 6,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 8),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = (
                                Biotype = list(nrow = if(length(row_key) > 6){2}else{1}, 
                                               title = "Biotype", 
                                               title_gp = gpar(fontsize = 8, 
                                                               fontface = 'plain'),
                                               title_position = "topleft",
                                               labels_gp = gpar(fontsize = 6),
                                               legend_direction = 'horizontal', 
                                               fontfamily = 'sans')),
                            width = unit(2.5,'mm')
            )#END rmap
            
            ##----Draw the combined heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FC',identifier,'.pdf'), width = 6, height = 5)
            draw(hmap + rmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##Create main ComplexHeatmap layer with dimensions for PowerPoint for the top 30 genes sorted by absolute fold change 
            hmap <- Heatmap(matrix = rlogdataFC30, 
                            name = 'top_30_FC',
                            column_title = paste0(cutoff_contrast, ': Top 30 fold change \nof all differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14.5),
                            cluster_rows = row_dend_FC,
                            row_dend_side = 'left',
                            row_dend_width = unit(15, "mm"),
                            show_row_dend = TRUE,
                            cluster_columns = col_dend_FC,
                            column_dend_side = 'top',
                            column_dend_height = unit(15,'mm'),
                            show_column_dend = TRUE,
                            show_row_names = FALSE,
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
            
            ##Create row annotation ComplexHeatmap layer for PowerPoint for the top 30 genes sorted by absolute fold change
            rmap <- Heatmap(matrix = as.matrix(phenodataFC30),
                            col = Biotype_colors_de,
                            name = 'Biotype', 
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 10.5,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 12.5),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = (
                                Biotype = list(nrow = if(length(row_key) > 6){2}else{1}, 
                                               title = "Biotype", 
                                               title_gp = gpar(fontsize = 12.5, 
                                                               fontface = 'plain'),
                                               title_position = "topleft",
                                               labels_gp = gpar(fontsize = 10.5, 
                                                                fontfamily = 'sans'),
                                               by_row = TRUE,
                                               legend_direction = 'horizontal')),
                            width = unit(3,'mm')
            )#END rmap
            
            ##----Draw the combined heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FC_PPT',identifier,'.pdf'), width = 9, height = 7.5)
            draw(hmap + rmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##----Create heatmap of top 30 genes sorted by FDR-adjusted p-value----
            #Select the top 30 genes
            heatmap_selector_FDR30 <- heatmap_selector_FDR[1:30,'genes']
            #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 30 genes
            heatmap_selector_FDR30 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FDR30
            #Select data from the TMM-normalized, rlog transformed count table matching the top 30 genes
            rlogdataFDR30 <- data.frame(DGE$rlog_transformed[heatmap_selector_FDR30,])
            #Scale the transformed count data
            rlogdataFDR30 <- t(scale(t(rlogdataFDR30)))
            #Select data from the annotation table matching the top 30 genes
            phenodataFDR30 <- DGE$genes[heatmap_selector_FDR30, 'Biotype',drop = FALSE]
            #Label phenotype data with gene symbols
            row.names(phenodataFDR30) <- DGE$genes$Symbol[heatmap_selector_FDR30]
            
            #Set up row annotation legend information
            row_key <- names(table(as.character(phenodataFDR30$Biotype)))
            #Set row dendrogram line thickness
            row_dend_FDR = as.dendrogram(hclust(dist(rlogdataFDR30))) %>% set('branches_lwd', 1.5)
            #Set column dendrogram line thickness
            col_dend_FDR = as.dendrogram(hclust(dist(t(rlogdataFDR30)))) %>% set('branches_lwd', 1.5)
            
            ##Create main ComplexHeatmap layer of the top 30 genes sorted by FDR-adjusted p-value
            hmap <- Heatmap(matrix = as.matrix(rlogdataFDR30), 
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
            
            ##Create row annotation ComplexHeatmap layer of the top 30 genes sorted by FDR-adjusted p-value
            rmap <- Heatmap(matrix = as.matrix(phenodataFDR30), 
                            col = Biotype_colors_de,
                            name = 'Biotype', 
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 6,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 8),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = (
                                Biotype = list(nrow = if(length(row_key) > 6){2}else{1}, 
                                               title = "Biotype", 
                                               title_gp = gpar(fontsize = 8, 
                                                               fontface = 'plain'),
                                               title_position = "topleft",
                                               labels_gp = gpar(fontsize = 6),
                                               legend_direction = 'horizontal', 
                                               fontfamily = 'sans')),
                            width = unit(2.5,'mm')
            )#END rmap
            
            ##----Draw the combined heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FDR',identifier,'.pdf'), width = 6, height = 5)
            draw(hmap + rmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            
            ##Create main ComplexHeatmap layer with dimensions for PowerPoint for the top 30 genes sorted by FDR-adjusted p-value
            hmap <- Heatmap(matrix = as.matrix(rlogdataFDR30), 
                            name = 'top_30_FDR',
                            column_title = paste0(cutoff_contrast, ': Top 30 FDR-adjusted p-value \nof all differentially expressed genes',appender), 
                            column_title_side = 'top',
                            column_title_gp = gpar(fontsize = 14.5),
                            cluster_rows = row_dend_FDR,
                            row_dend_side = 'left',
                            row_dend_width = unit(15, "mm"),
                            show_row_dend = TRUE,
                            cluster_columns = col_dend_FDR,
                            column_dend_side = 'top',
                            column_dend_height = unit(15,'mm'),
                            show_column_dend = TRUE,
                            show_row_names = FALSE,
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
            
            ##Create row annotation ComplexHeatmap layer for PowerPoint for the top 30 genes sorted by FDR-adjusted p-value
            rmap <- Heatmap(matrix = as.matrix(phenodataFDR30), 
                            col = Biotype_colors_de,
                            name = 'Biotype', 
                            row_names_side = 'right',
                            show_row_names = TRUE,
                            row_names_gp = gpar(fontsize = 10.5,
                                                fontface = 'italic'),
                            column_names_gp = gpar(fontsize = 12.5),
                            show_heatmap_legend = TRUE,
                            heatmap_legend_param = (
                                Biotype = list(nrow = if(length(row_key) > 6){2}else{1}, 
                                               title = "Biotype", 
                                               title_gp = gpar(fontsize = 12.5, 
                                                               fontface = 'plain'),
                                               title_position = "topleft",
                                               labels_gp = gpar(fontsize = 10.5, 
                                                                fontfamily = 'sans'),
                                               by_row = TRUE,
                                               legend_direction = 'horizontal')),
                            width = unit(3,'mm')
            )#END rmap
            
            ##----Draw the combined heatmap and save as a PDF---
            pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top30FDR_PPT',identifier,'.pdf'), width = 9, height = 7.5)
            draw(hmap + rmap, padding = unit(c(2, 2, 2, 2), "mm"), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
            dev.off()
            
            ##----Create heatmap of top 100 genes sorted by absolute fold change----
            #Select the top 100 genes
            heatmap_selector_FC100 <- heatmap_selector_FC[1:100,'genes']
            #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 100 genes
            heatmap_selector_FC100 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FC100
            #Select data from the TMM-normalized, rlog transformed count table matching the top 100 genes
            rlogdataFC100 <- data.frame(DGE$rlog_transformed[heatmap_selector_FC100,])
            #Scale the transformed count data
            rlogdataFC100 <- t(scale(t(rlogdataFC100)))
            #Select data from the annotation table matching the top 100 genes
            phenodataFC100 <- DGE$genes[heatmap_selector_FC100, 'Biotype',drop = FALSE]
            #Label phenotype data with gene symbols
            row.names(phenodataFC100) <- DGE$genes$Symbol[heatmap_selector_FC100]
            
            #Set up row annotation legend information
            row_key100 <- names(table(as.character(phenodataFC100$Biotype)))
            ##If there are 100 or fewer total significant differentially expressed genes, the following heatmaps are created
            if(nrow(heatmap_selector_FC) <= 100){
                #Identify the number of significant differentially expressed genes
                number_of_genes <- length(phenodataFC100$Biotype)
                ##Create main ComplexHeatmap layer of all significant differentially expressed genes
                hmap <- Heatmap(matrix = as.matrix(rlogdataFC100), 
                                name = 'top_100_fc',
                                column_title = paste0(cutoff_contrast, ':\nAll ',number_of_genes,' differentially expressed genes',appender), 
                                column_title_side = 'top',
                                column_title_gp = gpar(fontsize = 14),
                                row_dend_side = 'left',
                                row_dend_width = unit(13, "mm"),
                                show_row_dend = TRUE,
                                column_dend_side = 'top',
                                column_dend_height = unit(8,'mm'),
                                show_column_dend = TRUE,
                                show_row_names = FALSE,
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
                
                ##Create row annotation ComplexHeatmap layer
                rmap <- Heatmap(matrix = as.matrix(phenodataFC100), 
                                col = Biotype_colors_de,
                                name = 'Biotype', 
                                row_names_side = 'right',
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 5,
                                                    fontface = 'italic'),
                                column_names_gp = gpar(fontsize = 11),
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = (
                                    Biotype = list(nrow = if(length(row_key100) > 6){2}else{1}, 
                                                   title = "Biotype", 
                                                   title_gp = gpar(fontsize = 11, 
                                                                   fontface = 'plain'),
                                                   title_position = "topleft",
                                                   labels_gp = gpar(fontsize = 9, 
                                                                    fontfamily = 'sans'),
                                                   by_row = TRUE,
                                                   legend_direction = 'horizontal')),
                                width = unit(3,'mm')
                )#END rmap
                
                ##----Draw the combined heatmap and save as a PDF---
                pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_all',number_of_genes,identifier,'.pdf'), width = 8.5, height = 11)
                draw(hmap + rmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
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
                                show_row_names = FALSE,
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
                
                ##Create row annotation ComplexHeatmap layer of the top 100 genes sorted by absolute fold change
                rmap <- Heatmap(matrix = as.matrix(phenodataFC100), 
                                col = Biotype_colors_de,
                                name = 'Biotype', 
                                row_names_side = 'right',
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 5,
                                                    fontface = 'italic'),
                                column_names_gp = gpar(fontsize = 11),
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = (
                                    Biotype = list(nrow = if(length(row_key100) > 6){2}else{1}, 
                                                   title = "Biotype", 
                                                   title_gp = gpar(fontsize = 11, 
                                                                   fontface = 'plain'),
                                                   title_position = "topleft",
                                                   labels_gp = gpar(fontsize = 9, 
                                                                    fontfamily = 'sans'),
                                                   by_row = TRUE,
                                                   legend_direction = 'horizontal')),
                                width = unit(3,'mm')
                )#END rmap
                
                ##----Draw the combined heatmap and save as a PDF---
                pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top100FC',identifier,'.pdf'), width = 8.5, height = 11)
                draw(hmap + rmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
                dev.off()
                
                ##----Create heatmap of top 100 genes sorted by FDR-adjusted p-value----
                #Select the top 100 genes
                heatmap_selector_FDR100 <- heatmap_selector_FDR[1:100,'genes']
                #Determine the data from the TMM-normalized, rlog transformed count table that matches the top 100 genes
                heatmap_selector_FDR100 <- row.names(DGE$rlog_transformed) %in% heatmap_selector_FDR100
                #Select data from the TMM-normalized, rlog transformed count table matching the top 100 genes
                rlogdataFDR100 <- data.frame(DGE$rlog_transformed[heatmap_selector_FDR100,])
                #Scale the transformed count data
                rlogdataFDR100 <- t(scale(t(rlogdataFDR100)))
                #Select data from the annotation table matching the top 100 genes
                phenodataFDR100 <- DGE$genes[heatmap_selector_FDR100, 'Biotype',drop = FALSE]
                #Label phenotype data with gene symbols
                row.names(phenodataFDR100) <- DGE$genes$Symbol[heatmap_selector_FDR100]
                
                ##Set up row annotation legend information
                row_key100 <- names(table(as.character(phenodataFDR100$Biotype)))
                
                ##Create main ComplexHeatmap layer of the top 100 genes sorted by FDR-adjusted p-value
                hmap <- Heatmap(matrix = as.matrix(rlogdataFDR100), 
                                name = 'top_100_FDR',
                                column_title = paste0(cutoff_contrast, ': Top 100 FDR-adjusted p-value \nof all differentially expressed genes',appender), 
                                column_title_side = 'top',
                                column_title_gp = gpar(fontsize = 14),
                                row_dend_side = 'left',
                                row_dend_width = unit(13, "mm"),
                                show_row_dend = TRUE,
                                column_dend_side = 'top',
                                column_dend_height = unit(8,'mm'),
                                show_column_dend = TRUE,
                                show_row_names = FALSE,
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
                
                
                ##Create row annotation ComplexHeatmap layer of the top 100 genes sorted by FDR-adjusted p-value
                rmap <- Heatmap(matrix = as.matrix(phenodataFDR100), 
                                col = Biotype_colors_de,
                                name = 'Biotype', 
                                row_names_side = 'right',
                                show_row_names = TRUE,
                                row_names_gp = gpar(fontsize = 5,
                                                    fontface = 'italic'),
                                column_names_gp = gpar(fontsize = 11),
                                show_heatmap_legend = TRUE,
                                heatmap_legend_param = (
                                    Biotype = list(nrow = if(length(row_key100) > 6){2}else{1}, 
                                                   title = "Biotype", 
                                                   title_gp = gpar(fontsize = 11, 
                                                                   fontface = 'plain'),
                                                   title_position = "topleft",
                                                   labels_gp = gpar(fontsize = 9, 
                                                                    fontfamily = 'sans'),
                                                   by_row = TRUE,
                                                   legend_direction = 'horizontal')),
                                width = unit(3,'mm')
                )#END rmap
                
                ##----Draw the combined heatmap and save as a PDF---
                pdf(file = paste0(output_filepath,'/',cutoff_contrast,'/', cutoff_contrast,'_heatmap_top100FDR',identifier,'.pdf'), width = 8.5, height = 11)
                draw(hmap + rmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
                dev.off()
                
            }
        }
        

        ##----Generate MD plots of differentially expressed genes----
        #Create directory for Interactive Glimma plots
        if(!dir.exists(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))){dir.create(paste0(output_filepath,'/', cutoff_contrast,'/Interactive_Plots'))}
        #Determine which genes are considered significant based on the entered FDR_cutoff
        dt <- decideTestsDGE(comparisons, p.value = FDR_cutoff)
        #Add GeneID values to the differential expression table 
        comparisons$table$GeneID <- row.names(comparisons$table)
        #Add a gene symbols column to the comparisons differential expression table
        comparisons$table$Symbol <- comparisons$genes$Symbol
        #Add the FDR-adjusted p-values to the differential expression table
        comparisons$table$FDR <- FDR$FDR
        #Add the -log10(FDR) to the differential expression table
        comparisons$table$negLogPval <- FDR$negLogPval
        #Add the p-values to the differential expression table
        comparisons$table$pValue <- FDR$PValue
        #Add the biotype annotation information to the differential expression table
        comparisons$table$Biotype <- comparisons$genes$Biotype
        #Add human Ensembl ID annotation information to the differential expression table
        comparisons$table$Hsap_IDs <- comparisons$genes$Hsap_IDs
        #Add human symbol annotation information to the differential expression table
        comparisons$table$Hsap_symbols <- comparisons$genes$Hsap_symbols
        #Add gene description annotation information to the differential expression table
        comparisons$table$Description <- comparisons$genes$Description
        #Add Ensembl gene phenotype annotation information to the differential expression table
        comparisons$table$Phenotype <- comparisons$genes$Phenotype
        #Add Gene Ontology ID annotation information to the differential expression table
        comparisons$table$GO_IDs <- comparisons$genes$GO_IDs
        #Add GO term annotation information to the differential expression table
        comparisons$table$GO_Terms <- comparisons$genes$GO_Terms
        #Add GO definition annotation information to the differential expression table
        comparisons$table$GO_Definitions <- comparisons$genes$GO_Definitions
        #Add GO domain annotation information to the differential expression table
        comparisons$table$GO_Domains <- comparisons$genes$GO_Domains
        
        if(fold_change > 1){
            write.table(x = comparisons$table, file = paste0(output_filepath,'/',  cutoff_contrast,'/', cutoff_contrast, '_Glimma_log_ratios',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
        }
        
        #Create interactive mean-difference (MD) plot comparing the mean log2(CPM) to the log2(fold change)
        glXYPlot(x = comparisons$table$logCPM,
                 y = comparisons$table$logFC,
                 counts = cpm(DGE$counts, log = TRUE, normalized.lib.sizes = TRUE),
                 groups = DGE$samples$group,
                 samples = DGE$samples$subID,
                 status = dt,
                 anno = comparisons$table,
                 display.columns = c('GeneID','Symbol','pValue', 'FDR', 'logFC','logCPM','Biotype','Hsap_IDs','Hsap_symbols','Description','Phenotype','GO_Terms','GO_Definitions','GO_IDs'),
                 xlab = 'logCPM',
                 ylab = 'logFC',
                 launch = FALSE,
                 path = paste0(output_filepath,'/', cutoff_contrast),
                 folder = 'Interactive_Plots',
                 html = paste0('MD_Plot_', cutoff_contrast,identifier),
                 main = paste0(comparison, ': (', comparisons$comparison, ')<br> fold change cutoff = log<sub>2</sub>(',fold_change,')',appender))
        .jcall("java/lang/System", method = "gc")
        
        #Create interactive Volcano plot comparing the log2(fold change) and the -log10(FDR-adjusted p-value)
        glXYPlot(x = comparisons$table$logFC,
                 y = comparisons$table$negLogPval,
                 counts = cpm(DGE$counts, log = TRUE, normalized.lib.sizes = TRUE),
                 groups = DGE$samples$group,
                 samples = DGE$samples$subID,
                 status = dt,
                 anno = comparisons$table,
                 display.columns = c('GeneID','Symbol','pValue', 'FDR', 'logFC','logCPM','Biotype','Hsap_IDs','Hsap_symbols','Description','Phenotype','GO_Terms','GO_Definitions','GO_IDs'),
                 xlab = 'logFC',
                 ylab = '-log10(p-value)',
                 launch = FALSE,
                 path = paste0(output_filepath,'/', cutoff_contrast),
                 folder = 'Interactive_Plots',
                 html = paste0('Volcano-Plot_', cutoff_contrast, identifier),
                 main = paste0(comparison, ': (', comparisons$comparison, ')<br> fold change cutoff = log<sub>2</sub>(',fold_change,')',appender))
        .jcall("java/lang/System", method = "gc")
    }
}


#Determine the number of empty rows in the combined significant differential expression table
sum(rowSums(is.na(sig_gene_table)) == ncol(sig_gene_table))
#Keep rows that contain at least 1 value
sig_gene_table <- sig_gene_table[(rowSums(is.na(sig_gene_table)) < ncol(sig_gene_table)),]
#Replace NA values with ''
sig_gene_table[is.na(sig_gene_table)] <- ''

#Write the combined significant differential expression table to text file
write.table(x = sig_gene_table, file = paste0(output_filepath,'/All_chemicals_sig_log_ratios','_log2(',fold_change,')',identifier,'.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
#Write the combined significant differential expression table to Excel file
xlsx::write.xlsx(x = sig_gene_table, file = paste0(output_filepath,'/All_chemicals_sig_log_ratios','_log2(',fold_change,')',identifier,'.xlsx'), row.names = FALSE)
.jcall("java/lang/System", method = "gc")

all_sig_genes.o <- all_sig_genes
#Order all significantly differentially expressed genes by their FDR-adjusted p-value
all_sig_genes <- all_sig_genes[order(all_sig_genes$FDR, decreasing = FALSE),]
#Remove duplicated values from the significant gene table
all_sig_genes <- all_sig_genes[!duplicated(all_sig_genes$genes),]


##----16) Heatmap of gene count values from all combined significantly DE genes sorted by FDR-adjusted p-value----
#Select the top 100 genes
heatmap_selectorFDR100 <- row.names(DGE$rlog_transformed) %in% all_sig_genes$genes[1:100]
#Select data from the TMM-normalized, rlog transformed count table matching the top 100 genes
rlogdataFDR100 <- data.frame(DGE$rlog_transformed[heatmap_selectorFDR100,])
#Scale the transformed count data
rlogdataFDR100 <- t(scale(t(rlogdataFDR100)))
#Select data from the annotation table matching the top 100 genes
phenodataFDR100 <- DGE$genes[heatmap_selectorFDR100, 'Biotype',drop = FALSE]
#Label phenotype data with gene symbols
row.names(phenodataFDR100) <- DGE$genes$Symbol[heatmap_selector_FDR100]

##Set up row annotation information
row_key100 <- names(table(as.character(phenodataFDR100$Biotype)))

##Create main ComplexHeatmap layer of the top 100 combined significantly differentially expressed genes sorted by FDR-adjusted p-value
hmap <- Heatmap(matrix = as.matrix(rlogdataFDR100), 
                name = 'top_var_100',
                column_title = paste0('TMM-normalized, rlog-transformed, average count values from top 100 FDR−adjusted p−value\nof combined differentially expressed genes at fold change cutoff of log2(',fold_change,')',appender),

                column_title_side = 'top',
                column_title_gp = gpar(fontsize = 13),
                row_dend_side = 'left',
                row_dend_width = unit(13, "mm"),
                show_row_dend = TRUE,
                column_dend_side = 'top',
                column_dend_height = unit(8,'mm'),
                show_column_dend = TRUE,
                show_row_names = FALSE,
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


##Create row annotation ComplexHeatmap layer of the top 100 combined significantly differentially expressed genes sorted by FDR-adjusted p-value
rmap <- Heatmap(matrix = as.matrix(phenodataFDR100), 
                col = Biotype_colors_de,
                name = 'Biotype', 
                row_names_side = 'right',
                show_row_names = TRUE,
                row_names_gp = gpar(fontsize = 5,
                                    fontface = 'italic'),
                column_names_gp = gpar(fontsize = 11),
                show_heatmap_legend = TRUE,
                heatmap_legend_param = (
                    Biotype = list(nrow = if(length(row_key100) > 6){2}else{1}, 
                                   title = "Biotype", 
                                   title_gp = gpar(fontsize = 11, 
                                                   fontface = 'plain'),
                                   title_position = "topleft",
                                   labels_gp = gpar(fontsize = 9, 
                                                    fontfamily = 'sans'),
                                   by_row = TRUE,
                                   legend_direction = 'horizontal')),
                width = unit(3,'mm')
)#END rmap

##----Draw the combined heatmap and save as a PDF---
pdf(file = paste0(output_filepath, '/Heatmap_100FDR_all_sigGenes_log2(',fold_change,')',identifier,'.pdf'), width = 8.5, height = 11)
draw(hmap + rmap, padding = unit(c(2,2,2,2), 'mm'), gap = unit(1, "mm"), heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom')
dev.off()


##----17) Write session information----
#Name to add to the beginning of file names
name <- 'FRC_RNAseq_analysis'
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
    
}), paste0(output_filepath,'/',name,'_sessionInfo_',Sys.Date(),'.txt'))


#Save a copy of the R Global Environment to work with later
save.image(paste0(output_filepath,'/',name,'_RData_',Sys.Date(),'.RData'))


#Removes all objects from the global environment (can reload them from the .RData file created above)
rm(list = ls())
