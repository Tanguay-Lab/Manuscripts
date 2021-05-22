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
##  (541) 737-6500                                                       ##
##  dunhamcg@gmail.com                                                   ##
###########################################################################


##----1) Install any missing packages----
#Use to install multiple CRAN packages
#install.packages(c('xlsx','dplyr')) 
#Use to install a single CRAN package
#install.packages('xlsx') 


##----2) Load packages and specify file paths----
#List all required packages to load in loop

#List all required packages to load in loop
all_packages <- c('dplyr', 
                  'xlsx',
                  'ggplot2')
#Load all packages at once, ignore the screen print out
sapply(all_packages, require, character.only = TRUE)


#I append this to the beginning of all files to ensure all file names go into whatever folder you specify
output_filepath <- 'C:/Users/dunhamc/Google Drive/R_scripts/Manuscripts/SlincR (2nd paper)'
input_filepath <- 'C:/Users/dunhamc/Google Drive/Tanguay_Lab/Phenotype_Scrn/Cart_stain/08-11-2017/Cart_data_8_11_17_b.xlsx'

##----3) Read in table----
#read in file, specify path to file, tab of Excel workbook, and column header is names in first row
analysis_table <- read.xlsx(file = input_filepath, sheetIndex = 'R_input', header = TRUE, stringsAsFactors = FALSE)

#turn morphant status, treatment, and location into factors and specify control with levels
analysis_table$Morphant <- factor(analysis_table$Morphant, levels = c('ConMO','SlincR_MO'))
analysis_table$Treatment <- factor(analysis_table$Treatment, levels = c('DMSO','TCDD'))
analysis_table$Location <- factor(analysis_table$Location)

#make new column defining groups required for statistical testing below
analysis_table$Group <- factor(paste(analysis_table$Morphant, analysis_table$Treatment, sep = '_'), levels = c('ConMO_DMSO','ConMO_TCDD','SlincR_MO_DMSO','SlincR_MO_TCDD'))


##----4) Scatterplots w/ x and y SE bars----
cart_summary <- analysis_table %>%
    group_by(.dots= c('Morphant','Treatment','Location')) %>%
    summarize(N = n(),
              mean_x = signif(mean(X_coord),digits = 3),
              sd_x = signif(sd(X_coord), digits = 3),
              SE_x = signif(sd(X_coord)/sqrt(n()), digits = 3),
              mean_y = signif(mean(Y_coord), digits =3),
              sd_y = signif(sd(Y_coord), digits = 3),
              SE_y = signif(sd(Y_coord)/sqrt(n()), digits = 3))


##----Position A---- 
pdf(paste0(output_filepath,'/Cart_data_A.pdf'), width = 10, height = 8)
print(ggplot(data = subset(cart_summary, Location == 'A'), aes(x = mean_x, y = mean_y)) +
          geom_point(aes(fill = Treatment, shape = Morphant, size = Morphant), color = 'black', stroke = 0.05) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          scale_y_continuous(name = 'y-axis', expand = c(0,0), limits = c(0,165), position = 'right') +
          scale_x_continuous(name = 'x-axis', expand = c(0,0), limits = c(-255,0), position = 'bottom') +
          scale_shape_manual(values = c(22,23)) +
          scale_size_manual(values = c(9,9)) +
          scale_fill_manual(values = c('blue','red')) +
          geom_errorbar(aes(ymin = mean_y-SE_y, ymax = mean_y + SE_y),
                        width = 5, size = 0.75) +
          geom_errorbarh(aes(xmin = mean_x-SE_x, xmax = mean_x + SE_x),
                         height = 5, size = 0.75) +
          guides(fill = guide_legend(override.aes = list(size = 8, color = c('blue','red'), order = 1),
                                     shape = guide_legend(override.aes = list(size = c(8,8),fill ='black')), order = 2)) +
          ggtitle('Position A') +
          theme(text = element_text(family = 'sans'),
                plot.title = element_text(size = 24, hjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.key = element_rect(fill = 'white'),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 14),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_rect(fill = 'white', color = 'black'),
                panel.border = element_blank(),
                panel.grid = element_blank()))
dev.off()

##----Position B----
pdf(paste0(output_filepath,'/Cart_data_B.pdf'), width = 10, height = 8)
print(ggplot(data = subset(cart_summary, Location == 'B'), aes(x = mean_x, y = mean_y)) +
          geom_point(aes(fill = Treatment, shape = Morphant, size = Morphant), color = 'black', stroke = 0.05) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          scale_y_continuous(name = 'y-axis', expand = c(0,0), limits = c(0,165), position = 'left') +
          scale_x_continuous(name = 'x-axis', expand = c(0,0), limits = c(0,255), position = 'bottom') +
          scale_shape_manual(values = c(22,23)) +
          scale_size_manual(values = c(9,9)) +
          scale_fill_manual(values = c('blue','red')) +
          geom_errorbar(aes(ymin = mean_y-SE_y, ymax = mean_y + SE_y),
                        width = 5, size = 0.75) +
          geom_errorbarh(aes(xmin = mean_x-SE_x, xmax = mean_x + SE_x),
                         height = 5, size = 0.75) +
          guides(fill = guide_legend(override.aes = list(size = 8, color = c('blue','red'), order = 1),
                                     shape = guide_legend(override.aes = list(size = c(8,8),fill ='black')), order = 2)) +
          ggtitle('Position B') +
          theme(text = element_text(family = 'sans'),
                plot.title = element_text(size = 24, hjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.key = element_rect(fill = 'white'),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 14),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_rect(fill = 'white', color = 'black'),
                panel.border = element_blank(),
                panel.grid = element_blank()))
dev.off()

##----Position C----
pdf(paste0(output_filepath,'/Cart_data_C.pdf'), width = 10, height = 8)
print(ggplot(data = subset(cart_summary, Location == 'C'), aes(x = mean_x, y = mean_y)) +
          geom_point(aes(fill = Treatment, shape = Morphant, size = Morphant), color = 'black', stroke = 0.05) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          scale_y_continuous(name = 'y-axis', expand = c(0,0), limits = c(-165,0), position = 'right') +
          scale_x_continuous(name = 'x-axis', expand = c(0,0), limits = c(-255,0), position = 'top') +
          scale_shape_manual(values = c(22,23)) +
          scale_size_manual(values = c(9,9)) +
          scale_fill_manual(values = c('blue','red')) +
          geom_errorbar(aes(ymin = mean_y-SE_y, ymax = mean_y + SE_y),
                        width = 5, size = 0.75) +
          geom_errorbarh(aes(xmin = mean_x-SE_x, xmax = mean_x + SE_x),
                         height = 5, size = 0.75) +
          guides(fill = guide_legend(override.aes = list(size = 8, color = c('blue','red'), order = 1),
                                     shape = guide_legend(override.aes = list(size = c(8,8),fill ='black')), order = 2)) +
          ggtitle('Position C') +
          theme(text = element_text(family = 'sans'),
                plot.title = element_text(size = 24, hjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.key = element_rect(fill = 'white'),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 14),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_rect(fill = 'white', color = 'black'),
                panel.border = element_blank(),
                panel.grid = element_blank()))
dev.off()


##----Position D----
pdf(paste0(output_filepath,'/Cart_data_D.pdf'), width = 10, height = 8)
print(ggplot(data = subset(cart_summary, Location == 'D'), aes(x = mean_x, y = mean_y)) +
          geom_point(aes(fill = Treatment, shape = Morphant, size = Morphant), color = 'black', stroke = 0.05) +
          geom_hline(yintercept = 0) +
          geom_vline(xintercept = 0) +
          scale_y_continuous(name = 'y-axis', expand = c(0,0), limits = c(-165,0), position = 'top') +
          scale_x_continuous(name = 'x-axis', expand = c(0,0), limits = c(0,255), position = 'left') +
          scale_shape_manual(values = c(22,23)) +
          scale_size_manual(values = c(9,9)) +
          scale_fill_manual(values = c('blue','red')) +
          geom_errorbar(aes(ymin = mean_y-SE_y, ymax = mean_y + SE_y),
                        width = 5, size = 0.75) +
          geom_errorbarh(aes(xmin = mean_x-SE_x, xmax = mean_x + SE_x),
                         height = 5, size = 0.75) +
          guides(fill = guide_legend(override.aes = list(size = 8, color = c('blue','red'), order = 1),
                                     shape = guide_legend(override.aes = list(size = c(8,8),fill ='black')), order = 2)) +
          ggtitle('Position D') +
          theme(text = element_text(family = 'sans'),
                plot.title = element_text(size = 24, hjust = 0.5),
                legend.title = element_text(size = 14),
                legend.text = element_text(size = 12),
                legend.key = element_rect(fill = 'white'),
                axis.title = element_text(size = 18),
                axis.text = element_text(size = 14),
                axis.ticks = element_line(color = 'black'),
                panel.background = element_rect(fill = 'white', color = 'black'),
                panel.border = element_blank(),
                panel.grid = element_blank()))
dev.off()


##----5) Write session information----
writeLines(capture.output({
    cat('###########################################################################\n##  Script and Data Analysis by Cheryl Dunham                            ##\n##                                                                       ##\n##  Budding Bioinformatician                                             ##\n##  Tanguay Laboratory | Sinnhuber Aquatic Research Laboratory (SARL)    ##\n##  28645 E Hwy 34                                                       ##\n##  Corvallis, OR 97333                                                  ##\n##  Oregon State University                                              ##\n##                                                                       ##\n##  (541) 737-6500                                                       ##\n##  dunhamcg@gmail.com                                                   ##\n###########################################################################\n\n')
    cat('Input filepath:', input_filepath, '\n\n')
    cat('Output filepath:', output_filepath, '\n\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n')
    print(sessionInfo())
    cat('\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nR Citation Information\n')
    print(citation(package = 'base'))
    cat('\n_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nRStudio Citation Information\n\n')
    print(RStudio.Version())
    cat('_____________________________________________________________________________________\n_____________________________________________________________________________________\n\n\nAttached Package Citation Information\n\n')
    print(lapply(all_packages, citation))
    
}), paste0(output_filepath,'/Morphometric_analysis_sessionInfo_',Sys.Date(),'.txt'))
