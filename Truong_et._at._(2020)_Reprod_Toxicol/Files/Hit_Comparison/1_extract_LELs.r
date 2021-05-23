rm(list=ls())

#-----------------------------------------------------------
# Noyes
#-----------------------------------------------------------
df <- read.csv('data_in/Mapping of Noyes LEL for 44 FRCs 2020MAR30.csv', skip = 1, stringsAsFactors = F)
names(df)[c(2, 27:31)] <- c('Chemical.ID', 'Background', 'Excitatory', 'Refractory', 'Light', 'Dark')

idx <- df[, 5] != ''
noyes <- list(
  'morph'    = df[idx, c(2,  5:26)],
  'behav24h' = df[idx, c(2, 27:29)],
  'behav5d'  = df[idx, c(2, 30:31)]
)

toLEL <- function(x) {
  x <- sapply(unlist(x, use.names = F), trimws)
  x[x == '-'] <- NA
  idx <- which(grepl('hyper', x, ignore.case = T))
  if (length(idx > 0)) {
    x[idx] <- sub('\\s?\\(hyper.*', '', x[idx], ignore.case = T)
  }
  idx <- which(grepl('hypo', x, ignore.case = T))
  if (length(idx > 0)) {
    x[idx] <- paste0('-', sub('\\s?\\(hypo.*', '', x[idx], ignore.case = T))
  }
  return(x)
}

for (assay in c('morph', 'behav24h', 'behav5d')) {
  df <- noyes[[assay]]
  for (i in 2:ncol(df)) {
    noyes[[assay]][, i] <- as.numeric(toLEL(df[, i]))
  }
}

#-----------------------------------------------------------
# Truong
#-----------------------------------------------------------
truong <- list(
  'morph'    = read.csv('data_in/LEL_morph.csv', stringsAsFactors = F),
  'behav24h' = read.csv('data_in/LEL_24h.csv', stringsAsFactors = F),
  'behav5d'  = read.csv('data_in/LEL_5d.csv', stringsAsFactors = F)
)

truong[['morph']] <- truong[['morph']][, 1:(ncol(truong[['morph']]) - 7)]
for (i in 2:ncol(truong[['morph']])) {
  truong[['morph']][, i] <- as.numeric(sub('\\s?uM', '', truong[['morph']][, i]))
}

getLEL <- function(df, intervals) {
  out <- NULL
  for (interval in intervals) {
    x <- df[df$Interval == interval, c(1, 3, 5)]
    x <- data.frame(x$Chemical.ID, sign(x$potencyLEL) * as.numeric(sub('\\s?uM', '', x$CONC)))
    names(x) <- c('Chemical.ID', interval)
    out <- if (is.null(out)) x else merge(out, x)
  }
  return(out)
}

for (assay in c('behav24h', 'behav5d')) {
  intervals <- if (assay == 'behav24h') c('Background', 'Excitatory', 'Refractory') else c('Light', 'Dark')
  truong[[assay]] <- getLEL(truong[[assay]], intervals)
}

#-----------------------------------------------------------
# Write overlap
#-----------------------------------------------------------
chems <- sort(noyes[['morph']]$Chemical.ID[noyes[['morph']]$Chemical.ID %in% truong[['morph']]$Chemical.ID])

for (assay in c('morph', 'behav24h', 'behav5d')) {
  write.csv(noyes[[assay]][match(chems, noyes[[assay]][, 1]), ], paste0('data_out/noyes_', assay, '.csv'), row.names = F)
  write.csv(truong[[assay]][match(chems, truong[[assay]][, 1]), ], paste0('data_out/truong_', assay, '.csv'), row.names = F)
}

