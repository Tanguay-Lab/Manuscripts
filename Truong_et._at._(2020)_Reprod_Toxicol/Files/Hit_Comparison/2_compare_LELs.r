rm(list=ls())

save <- function(df, margin, fn) {
  if (margin == 1) {
    x <- data.frame(df1[, 1], apply(df, margin, function(x) round(100 * sum(x) / length(x), 2)))
    names(x) <- c('Chemical.ID', '% concordance')
    write.csv(x, fn, row.names = F)
  } else {
    df <- data.frame(names(df1)[-1], apply(hits, 2, function(x) round(100 * sum(x) / length(x), 2)))
    names(df) <- c('Endpoint', '% concordance')
    write.csv(df, fn, row.names = F)
  }
}

#-----------------------------------------------------------
# Morph
#-----------------------------------------------------------
df1 <- read.csv('data_out/noyes_morph.csv', stringsAsFactors = F)
df2 <- read.csv('data_out/truong_morph.csv', stringsAsFactors = F)

# Hits
hits <- !is.na(df1[, -1]) == !is.na(df2[, -1])
save(hits, 1, 'data_out/morph_chem_hits.csv')
save(hits, 2, 'data_out/morph_endpoint_hits.csv')
cat('Overall morph hit concorance = ', round(100 * sum(hits) / (nrow(hits) * ncol(hits)), 2), "\n")

# Potency - check if roughly same order of magnitude
pots <- round(log10(df1[, -1] / df2[, -1]))
pots <- is.na(pots) | pots == 0
pots[!hits] <- F
save(pots, 1, 'data_out/morph_chem_potency_hits.csv')
save(pots, 2, 'data_out/morph_endpoint_potency_hits.csv')
cat('Overall morph potency hit concorance = ', round(100 * sum(pots) / (nrow(pots) * ncol(pots)), 2), "\n")

#-----------------------------------------------------------
# Behav24h
#-----------------------------------------------------------
df1 <- read.csv('data_out/noyes_behav24h.csv', stringsAsFactors = F)
df2 <- read.csv('data_out/truong_behav24h.csv', stringsAsFactors = F)

# Hits
hits <- !is.na(df1[, -1]) == !is.na(df2[, -1])
save(hits, 1, 'data_out/behav24h_chem_hits.csv')
save(hits, 2, 'data_out/behav24h_endpoint_hits.csv')
cat('Overall behav24h hit concorance = ', round(100 * sum(hits) / (nrow(hits) * ncol(hits)), 2), "\n")

# Potency - check if roughly same order of magnitude and same direction
suppressWarnings(pots <- round(log10(df1[, -1] / df2[, -1])))
pots <- (is.na(pots) | pots == 0) & !sapply(pots, is.nan)
pots[!hits] <- F
save(pots, 1, 'data_out/behav24h_chem_potency_hits.csv')
save(pots, 2, 'data_out/behav24h_endpoint_potency_hits.csv')
cat('Overall behav24h potency hit concorance = ', round(100 * sum(pots) / (nrow(pots) * ncol(pots)), 2), "\n")

#-----------------------------------------------------------
# Behav5d
#-----------------------------------------------------------
df1 <- read.csv('data_out/noyes_behav5d.csv', stringsAsFactors = F)
df2 <- read.csv('data_out/truong_behav5d.csv', stringsAsFactors = F)

# Hits
hits <- !is.na(df1[, -1]) == !is.na(df2[, -1])
save(hits, 1, 'data_out/behav5d_chem_hits.csv')
save(hits, 2, 'data_out/behav5d_endpoint_hits.csv')
cat('Overall behav5d hit concorance = ', round(100 * sum(hits) / (nrow(hits) * ncol(hits)), 2), "\n")

# Potency - check if roughly same order of magnitude and same direction
suppressWarnings(pots <- round(log10(df1[, -1] / df2[, -1])))
pots <- (is.na(pots) | pots == 0) & !sapply(pots, is.nan)
pots[!hits] <- F
save(pots, 1, 'data_out/behav5d_chem_potency_hits.csv')
save(pots, 2, 'data_out/behav5d_endpoint_potency_hits.csv')
cat('Overall behav5d potency hit concorance = ', round(100 * sum(pots) / (nrow(pots) * ncol(pots)), 2), "\n")
