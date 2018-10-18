library(PAIRADISE)

input_files <- dir(path='rMATS_out', pattern='*/*.MATS.JCEC.txt', recursive=T, full.names=T)
output_files <- sub('rMATS_out', 'pairadise', input_files)
output_files <- sub('.txt', '.paired.compared.txt', output_files)

for (i in 1:length(input_files)) {
  print(input_files[i])
  AS <- read.table(input_files[i], header = T)
  AS_sub <- AS[c('ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncFormLen', 'SkipFormLen')]
  results <- pairadise(AS_sub, equal.variance = FALSE, numCluster = 8, sig.level = 0.999)
  dfm <- results$sig.results.FDR
  colnames(dfm) <- paste0(colnames(dfm), '(Paired-test)')
  output <- merge(AS, dfm, by.x=1, by.y=1)
  write.table(output, output_files[i], sep='\t', row.names = F, quote = F)
}
