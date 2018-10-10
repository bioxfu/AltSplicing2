AS_lst <- c('SE', 'RI', 'MXE', 'A3SS', 'A5SS')
#VS <- c('mdg10nramp1_Mn_vs_nramp1_Mn', 'mdg10nramp1_vs_mdg10nramp1_Mn', 'mdg10nramp1_vs_nramp1', 'nramp1_vs_nramp1_Mn')

VS <- commandArgs(T)

for (AS in AS_lst) {

  dfms <- list()
  for (i in 1:length(VS)) {
    dfm <- read.table(paste0('rMATS_out_reformat/', VS[i], '.', AS), header = T)
    dfm$significant <- 0 
    dfm$significant[(dfm$PValue < 0.05 & abs(dfm$IncLevelDifference) > 0.1 & dfm$IncLevelDifference > 0)] <- 1
    dfm$significant[(dfm$PValue < 0.05 & abs(dfm$IncLevelDifference) > 0.1 & dfm$IncLevelDifference < 0)] <- -1
    colnames(dfm)[-c(1,2)] <- paste0(VS[i], '.', colnames(dfm)[-c(1,2)])
    dfms[[i]] <- dfm
  }
  
  merge_dfm <- dfms[[1]]
  for (i in 2:length(VS)) {
    merge_dfm <- merge(merge_dfm, dfms[[i]][-2], by.x = 1, by.y = 1)
  }
  
  merge_dfm_sig <- merge_dfm[rowSums(merge_dfm[grep('significant', colnames(merge_dfm))]) != 0, ]
  
  write.table(merge_dfm_sig, paste0('tables/merge_', AS, '_PValue0.05_deltaPSI0.1'), row.names=F, quote=F, sep='\t')
              
}

