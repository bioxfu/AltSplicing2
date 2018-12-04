

for (i in c('12', '13', '14', '23', '24', '34')) {
  pdf(paste0('figures/PSI_correlation_MeCP2_KO.', i, '.pdf'))
  dfm1 <- read.table(paste0('tables.', i, '/merge_SE_FDR0.05_deltaPSI0.1'), header = 1, sep = '\t')
  dfm2 <- read.table(paste0('tables.', i, '/merge_RI_FDR0.05_deltaPSI0.1'), header = 1, sep = '\t')
  dfm3 <- read.table(paste0('tables.', i, '/merge_MXE_FDR0.05_deltaPSI0.1'), header = 1, sep = '\t')
  dfm4 <- read.table(paste0('tables.', i, '/merge_A3SS_FDR0.05_deltaPSI0.1'), header = 1, sep = '\t')
  dfm5 <- read.table(paste0('tables.', i, '/merge_A5SS_FDR0.05_deltaPSI0.1'), header = 1, sep = '\t')

  dfm1 <- dfm1[dfm1[,paste0('WT_vs_MeCP2_KO.', i, '.significant')] != 0 & dfm1$WT_vs_RBFOX2_KO.significant != 0 & dfm1$WT_vs_RBFOX2_KI.significant != 0,]
  dfm2 <- dfm2[dfm2[,paste0('WT_vs_MeCP2_KO.', i, '.significant')] != 0 & dfm2$WT_vs_RBFOX2_KO.significant != 0 & dfm2$WT_vs_RBFOX2_KI.significant != 0,]
  dfm3 <- dfm3[dfm3[,paste0('WT_vs_MeCP2_KO.', i, '.significant')] != 0 & dfm3$WT_vs_RBFOX2_KO.significant != 0 & dfm3$WT_vs_RBFOX2_KI.significant != 0,]
  dfm4 <- dfm4[dfm4[,paste0('WT_vs_MeCP2_KO.', i, '.significant')] != 0 & dfm4$WT_vs_RBFOX2_KO.significant != 0 & dfm4$WT_vs_RBFOX2_KI.significant != 0,]
  dfm5 <- dfm5[dfm5[,paste0('WT_vs_MeCP2_KO.', i, '.significant')] != 0 & dfm5$WT_vs_RBFOX2_KO.significant != 0 & dfm5$WT_vs_RBFOX2_KI.significant != 0,]
  
  x <- c(dfm1[,paste0('WT_vs_MeCP2_KO.', i, '.IncLevelDifference')],
         dfm2[,paste0('WT_vs_MeCP2_KO.', i, '.IncLevelDifference')],
         dfm3[,paste0('WT_vs_MeCP2_KO.', i, '.IncLevelDifference')],
         dfm4[,paste0('WT_vs_MeCP2_KO.', i, '.IncLevelDifference')],
         dfm5[,paste0('WT_vs_MeCP2_KO.', i, '.IncLevelDifference')])
         
  y <- c(dfm1$WT_vs_RBFOX2_KO.IncLevelDifference,
         dfm2$WT_vs_RBFOX2_KO.IncLevelDifference,
         dfm3$WT_vs_RBFOX2_KO.IncLevelDifference,
         dfm4$WT_vs_RBFOX2_KO.IncLevelDifference,
         dfm5$WT_vs_RBFOX2_KO.IncLevelDifference)
         
  z <- c(dfm1$WT_vs_RBFOX2_KI.IncLevelDifference,
         dfm2$WT_vs_RBFOX2_KI.IncLevelDifference,
         dfm3$WT_vs_RBFOX2_KI.IncLevelDifference,
         dfm4$WT_vs_RBFOX2_KI.IncLevelDifference,
         dfm5$WT_vs_RBFOX2_KI.IncLevelDifference)
  
  par(mfrow=c(2,2))
  plot(x, y, xlab=paste0('MeCP2_KO.', i), ylab='RBFOX2_KO')
  text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))
  plot(x, z, xlab=paste0('MeCP2_KO.', i), ylab='RBFOX2_KI')
  text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,z),2))))
  plot(y, z, xlab='RBFOX2_KO', ylab='RBFOX2_KI')
  text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(y,z),2))))
  dev.off()
  
}

## 14 is the best
