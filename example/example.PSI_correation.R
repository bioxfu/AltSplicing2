

for (i in c('12', '13', '14', '23', '24', '34')) {
  pdf(paste0('figures/PSI_correlation_MeCP2_KO.', i, '.pdf'))
  dfm <- read.table(paste0('tables.', i, '/merge_SE_FDR0.05_deltaPSI0.1'), header = 1, sep = '\t')
  dfm <- dfm[dfm[,paste0('WT_vs_MeCP2_KO.', i, '.significant')] != 0 & dfm$WT_vs_RBFOX2_KO.significant != 0 & dfm$WT_vs_RBFOX2_KI.significant != 0,]
  x <- dfm[,paste0('WT_vs_MeCP2_KO.', i, '.IncLevelDifference')]
  y <- dfm$WT_vs_RBFOX2_KO.IncLevelDifference
  z <- dfm$WT_vs_RBFOX2_KI.IncLevelDifference
  
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
