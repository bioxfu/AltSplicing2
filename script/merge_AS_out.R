library(RColorBrewer)
cols <- brewer.pal(3, 'Set2')

#AS_lst <- c('SE', 'RI', 'MXE', 'A3SS', 'A5SS')
AS_lst <- c('SE')
#VS <- c('loc_EV_vs_loc_C3_mut', 'loc_EV_vs_loc_TYLCV', 'sys_EV_vs_sys_TYLCV', 'RL1_vs_RL5', 'RL1_vs_RL7')

VS <- commandArgs(T)

merge_tables_significant <- function(P) {
  for (AS in AS_lst) {
    dfms <- list()
    for (i in 1:length(VS)) {
      dfm <- read.table(paste0('rMATS_out_reformat/', VS[i], '.', AS), header = T)
      dfm$significant <- 0 
      dfm$significant[(dfm[, P] < 0.05 & abs(dfm$IncLevelDifference) > 0.1 & dfm$IncLevelDifference > 0)] <- 1
      dfm$significant[(dfm[, P] < 0.05 & abs(dfm$IncLevelDifference) > 0.1 & dfm$IncLevelDifference < 0)] <- -1
      colnames(dfm)[-c(1,2,3)] <- paste0(VS[i], '.', colnames(dfm)[-c(1,2,3)])
      dfms[[i]] <- dfm
    }
    
    merge_dfm <- dfms[[1]]
    for (i in 2:length(VS)) {
      merge_dfm <- merge(merge_dfm, dfms[[i]][c(-2,-3)], by.x = 1, by.y = 1)
    }
    
    merge_dfm_sig <- merge_dfm[rowSums(merge_dfm[grep('significant', colnames(merge_dfm))]) != 0, ]
    write.table(merge_dfm_sig, paste0('tables/merge_', AS, '_', P, '0.05_deltaPSI0.1'), row.names=F, quote=F, sep='\t')
    
    stat <- sapply(merge_dfm_sig[, grep('significant', colnames(merge_dfm_sig))], table)
    stat <- stat[c('-1', '1'),]
    rownames(stat) <- c('exon_inclusion_up', 'exon_inclusion_down')
    colnames(stat) <- sub('.significant', '', colnames(stat))

    pdf(paste0('tables/merge_', AS, '_', P, '0.05_deltaPSI0.1.pdf'))
    par(mar=c(10,5,2,2))
    bp <- barplot(stat, beside = T, las=2, col = cols[1:2], ylim=c(0, max(stat)*1.2), border = 'white', ylab = 'Number of events')
    text(bp[1, ], stat[1, ], stat[1, ], pos = 3)
    text(bp[2, ], stat[2, ], stat[2, ], pos = 3)
    legend('topleft', rownames(stat), fill = cols[1:2], bty='n', border = NA)
    dev.off()
  }
}

merge_tables_significant('PValue')
merge_tables_significant('FDR')

