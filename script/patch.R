x <- read.table('bam1/Result/Diff_compare1_intron_PI.txt', head=T)

y <- x[(x$FDR_rMATS<0.05 | x$FDR_DEXSeq<0.05) & abs(x$Diff_PI_Junction)>0.1, ]

write.table(y, 'bam1/Result/Diff_compare1_intron_PI_FDR0.05_DiffPIJunc0.1.txt', sep='\t', row.names=F, quote=F)
