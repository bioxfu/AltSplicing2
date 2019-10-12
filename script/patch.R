x <- read.table('bam1/Result/Diff_compare1_intron_PI.txt', head=T, stringsAsFactors = F)
gene <- read.table('/cluster/home/xfu/Gmatic7/genome/human/GRCh38_gene2name.tsv', stringsAsFactors = F)
colnames(gene) <- c('Gene_id', 'GeneName')

y <- x[(x$FDR_rMATS<0.05 | x$FDR_DEXSeq<0.05) & abs(x$Diff_PI_Junction)>0.1, ]
y <- merge(gene, y, by.x=1, by.y=2, all.y=T)

inc_1 <- do.call(rbind, strsplit(y$Inclusion_counts_SAMPLE1, ','))
colnames(inc_1) <- paste0('Inclusion_counts_SAMPLE1_rep', 1:ncol(inc_1))
skp_1 <- do.call(rbind, strsplit(y$Skipping_counts_SAMPLE1, ','))
colnames(skp_1) <- paste0('Skipping_counts_SAMPLE1_rep', 1:ncol(skp_1))
inc_2 <- do.call(rbind, strsplit(y$Inclusion_counts_SAMPLE2, ','))
colnames(inc_2) <- paste0('Inclusion_counts_SAMPLE2_rep', 1:ncol(inc_2))
skp_2 <- do.call(rbind, strsplit(y$Skipping_counts_SAMPLE2, ','))
colnames(skp_2) <- paste0('Skipping_counts_SAMPLE2_rep', 1:ncol(skp_2))

y <- cbind(y[,1:9], inc_1, skp_1, inc_2, skp_2, y[,14:ncol(y)])

write.table(y, 'bam1/Result/Diff_compare1_intron_PI_FDR0.05_DiffPIJunc0.1.txt', sep='\t', row.names=F, quote=F)
