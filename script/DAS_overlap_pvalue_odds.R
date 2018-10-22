library(magrittr)
library(GeneOverlap)
library(gplots)
library(plotrix)

argv <- commandArgs(T)
psi_cutoff <- argv[1]
KO_path <- argv[2]
MT_path <- argv[3]
KO_all <- NULL
MT_all <- NULL

for (AS in c('SE', 'RI', 'MXE', 'A3SS', 'A5SS')) {
  KO <- read.table(paste0(KO_path, '/', AS, '.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t')
  MT <- read.table(paste0(MT_path, '/', AS, '.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t')
  if (AS == 'SE') {
    KO$event <- paste(KO$chr, KO$strand, KO$exonStart_0base, KO$exonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, sep='|')
    MT$event <- paste(MT$chr, MT$strand, MT$exonStart_0base, MT$exonEnd, MT$upstreamES, MT$upstreamEE, MT$downstreamES, MT$downstreamEE, sep='|')
  }
  if (AS == 'RI') {
    KO$event <- paste(KO$chr, KO$strand, KO$riExonStart_0base, KO$riExonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, sep='|')
    MT$event <- paste(MT$chr, MT$strand, MT$riExonStart_0base, MT$riExonEnd, MT$upstreamES, MT$upstreamEE, MT$downstreamES, MT$downstreamEE, sep='|')
  }
  if (AS == 'MXE') {
    KO$event <- paste(KO$chr, KO$strand, KO$X1stExonStart_0base, KO$X1stExonEnd, KO$X2ndExonStart_0base, KO$X2ndExonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, sep='|')
    MT$event <- paste(MT$chr, MT$strand, MT$X1stExonStart_0base, MT$X1stExonEnd, MT$X2ndExonStart_0base, MT$X2ndExonEnd, MT$upstreamES, MT$upstreamEE, MT$downstreamES, MT$downstreamEE, sep='|')
  }
  if (AS == 'A3SS') {
    KO$event <- paste(KO$chr, KO$strand, KO$longExonStart_0base, KO$longExonEnd, KO$shortES, KO$shortEE, KO$flankingES, KO$flankingEE, sep='|')
    MT$event <- paste(MT$chr, MT$strand, MT$longExonStart_0base, MT$longExonEnd, MT$shortES, MT$shortEE, MT$flankingES, MT$flankingEE, sep='|')
  }
  if (AS == 'A5SS') {
    KO$event <- paste(KO$chr, KO$strand, KO$longExonStart_0base, KO$longExonEnd, KO$shortES, KO$shortEE, KO$flankingES, KO$flankingEE, sep='|')
    MT$event <- paste(MT$chr, MT$strand, MT$longExonStart_0base, MT$longExonEnd, MT$shortES, MT$shortEE, MT$flankingES, MT$flankingEE, sep='|')
  }
  
  KO_sig <- abs(KO$IncLevelDifference) > psi_cutoff & KO$Raw.p.values.Paired.test. < 0.05
  MT_sig <- abs(MT$IncLevelDifference) > psi_cutoff & MT$Raw.p.values.Paired.test. < 0.05
  
  KO_all <- c(KO_all, KO[KO_sig, 'event'])
  MT_all <- c(MT_all, MT[MT_sig, 'event'])
  
}

gene_lst <- list(KO=KO_all, Mutant=MT_all)
genome <- read.table('all_possible_AS_events')

mat <- NULL
for (i in 1:(length(gene_lst)-1)) {
  for (j in (i+1):length(gene_lst)) {
    go.obj <- newGeneOverlap(gene_lst[[i]], gene_lst[[j]], genome.size=nrow(genome)) %>% testGeneOverlap()
    num <- getIntersection(go.obj) %>% length()
    p <- getPval(go.obj) %>% format(digits=2)
    odds <- getOddsRatio(go.obj) %>% format(digits=2)
    mat <- rbind(mat, c(names(gene_lst)[i], names(gene_lst)[j], num, p, odds))
  }
}
colnames(mat) <- c('List1', 'List2', 'Number', 'P-value', 'OddsRatio')

pdf(paste0('figure/KO_T158M_AS_overlap_PSI', psi_cutoff, '.pdf'), wid=10)
layout(matrix(c(1,2),nrow=1), wid=c(2,1))
par(mar=c(2,2,2,2))
par(xpd=T)
venn(gene_lst)
plot(0, type='n', xaxt='n', yaxt='n', bty='n')
addtable2plot(0.45, 0, mat, display.rownames=F, display.colnames=T, bty='n')
dev.off()

