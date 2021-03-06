library(magrittr)
library(GeneOverlap)
library(gplots)
library(plotrix)
library(topGO)
library(GO.db)

topGO <- function(myGenes, category='BP', p_cutoff=0.05, gomap, geneid){
  geneID2GO <- readMappings(file = gomap)
  geneNames <- read.table(geneid, stringsAsFactors=F)$V1
  geneList <- factor(as.integer(geneNames %in% myGenes))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata", ontology=category, allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
  allRes <- GenTable(GOdata,pvalue=resultFisher,topNodes=100)
  allRes$pvalue[grep('<',allRes$pvalue)] <- "1e-30"
  allRes$pvalue <- as.numeric(allRes$pvalue)
  allRes <- allRes[order(allRes$pvalue,decreasing=F),]
  allRes$catagory <- category
  allRes <- allRes[allRes$pvalue < p_cutoff, ]
  return(allRes)
}

argv <- commandArgs(T)
fdr_cutoff <- as.numeric(argv[1])
psi_cutoff <- as.numeric(argv[2])
Sample1_path <- argv[3]
Sample2_path <- argv[4]
Sample3_path <- argv[5]
Sample1_name <- argv[6]
Sample2_name <- argv[7]
Sample3_name <- argv[8]
rpkm_path <- argv[9]
rpkm_cutoff <- as.numeric(argv[10])
# rpkm_path <- '/data1/HJY/Project/20181121JY/RNA-Seq/table_8samples/RPKM_table_FDR0.05_FC1.5_all.tsv'
# psi_cutoff <- 0.1
# fdr_cutoff <- 0.01
# rpkm_cutoff <- 2
# Sample1_path <- '/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_MeCP2_KO.14/'
# Sample2_path <- '/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KO/'
# Sample3_path <- '/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KI/'
# Sample1_name <- 'MeCP2_KO'
# Sample2_name <- 'RBFOX2_KO'
# Sample3_name <- 'RBFOX2_KI'
rpkm <- read.table(rpkm_path, header = T, row.names = 1)
rpkm <- rpkm[1:(grep('_vs_', colnames(rpkm))[1]-1)]
selected_gene <- rownames(rpkm)[rowSums(rpkm)/ncol(rpkm) >= rpkm_cutoff]

Sample1_all <- NULL
Sample2_all <- NULL
Sample3_all <- NULL

for (AS in c('SE', 'RI', 'MXE', 'A3SS', 'A5SS')) {
  Sample1 <- read.table(paste0(Sample1_path, '/', AS, '.MATS.JCEC.txt'), header = T, sep = '\t')
  Sample2 <- read.table(paste0(Sample2_path, '/', AS, '.MATS.JCEC.txt'), header = T, sep = '\t')
  Sample3 <- read.table(paste0(Sample3_path, '/', AS, '.MATS.JCEC.txt'), header = T, sep = '\t')
  Sample1 <- Sample1[Sample1$GeneID %in% selected_gene, ]
  Sample2 <- Sample2[Sample2$GeneID %in% selected_gene, ]
  Sample3 <- Sample3[Sample3$GeneID %in% selected_gene, ]
  if (AS == 'SE') {
    Sample1$event <- paste('SE', Sample1$chr, Sample1$strand, Sample1$exonStart_0base, Sample1$exonEnd, Sample1$upstreamES, Sample1$upstreamEE, Sample1$downstreamES, Sample1$downstreamEE, Sample1$geneSymbol, Sample1$GeneID, sep='|')
    Sample2$event <- paste('SE', Sample2$chr, Sample2$strand, Sample2$exonStart_0base, Sample2$exonEnd, Sample2$upstreamES, Sample2$upstreamEE, Sample2$downstreamES, Sample2$downstreamEE, Sample2$geneSymbol, Sample2$GeneID, sep='|')
    Sample3$event <- paste('SE', Sample3$chr, Sample3$strand, Sample3$exonStart_0base, Sample3$exonEnd, Sample3$upstreamES, Sample3$upstreamEE, Sample3$downstreamES, Sample3$downstreamEE, Sample3$geneSymbol, Sample3$GeneID, sep='|')
  }
  if (AS == 'RI') {
    Sample1$event <- paste('RI', Sample1$chr, Sample1$strand, Sample1$riExonStart_0base, Sample1$riExonEnd, Sample1$upstreamES, Sample1$upstreamEE, Sample1$downstreamES, Sample1$downstreamEE, Sample1$geneSymbol, Sample1$GeneID, sep='|')
    Sample2$event <- paste('RI', Sample2$chr, Sample2$strand, Sample2$riExonStart_0base, Sample2$riExonEnd, Sample2$upstreamES, Sample2$upstreamEE, Sample2$downstreamES, Sample2$downstreamEE, Sample2$geneSymbol, Sample2$GeneID, sep='|')
    Sample3$event <- paste('RI', Sample3$chr, Sample3$strand, Sample3$riExonStart_0base, Sample3$riExonEnd, Sample3$upstreamES, Sample3$upstreamEE, Sample3$downstreamES, Sample3$downstreamEE, Sample3$geneSymbol, Sample3$GeneID, sep='|')
  }
  if (AS == 'MXE') {
    Sample1$event <- paste('MXE', Sample1$chr, Sample1$strand, Sample1$X1stExonStart_0base, Sample1$X1stExonEnd, Sample1$X2ndExonStart_0base, Sample1$X2ndExonEnd, Sample1$upstreamES, Sample1$upstreamEE, Sample1$downstreamES, Sample1$downstreamEE, Sample1$geneSymbol, Sample1$GeneID, sep='|')
    Sample2$event <- paste('MXE', Sample2$chr, Sample2$strand, Sample2$X1stExonStart_0base, Sample2$X1stExonEnd, Sample2$X2ndExonStart_0base, Sample2$X2ndExonEnd, Sample2$upstreamES, Sample2$upstreamEE, Sample2$downstreamES, Sample2$downstreamEE, Sample2$geneSymbol, Sample2$GeneID, sep='|')
    Sample3$event <- paste('MXE', Sample3$chr, Sample3$strand, Sample3$X1stExonStart_0base, Sample3$X1stExonEnd, Sample3$X2ndExonStart_0base, Sample3$X2ndExonEnd, Sample3$upstreamES, Sample3$upstreamEE, Sample3$downstreamES, Sample3$downstreamEE, Sample3$geneSymbol, Sample3$GeneID, sep='|')
  }
  if (AS == 'A3SS') {
    Sample1$event <- paste('A3SS', Sample1$chr, Sample1$strand, Sample1$longExonStart_0base, Sample1$longExonEnd, Sample1$shortES, Sample1$shortEE, Sample1$flankingES, Sample1$flankingEE, Sample1$geneSymbol, Sample1$GeneID, sep='|')
    Sample2$event <- paste('A3SS', Sample2$chr, Sample2$strand, Sample2$longExonStart_0base, Sample2$longExonEnd, Sample2$shortES, Sample2$shortEE, Sample2$flankingES, Sample2$flankingEE, Sample2$geneSymbol, Sample2$GeneID, sep='|')
    Sample3$event <- paste('A3SS', Sample3$chr, Sample3$strand, Sample3$longExonStart_0base, Sample3$longExonEnd, Sample3$shortES, Sample3$shortEE, Sample3$flankingES, Sample3$flankingEE, Sample3$geneSymbol, Sample3$GeneID, sep='|')
  }
  if (AS == 'A5SS') {
    Sample1$event <- paste('A5SS', Sample1$chr, Sample1$strand, Sample1$longExonStart_0base, Sample1$longExonEnd, Sample1$shortES, Sample1$shortEE, Sample1$flankingES, Sample1$flankingEE, Sample1$geneSymbol, Sample1$GeneID, sep='|')
    Sample2$event <- paste('A5SS', Sample2$chr, Sample2$strand, Sample2$longExonStart_0base, Sample2$longExonEnd, Sample2$shortES, Sample2$shortEE, Sample2$flankingES, Sample2$flankingEE, Sample2$geneSymbol, Sample2$GeneID, sep='|')
    Sample3$event <- paste('A5SS', Sample3$chr, Sample3$strand, Sample3$longExonStart_0base, Sample3$longExonEnd, Sample3$shortES, Sample3$shortEE, Sample3$flankingES, Sample3$flankingEE, Sample3$geneSymbol, Sample3$GeneID, sep='|')
  }
  
  Sample1_sig <- abs(Sample1$IncLevelDifference) > psi_cutoff & Sample1$FDR < fdr_cutoff
  Sample2_sig <- abs(Sample2$IncLevelDifference) > psi_cutoff & Sample2$FDR < fdr_cutoff
  Sample3_sig <- abs(Sample3$IncLevelDifference) > psi_cutoff & Sample3$FDR < fdr_cutoff
  
  Sample1_all <- rbind(Sample1_all, Sample1[Sample1_sig, c('event', 'IncLevelDifference')])
  Sample2_all <- rbind(Sample2_all, Sample2[Sample2_sig, c('event', 'IncLevelDifference')])
  Sample3_all <- rbind(Sample3_all, Sample3[Sample3_sig, c('event', 'IncLevelDifference')])
  
}

gene_lst <- list(Sample1=Sample1_all[,'event'], Sample2=Sample2_all[,'event'], Sample3=Sample3_all[,'event'])
names(gene_lst) <- c(Sample1_name, Sample2_name, Sample3_name)
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

pdf(paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '.pdf'), wid=10)
layout(matrix(c(1,2),nrow=1), wid=c(1,1))
par(mar=c(2,0,2,0))
par(xpd=T)
v <- venn(gene_lst)
plot(0, type='n', xaxt='n', yaxt='n', bty='n')
addtable2plot(0.65, 0, mat, display.rownames=F, display.colnames=T, bty='n')
dev.off()

pdf(paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_', Sample1_name, '_vs_', Sample2_name, '.pdf'))
v12 <- venn(gene_lst[c(Sample1_name, Sample2_name)])
write.table(attr(v12, 'intersections')[[paste0(Sample1_name, ':', Sample2_name)]], paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_', Sample1_name, '_vs_', Sample2_name, '_olp.txt'), row.names = F, col.names = F, quote = F)
dev.off()

pdf(paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_', Sample1_name, '_vs_', Sample3_name, '.pdf'))
v13 <- venn(gene_lst[c(Sample1_name, Sample3_name)])
write.table(attr(v13, 'intersections')[[paste0(Sample1_name, ':', Sample3_name)]], paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_', Sample1_name, '_vs_', Sample3_name, '_olp.txt'), row.names = F, col.names = F, quote = F)
dev.off()

pdf(paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_', Sample2_name, '_vs_', Sample3_name, '.pdf'))
v23 <- venn(gene_lst[c(Sample2_name, Sample3_name)])
write.table(attr(v23, 'intersections')[[paste0(Sample2_name, ':', Sample3_name)]], paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_', Sample2_name, '_vs_', Sample3_name, '_olp.txt'), row.names = F, col.names = F, quote = F)
dev.off()

Sample12_all <- merge(Sample1_all, Sample2_all, by.x = 1, by.y = 1)
Sample13_all <- merge(Sample1_all, Sample3_all, by.x = 1, by.y = 1)
Sample23_all <- merge(Sample2_all, Sample3_all, by.x = 1, by.y = 1)

intersections <- attr(v, 'intersections')
s123 <- intersections[[paste0(Sample1_name, ':', Sample2_name, ':', Sample3_name)]]
s12 <- intersections[[paste0(Sample1_name, ':', Sample2_name)]]
s23 <- intersections[[paste0(Sample2_name, ':', Sample3_name)]]
s13 <- intersections[[paste0(Sample1_name, ':', Sample3_name)]]

pdf(paste0('figures/AS_overlap_FDR', fdr_cutoff, '_PSI', psi_cutoff, '_RPKM', rpkm_cutoff, '_correlation', '.pdf'), wid=10)
par(mfrow=c(2,3))
# Sample1 overlap Sample2
x <- Sample12_all[Sample12_all$event %in% c(s12, s123), 2]
y <- Sample12_all[Sample12_all$event %in% c(s12, s123), 3]
plot(x, y, xlab=paste0('IncLevelDifference of ', Sample1_name), ylab=paste0('IncLevelDifference of ', Sample2_name), main = paste0(Sample1_name, ' overlap with ', Sample2_name, ' (', length(x), ')'))
text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))

# Sample1 overlap Sample3
x <- Sample13_all[Sample13_all$event %in% c(s13, s123), 2]
y <- Sample13_all[Sample13_all$event %in% c(s13, s123), 3]
plot(x, y, xlab=paste0('IncLevelDifference of ', Sample1_name), ylab=paste0('IncLevelDifference of ', Sample3_name), main = paste0(Sample1_name, ' overlap with ', Sample3_name, ' (', length(x), ')'))
text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))

# Sample2 overlap Sample3
x <- Sample23_all[Sample23_all$event %in% c(s23, s123), 2]
y <- Sample23_all[Sample23_all$event %in% c(s23, s123), 3]
plot(x, y, xlab=paste0('IncLevelDifference of ', Sample2_name), ylab=paste0('IncLevelDifference of ', Sample3_name), main = paste0(Sample2_name, ' overlap with ', Sample3_name, ' (', length(x), ')'))
text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))

# Sample1 overlap Sample2 and Sample3
x <- Sample12_all[Sample12_all$event %in% s123, 2]
y <- Sample12_all[Sample12_all$event %in% s123, 3]
plot(x, y, xlab=paste0('IncLevelDifference of ', Sample1_name), ylab=paste0('IncLevelDifference of ', Sample2_name), main = paste0('three samples overlap (', length(x), ')'))
text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))

x <- Sample13_all[Sample13_all$event %in% s123, 2]
y <- Sample13_all[Sample13_all$event %in% s123, 3]
plot(x, y, xlab=paste0('IncLevelDifference of ', Sample1_name), ylab=paste0('IncLevelDifference of ', Sample3_name), main = paste0('three samples overlap (', length(x), ')'))
text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))

x <- Sample23_all[Sample23_all$event %in% s123, 2]
y <- Sample23_all[Sample23_all$event %in% s123, 3]
plot(x, y, xlab=paste0('IncLevelDifference of ', Sample2_name), ylab=paste0('IncLevelDifference of ', Sample3_name), main = paste0('three samples overlap (', length(x), ')'))
text(-0.8, 0.8, substitute(paste(italic(r), '=', x), list(x=round(cor(x,y),2))))

dev.off()


