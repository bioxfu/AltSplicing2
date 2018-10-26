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
psi_cutoff <- argv[1]
KO_path <- argv[2]
MT_path <- argv[3]
# psi_cutoff <- 0.1
# KO_path <- '/data1/HJY/Project/Mecp2/KO/AltSplicing2_124/pairadise/WT_vs_KO/'
# MT_path <- '/data1/HJY/Project/Mecp2/T158M/AltSplicing2_134/pairadise/WT_vs_T158/'
KO_all <- NULL
MT_all <- NULL

for (AS in c('SE', 'RI', 'MXE', 'A3SS', 'A5SS')) {
  KO <- read.table(paste0(KO_path, '/', AS, '.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t')
  MT <- read.table(paste0(MT_path, '/', AS, '.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t')
  if (AS == 'SE') {
    KO$event <- paste('SE', KO$chr, KO$strand, KO$exonStart_0base, KO$exonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, KO$geneSymbol, KO$GeneID, sep='|')
    MT$event <- paste('SE', MT$chr, MT$strand, MT$exonStart_0base, MT$exonEnd, MT$upstreamES, MT$upstreamEE, MT$downstreamES, MT$downstreamEE, MT$geneSymbol, MT$GeneID, sep='|')
  }
  if (AS == 'RI') {
    KO$event <- paste('RI', KO$chr, KO$strand, KO$riExonStart_0base, KO$riExonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, KO$geneSymbol, KO$GeneID, sep='|')
    MT$event <- paste('RI', MT$chr, MT$strand, MT$riExonStart_0base, MT$riExonEnd, MT$upstreamES, MT$upstreamEE, MT$downstreamES, MT$downstreamEE, MT$geneSymbol, MT$GeneID, sep='|')
  }
  if (AS == 'MXE') {
    KO$event <- paste('MXE', KO$chr, KO$strand, KO$X1stExonStart_0base, KO$X1stExonEnd, KO$X2ndExonStart_0base, KO$X2ndExonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, KO$geneSymbol, KO$GeneID, sep='|')
    MT$event <- paste('MXE', MT$chr, MT$strand, MT$X1stExonStart_0base, MT$X1stExonEnd, MT$X2ndExonStart_0base, MT$X2ndExonEnd, MT$upstreamES, MT$upstreamEE, MT$downstreamES, MT$downstreamEE, MT$geneSymbol, MT$GeneID, sep='|')
  }
  if (AS == 'A3SS') {
    KO$event <- paste('A3SS', KO$chr, KO$strand, KO$longExonStart_0base, KO$longExonEnd, KO$shortES, KO$shortEE, KO$flankingES, KO$flankingEE, KO$geneSymbol, KO$GeneID, sep='|')
    MT$event <- paste('A3SS', MT$chr, MT$strand, MT$longExonStart_0base, MT$longExonEnd, MT$shortES, MT$shortEE, MT$flankingES, MT$flankingEE, MT$geneSymbol, MT$GeneID, sep='|')
  }
  if (AS == 'A5SS') {
    KO$event <- paste('A5SS', KO$chr, KO$strand, KO$longExonStart_0base, KO$longExonEnd, KO$shortES, KO$shortEE, KO$flankingES, KO$flankingEE, KO$geneSymbol, KO$GeneID, sep='|')
    MT$event <- paste('A5SS', MT$chr, MT$strand, MT$longExonStart_0base, MT$longExonEnd, MT$shortES, MT$shortEE, MT$flankingES, MT$flankingEE, MT$geneSymbol, MT$GeneID, sep='|')
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
v <- venn(gene_lst)
plot(0, type='n', xaxt='n', yaxt='n', bty='n')
addtable2plot(0.45, 0, mat, display.rownames=F, display.colnames=T, bty='n')
dev.off()

intersections <- attr(v, 'intersections')
olp_gene <- sub('\\..+', '', sub('.+\\|', '', intersections$`KO:Mutant`))
go <- topGO(olp_gene, gomap='script/mouse_gene2GO.map', geneid='script/mouse_genewithGO')
go$Term <- apply(go, 1, function(x){Term(GOTERM[[x[1]]])})
write.table(go, paste0('table/KO_T158M_AS_overlap_PSI', psi_cutoff, '.GO.tsv'), sep = '\t', quote = F, row.names = F)
write.table(intersections$`KO:Mutant`, paste0('table/KO_T158M_AS_overlap_PSI', psi_cutoff, '.events.tsv'), sep = '\t', quote = F, row.names = F, col.names = F)

