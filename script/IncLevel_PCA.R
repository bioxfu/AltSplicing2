library(RColorBrewer)

vs <- commandArgs(T)
#vs <- 'WT_vs_T158'
groups <- unlist(strsplit(vs, '_vs_'))

A3SS <- read.table(paste0('pairadise/', vs, '/A3SS.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t', stringsAsFactors = F)
A5SS <- read.table(paste0('pairadise/', vs, '/A5SS.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t', stringsAsFactors = F)
MXE <- read.table(paste0('pairadise/', vs, '/MXE.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t', stringsAsFactors = F)
RI <- read.table(paste0('pairadise/', vs, '/RI.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t', stringsAsFactors = F)
SE <- read.table(paste0('pairadise/', vs, '/SE.MATS.JCEC.paired.compared.txt'), header = T, sep = '\t', stringsAsFactors = F)

g1_str <- c(A3SS$IncLevel1, A5SS$IncLevel1, MXE$IncLevel1, RI$IncLevel1, SE$IncLevel1)
g2_str <- c(A3SS$IncLevel2, A5SS$IncLevel2, MXE$IncLevel2, RI$IncLevel2, SE$IncLevel2)

g1_dfm <- as.data.frame(apply(do.call(rbind, strsplit(g1_str, ',')), 2, as.numeric))
g2_dfm <- as.data.frame(apply(do.call(rbind, strsplit(g2_str, ',')), 2, as.numeric))
colnames(g1_dfm) <- paste0('g1.', 1:ncol(g1_dfm))
colnames(g2_dfm) <- paste0('g2.', 1:ncol(g2_dfm))

dfm <- cbind(g1_dfm, g2_dfm)
dfm <- dfm[apply(dfm, 1, function(x){sum(is.na(x))}) == 0,]
dfm <- dfm[apply(dfm, 1, sd) != 0,]

pca <- prcomp(t(dfm), center = T, scale. = T)
x <- pca$x
prop_var <- round(summary(pca)$importance[2,1:2]*100,0)

set_cols <- brewer.pal(8, 'Set1')
cols <- rep(set_cols, each=ncol(g1_dfm))

pdf(paste0('figure/', vs, '.IncLevel_PCA.pdf'), hei=7, wid=7)
layout(matrix(c(1,2),nrow=1), wid=c(5, 2))
par(mar=c(5,4,4,0))
par(xpd=TRUE)
plot(x[,1], x[,2], xlim=range(x[,1])*1.1, ylim=range(x[,2])*1.1, col=cols, cex=1, type='n',
     xlab=paste0('PC1 (',prop_var[1],'% of Variance)'),
     ylab=paste0('PC2 (',prop_var[2],'% of Variance)')
)
text(x[,1], x[,2], colnames(dfm), col=cols)
par(mar=c(5,0,4,0))
plot.new()
legend('left', paste0(unique(sub('\\..+', '', colnames(dfm))), ': ', groups), pch = 1, col=set_cols, bty='n')
dev.off()


