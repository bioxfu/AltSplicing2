up <- read.table('test.up.phastCons')
dn <- read.table('test.dn.phastCons')

plot(1:180, up$V2, xlim=c(1, 370), ylim=c(0,0.8), type = 'l', lwd=3, ylab='Average phastCons score', xaxt='n', xlab='Exon')
lines(191:370, dn$V2, lwd=3)
axis(1, at= c(seq(1,150,10), 150), labels = c(seq(-150,-10, 10), -1), las=2)
axis(1, at= c(seq(1,150,10), 150)+220, labels = paste0('+', c(1, seq(10, 150, 10))), las=2)
rect(xleft = 150, ybottom = -0.05, xright = 221, ytop = 0, col = 'black', border = NA)
rect(xleft = 180, ybottom = -0.05, xright = 191, ytop = 0, col = 'white', border = NA)
