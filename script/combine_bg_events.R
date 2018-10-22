all_events <- NULL

paths <- commandArgs(T)

for (p in paths) {
  for (AS in c('SE', 'RI', 'MXE', 'A3SS', 'A5SS', 'novelEvents.SE', 'novelEvents.RI', 'novelEvents.MXE', 'novelEvents.A3SS', 'novelEvents.A5SS')) {
    KO <- read.table(paste0(p, '/fromGTF.', AS, '.txt'), header = T, sep = '\t')
    if (AS == 'SE' || AS == 'novelEvents.SE') {
      KO$event <- paste(KO$chr, KO$strand, KO$exonStart_0base, KO$exonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, sep='|')
    }
    if (AS == 'RI' || AS == 'novelEvents.RI') {
      KO$event <- paste(KO$chr, KO$strand, KO$riExonStart_0base, KO$riExonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, sep='|')
    }
    if (AS == 'MXE' || AS == 'novelEvents.MXE') {
      KO$event <- paste(KO$chr, KO$strand, KO$X1stExonStart_0base, KO$X1stExonEnd, KO$X2ndExonStart_0base, KO$X2ndExonEnd, KO$upstreamES, KO$upstreamEE, KO$downstreamES, KO$downstreamEE, sep='|')
    }
    if (AS == 'A3SS' || AS == 'novelEvents.A3SS') {
      KO$event <- paste(KO$chr, KO$strand, KO$longExonStart_0base, KO$longExonEnd, KO$shortES, KO$shortEE, KO$flankingES, KO$flankingEE, sep='|')
    }
    if (AS == 'A5SS' || AS == 'novelEvents.A5SS') {
      KO$event <- paste(KO$chr, KO$strand, KO$longExonStart_0base, KO$longExonEnd, KO$shortES, KO$shortEE, KO$flankingES, KO$flankingEE, sep='|')
    }
    all_events <- c(all_events, KO$event)
  }
}

write.table(unique(sort(all_events)), 'all_possible_AS_events', quote = F, col.names = F, row.names = F)
