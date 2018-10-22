## merge all possible alternative splicing (AS) events derived from GTF and RNA. 
fromGTF_PATH1=/data1/HJY/Project/Mecp2/KO/AltSplicing2/rMATS_out/WT_vs_KO/
fromGTF_PATH2=/data1/HJY/Project/Mecp2/T158M/AltSplicing2/rMATS_out/WT_vs_T158/
Rscript script/combine_bg_events.R $fromGTF_PATH1 $fromGTF_PATH2

## compare the DAS of KO and Mutant, PSI cutoff: 0.1, 0.15
KO=/data1/HJY/Project/Mecp2/KO/AltSplicing2/pairadise/WT_vs_KO/
MT=/data1/HJY/Project/Mecp2/T158M/AltSplicing2/pairadise/WT_vs_T158/

Rscript script/DAS_overlap_pvalue_odds.R 0.1 $KO $MT
Rscript script/DAS_overlap_pvalue_odds.R 0.15 $KO $MT
