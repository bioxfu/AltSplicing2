## merge all possible alternative splicing (AS) events derived from GTF and RNA. 
fromGTF_PATH1=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_MeCP2_KO.14
fromGTF_PATH2=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KO
fromGTF_PATH3=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KI
Rscript script/combine_bg_events.R $fromGTF_PATH1 $fromGTF_PATH2 $fromGTF_PATH3

## compare the DAS of KO and Mutant, PSI cutoff: 0.1, 0.15
SAM1=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_MeCP2_KO.14/
SAM2=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KO/
SAM3=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KI/

Rscript script/DAS_overlap_pvalue_odds_3sample.R 0.1 $SAM1 $SAM2 $SAM3 MeCP2_KO RBFOX2_KO RBFOX2_KI

