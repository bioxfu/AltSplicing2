## merge all possible alternative splicing (AS) events derived from GTF and RNA. 
fromGTF_PATH1=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_MeCP2_KO.14
fromGTF_PATH2=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KO
fromGTF_PATH3=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KI
Rscript script/combine_bg_events.R $fromGTF_PATH1 $fromGTF_PATH2 $fromGTF_PATH3

## compare the DAS of KO and Mutant, PSI cutoff: 0.1, 0.15
SAM1=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_MeCP2_KO.14/
SAM2=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KO/
SAM3=/data1/HJY/Project/20181121JY/AltSplicing2/rMATS_out/WT_vs_RBFOX2_KI/
FDR=0.01
PSI=0.1
RPKM=/data1/HJY/Project/20181121JY/RNA-Seq/table_8samples/RPKM_table_FDR0.05_FC1.5_all.tsv
RPKM_cutoff=2

Rscript script/DAS_overlap_pvalue_odds_3sample.R $FDR $PSI $SAM1 $SAM2 $SAM3 MeCP2_KO RBFOX2_KO RBFOX2_KI $RPKM $RPKM_cutoff

