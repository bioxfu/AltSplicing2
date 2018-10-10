SAM1=nramp1
SAM2=nramp1_Mn
SAM3=mdg10nramp1
SAM4=mdg10nramp1_Mn
SAM1_1=nramp1_1
SAM1_2=nramp1_2
SAM1_3=nramp1_3
SAM2_1=nramp1_Mn_1
SAM2_2=nramp1_Mn_2
SAM2_3=nramp1_Mn_3
SAM3_1=mdg10nramp1_1
SAM3_2=mdg10nramp1_2
SAM3_3=mdg10nramp1_3
SAM4_1=mdg10nramp1_Mn_1
SAM4_2=mdg10nramp1_Mn_2
SAM4_3=mdg10nramp1_Mn_3

BAM_DIR=~/Project/CFH0003_20180831/RNA-Seq/bam

## Arabidopsis
GTF=/cluster/home/xfu/Gmatic7/gene/tair10/tair10.gtf
ANNO=/cluster/home/xfu/Gmatic7/gene/tair10/tair10_gene_anno.tsv

## bam list
echo "$BAM_DIR/$SAM1_1.bam,$BAM_DIR/$SAM1_2.bam,$BAM_DIR/$SAM1_3.bam" > b1.txt
echo "$BAM_DIR/$SAM2_1.bam,$BAM_DIR/$SAM2_2.bam,$BAM_DIR/$SAM2_3.bam" > b2.txt
echo "$BAM_DIR/$SAM3_1.bam,$BAM_DIR/$SAM3_2.bam,$BAM_DIR/$SAM3_3.bam" > b3.txt
echo "$BAM_DIR/$SAM4_1.bam,$BAM_DIR/$SAM4_2.bam,$BAM_DIR/$SAM4_3.bam" > b4.txt

mkdir rMATS_out
rmats.py --b1 b1.txt --b2 b2.txt --gtf $GTF --od rMATS_out/${SAM1}_vs_${SAM2} -t paired --nthread 20 --readLength 101 --tstat 20 --libType fr-firststrand 
rmats.py --b1 b3.txt --b2 b4.txt --gtf $GTF --od rMATS_out/${SAM3}_vs_${SAM4} -t paired --nthread 20 --readLength 101 --tstat 20 --libType fr-firststrand 
rmats.py --b1 b3.txt --b2 b1.txt --gtf $GTF --od rMATS_out/${SAM3}_vs_${SAM1} -t paired --nthread 20 --readLength 101 --tstat 20 --libType fr-firststrand 
rmats.py --b1 b4.txt --b2 b2.txt --gtf $GTF --od rMATS_out/${SAM4}_vs_${SAM2} -t paired --nthread 20 --readLength 101 --tstat 20 --libType fr-firststrand 

## select some columns
mkdir rMATS_out_reformat
./script/reformat_rMATS_out.sh ${SAM1}_vs_${SAM2}
./script/reformat_rMATS_out.sh ${SAM3}_vs_${SAM4}
./script/reformat_rMATS_out.sh ${SAM3}_vs_${SAM1}
./script/reformat_rMATS_out.sh ${SAM4}_vs_${SAM2}

## merge tables
mkdir tables
Rscript script/merge_AS_out.R ${SAM1}_vs_${SAM2} ${SAM3}_vs_${SAM4} ${SAM3}_vs_${SAM1} ${SAM4}_vs_${SAM2}

## add gene annotation
Rscript script/AS_gene_anno.R $ANNO merge_A3SS_PValue0.05_deltaPSI0.1
Rscript script/AS_gene_anno.R $ANNO merge_A5SS_PValue0.05_deltaPSI0.1
Rscript script/AS_gene_anno.R $ANNO merge_MXE_PValue0.05_deltaPSI0.1
Rscript script/AS_gene_anno.R $ANNO merge_RI_PValue0.05_deltaPSI0.1
Rscript script/AS_gene_anno.R $ANNO merge_SE_PValue0.05_deltaPSI0.1

