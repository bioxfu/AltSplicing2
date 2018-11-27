SAM1=WT
SAM2=T158
SAM1_1=WT_1_18R370
SAM1_2=WT_2_18R372
SAM1_3=WT_3_18R374
SAM1_4=WT_4_18R376
SAM2_1=T158_1_18R371
SAM2_2=T158_2_18R373
SAM2_3=T158_3_18R375
SAM2_4=T158_4_18R377

BAM_DIR=/data1/HJY/Project/Mecp2/T158M/RNA-Seq/bam

## Human
#GTF=/cluster/home/xfu/Gmatic6/gene/hg38_v26/gencode_hg38_v26.gtf
## Mouse
GTF=/cluster/home/xfu/Gmatic6/gene/mm10_vM13/gencode_mm10_vM13.gtf

## bam list
echo "$BAM_DIR/$SAM1_1.bam,$BAM_DIR/$SAM1_2.bam,$BAM_DIR/$SAM1_3.bam,$BAM_DIR/$SAM1_4.bam" > b1.txt
echo "$BAM_DIR/$SAM2_1.bam,$BAM_DIR/$SAM2_2.bam,$BAM_DIR/$SAM2_3.bam,$BAM_DIR/$SAM2_4.bam" > b2.txt

mkdir rMATS_out
rmats.py --b1 b1.txt --b2 b2.txt --gtf $GTF --od rMATS_out/${SAM1}_vs_${SAM2} -t paired --nthread 20 --readLength 150 --tstat 20 --libType fr-firststrand 

## Do paired-test
mkdir -p pairadise/${SAM1}_vs_${SAM2}
Rscript ./script/AS_paired_compare.R

## Plot PCA using IncLevel
mkdir figure
Rscript script/IncLevel_PCA.R ${SAM1}_vs_${SAM2}

