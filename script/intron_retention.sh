source activate rMATS
git clone https://github.com/Xinglab/SQUID

#DIR=$HOME/Project/HJY/RBM17/HepG2/fastq
DIR=$HOME/Project/HJY/RBM17/K562/fastq

GTF=$HOME/Gmatic7/gene/human/GRCh38_v29.gtf
CTRL1=$DIR/ctrl_1
CTRL2=$DIR/ctrl_2
EXPT1=$DIR/expt_1
EXPT2=$DIR/expt_2
KALLISTO=$HOME/Gmatic7/genome/human/kallisto/GRCh38
STAR=$HOME/Gmatic7/genome/human/STAR/

nohup python ../SQUID.py --GTF $GTF --fastq ${EXPT1}_R1_paired.fastq.gz:${EXPT1}_R2_paired.fastq.gz,${EXPT2}_R1_paired.fastq.gz:${EXPT2}_R2_paired.fastq.gz,${CTRL1}_R1_paired.fastq.gz:${CTRL1}_R2_paired.fastq.gz,${CTRL2}_R1_paired.fastq.gz:${CTRL2}_R2_paired.fastq.gz \
--index_kallisto $KALLISTO --index_star $STAR --anchor 8 --length 100 --lib first --read P --Cal All --c1 0.05 --p 30 --Comparison ./Comparison --analysis U -o ./bam1 --resume true 

#nohup python ../SQUID.py --GTF $GTF --fastq ${EXPT1}_R1_paired.fastq.gz:${EXPT1}_R2_paired.fastq.gz,${EXPT2}_R1_paired.fastq.gz:${EXPT2}_R2_paired.fastq.gz,${CTRL1}_R1_paired.fastq.gz:${CTRL1}_R2_paired.fastq.gz,${CTRL2}_R1_paired.fastq.gz:${CTRL2}_R2_paired.fastq.gz \
#--index_kallisto $KALLISTO --index_star $STAR --anchor 8 --length 100 --lib first --read P --Cal All --c1 0.05 --p 30 --Comparison ./Comparison --analysis U -o ./bam1 --resume true \
#--align ./bam1/align/S0Aligned.sortedByCoord.out.bam,./bam1/align/S1Aligned.sortedByCoord.out.bam,./bam1/align/S2Aligned.sortedByCoord.out.bam,./bam1/align/S3Aligned.sortedByCoord.out.bam \
#--FPKM  ./bam1/FPKM/transcript_exp.txt &

# manually run to fix the bugs
cat ./bam1/FPKM/transcript_exp.txt|sed -r 's/\|.+\|//' > tmp; mv tmp ./bam1/FPKM/transcript_exp.txt
python /cluster/home/xfu/Project/HJY/RBM17/SQUID/bin/Transcript.py --FPKM  ./bam1/FPKM/transcript_exp.txt --intron ./bam1/gtf_files/Intron_transcript.txt --output ./bam1/DEXSeq_files/DEXSeq_compare1/FPKM_intron_transcript.txt
python patch.py
Rscript /cluster/home/xfu/Project/HJY/RBM17/SQUID/bin/RP_value.R ./bam1/DEXSeq_files/DEXSeq_compare1/Diff_compare1_intron_temp.txt 100 ./bam1/DEXSeq_files/DEXSeq_compare1/rank_product_test_compare1.txt ./bam1/Result/Diff_compare1_intron_PI.txt 0.05 0.1 ./bam1/Result/Increase_compare1_intron_PI.txt ./bam1/Result/Decrease_compare1_intron_PI.txt
Rscript patch.R
