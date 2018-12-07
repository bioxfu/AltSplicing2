# download phastCons and liftover chain files
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons20way/hg38.phastCons20way.bw
# wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

## install bwtool
# mkdir tools
# PREFIX=$PWD/tools
# git clone https://github.com/CRG-Barcelona/libbeato.git
# git clone https://github.com/CRG-Barcelona/bwtool.git
# cd libbeato/
# ./configure --prefix=$PREFIX CFLAGS="-g -O0 -I$PREFIX/include" LDFLAGS=-L$PREFIX/lib
# make
# make install
# cd ../bwtool/
# ./configure --prefix=$PREFIX CFLAGS="-g -O0 -I$PREFIX/include" LDFLAGS=-L$PREFIX/lib
# make
# make install
# cd ..

# plot phastCons for SE 
# grep '^SE|' figures/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp.txt|tr '|' '\t'|awk '{if($3=="+"){print $2"\t"$4"\t"$4+1"\tUP\t0\t+"} if($3=="-"){print $2"\t"$5-1"\t"$5"\tUP\t0\t-"}}' > test.up.bed
# grep '^SE|' figures/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp.txt|tr '|' '\t'|awk '{if($3=="-"){print $2"\t"$4"\t"$4+1"\tDN\t0\t-"} if($3=="+"){print $2"\t"$5-1"\t"$5"\tDN\t0\t+"}}' > test.dn.bed
# tools/bin/bwtool aggregate 150:30 test.up.bed hg38.phastCons7way.bw test.up.phastCons
# tools/bin/bwtool aggregate 30:150 test.dn.bed hg38.phastCons7way.bw test.dn.phastCons
# Rscrpt script/plot_phastCons.R

mkdir conservation
cat /data1/HJY/Project/Mecp2/KO/AltSplicing2_124/pairadise/WT_vs_KO/SE.MATS.JCEC.paired.compared.txt|awk '{print $4"\t"$6"\t"$7"\t"$2"\t"$3"\t"$5"\t"$0}'|grep -v 'exonStart_0base' > conservation/mm10_WT_vs_KO.bed

grep '^SE|' figures/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp.txt|tr '|' '\t'|awk '{print $2"\t"$4"\t"$5"\t"$10"\t"$11"\t"$3}' > tmp; paste tmp figures/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp.txt > conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp.bed
liftOver conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp.bed hg38ToMm10.over.chain.gz conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp_mm10.bed unMapped
bedtools intersect -a conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp_mm10.bed -b /home/xfu/Gmatic6/gene/mm10_vM13/gencode_mm10_vM13.gtf -wa -wb |awk '{if($10=="exon" && $2==($11-1) && $3==$12)print}'|cut -f1-14|sortBed|uniq > conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp_mm10_exons

head -1 /data1/HJY/Project/Mecp2/KO/AltSplicing2_124/pairadise/WT_vs_KO/SE.MATS.JCEC.paired.compared.txt > conservation/MeCP2_KO_vs_RBFOX2_KI_olp_Mecp2_KO_SE.MATS.JCEC.paired.compared.tsv
bedtools intersect -a conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp_mm10_exons -b conservation/mm10_WT_vs_KO.bed -wa -wb |cut -f7,22-  >> conservation/MeCP2_KO_vs_RBFOX2_KI_olp_Mecp2_KO_SE.MATS.JCEC.paired.compared.tsv
cut -f7,8,11,12,14 conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp_mm10_exons > conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KI_olp_mm10_exons.tsv


grep '^SE|' figures/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp.txt|tr '|' '\t'|awk '{print $2"\t"$4"\t"$5"\t"$10"\t"$11"\t"$3}' > tmp; paste tmp figures/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp.txt > conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp.bed
liftOver conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp.bed hg38ToMm10.over.chain.gz conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp_mm10.bed unMapped
bedtools intersect -a conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp_mm10.bed -b /home/xfu/Gmatic6/gene/mm10_vM13/gencode_mm10_vM13.gtf -wa -wb |awk '{if($10=="exon" && $2==($11-1) && $3==$12)print}'|cut -f1-14|sortBed|uniq > conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp_mm10_exons

head -1 /data1/HJY/Project/Mecp2/KO/AltSplicing2_124/pairadise/WT_vs_KO/SE.MATS.JCEC.paired.compared.txt > conservation/MeCP2_KO_vs_RBFOX2_KO_olp_Mecp2_KO_SE.MATS.JCEC.paired.compared.tsv
bedtools intersect -a conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp_mm10_exons -b conservation/mm10_WT_vs_KO.bed -wa -wb |cut -f7,22-  >> conservation/MeCP2_KO_vs_RBFOX2_KO_olp_Mecp2_KO_SE.MATS.JCEC.paired.compared.tsv
cut -f7,8,11,12,14 conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp_mm10_exons > conservation/AS_overlap_FDR0.01_PSI0.1_RPKM2_MeCP2_KO_vs_RBFOX2_KO_olp_mm10_exons.tsv
