output_S = 'bam1/DEXSeq_files/DEXSeq_compare1'
output_path = 'bam1/Result'
output = 'bam1'
ss1=[1,2]
ss2=[3,4]
fr1 = open("%s/DEXSeq.txt" %(output_S))
fr2 = open("%s/DEXSeq_counts.txt" %(output_S))
fr3 = open("%s/intron_PI.txt" % (output_path))
fr4 = open("%s/rMATS_files/rMATS_%s_Junction/rMATS_Result.txt" %(output,'compare1'))

fw = open("%s/Diff_%s_intron_temp.txt" %(output_S, 'compare1'),"w")
#fw = open("test_out.txt","w")
fw.write("Intron_id\tGene_id\tStrand\tChr\tStart\tEnd\tAnnotated\tAttributes\tInclusion_counts_SAMPLE1\tSkipping_counts_SAMPLE1\tInclusion_counts_SAMPLE2\tSkipping_counts_SAMPLE2\tInclusion_length\tSkipping_length\tPValue_rMATS\tFDR_rMATS\tPI_Junction_SAMPLE1\tPI_Junction_SAMPLE2\tDiff_PI_Junction\tObserved_counts_SAMPLE1\tExpected_counts_SAMPLE1\tObserved_counts_SAMPLE2\tExpected_counts_SAMPLE2\tPValue_DEXSeq\tFDR_DEXSeq\tPI_Density_SAMPLE1\tPI_Density_SAMPLE2\tDiff_PI_Density\n")		
info1 = {}
info2 = {}
info3 = {}
info4 = {}
len1 = 0
len2 = 0	

for line in fr1:
	lst = line.strip().split('\t')
	info1[lst[2]] = lst

for line in fr2:
	lst = line.strip().split('\t')
	info2[lst[0]] = lst

for line in fr3:
	lst = line.strip().split('\t')
	info3[lst[0]] = lst

for line in fr4:
	lst = line.strip().split('\t')
	info4[lst[0]] = lst
	len1 = lst[5]
	len2 = lst[6]	

fr1.close()
fr2.close()
fr3.close()
fr4.close()

inc1 = ['0'] * len(ss1)
skp1 = ['0'] * len(ss1)
inc2 = ['0'] * len(ss2)
skp2 = ['0'] * len(ss2)
PI1=["NA"] * len(ss1)
PI2= ["NA"] * len(ss1)

for x in sorted(info3.keys()):
	if x != 'Intron_id':
		a1 = info1.get(x, 'NA')
		a2 = info2.get(x, 'NA')
		a3 = info3.get(x, 'NA')
		a4 = info4.get(x, 'NA')
		if a4 =='NA':
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\t%s\t%s\tNA\t%s\t%s\t%s\n" %("\t".join(a3[0:8]), ",".join(inc1),",".join(skp1),",".join(inc2),",".join(skp2),len1,len2, ",".join(PI1),",".join(PI2),"\t".join(a2[1:5]), "\t".join(a1[6:8]),"\t".join(a2[7:10])))
		else:
			fw.write("%s\t%s\t%s\t%s\t%s\n" %("\t".join(a3[0:8]), "\t".join(a4[1:]),"\t".join(a2[1:5]), "\t".join(a1[6:8]),"\t".join(a2[7:10])))

fw.close()
