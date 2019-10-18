import sys

cutoff = 0.1
gene2name = {}
gene2type = {}

#with open('~/Gmatic7/gene/human/GRCh38_v29_id2name2biotype') as f:
with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip('\n').split('\t')
		gene2name[lst[0]] = lst[1]
		gene2type[lst[0]] = lst[2]

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip('\n').split('\t')
		inc_cnt, skp_cnt, inc_body = lst[8], lst[9], lst[10]
		pi_junc1, pi_junc2 = lst[-2:]
		if pi_junc1 != 'PI_Junction':
			if (not 'NA' in pi_junc1) and (not 'NA' in pi_junc2):
				s1 = [float(x) for x in pi_junc1.split(',')]
				s2 = [float(x) for x in pi_junc2.split(',')]
				delta_pi1 = (s1[2] + s1[3])/2 - (s1[0] + s1[1])/2
				delta_pi2 = (s2[2] + s2[3])/2 - (s2[0] + s2[1])/2
				if abs(delta_pi1) > cutoff:
					lst[1] = lst[1]+'\t'+','.join([gene2name[x] for x in lst[1].split(',')])+'\t'+','.join([gene2type[x] for x in lst[1].split(',')])
					lst[8] = inc_cnt.replace(',', '\t')
					lst[9] = skp_cnt.replace(',', '\t')
					lst[10] = inc_body.replace(',', '\t')
					lst[-2] = pi_junc1.replace(',', '\t')
					lst[-1] = pi_junc2.replace(',', '\t')
					print('%s\t%s\t%s' % ('\t'.join(lst), delta_pi1, delta_pi2))
		else:
			lst[1] = '\t'.join(('Gene_id', 'Gene_name', 'Biotype'))
			lst[8] = '\t'.join(('Inclusion_counts.ctr1','Inclusion_counts.ctr2','Inclusion_counts.expt1','Inclusion_counts.expt2'))
			lst[9] = '\t'.join(('Skipping_counts.ctr1','Skipping_counts.ctr2','Skipping_counts.expt1','Skipping_counts.expt2'))
			lst[10] = '\t'.join(('Inclusion_counts_with_intron_body.ctr1','Inclusion_counts_with_intron_body.ctr2','Inclusion_counts_with_intron_body.expt1','Inclusion_counts_with_intron_body.expt2'))
			lst[-2] = '\t'.join(('PI_Junction.ctr1','PI_Junction.ctr2','PI_Junction.expt1','PI_Junction.expt2'))
			lst[-1] = '\t'.join(('PI_JunctionIntron.ctr1','PI_JunctionIntron.ctr2','PI_JunctionIntron.expt1','PI_JunctionIntron.expt2'))
			print('%s\t%s\t%s' % ('\t'.join(lst), 'delta_PI_Junction', 'delta_PI_JunctionIntron'))
