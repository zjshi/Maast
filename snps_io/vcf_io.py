class SNP:
	def __init__(self, chrom, variant_id, pos, ref, alt, third=None, forth=None, avail_alleles=None, info=None, fmt=None, sample_ids=None):
		self.chrom = chrom
		self.var_id = variant_id
		self.pos = pos

		self.ref_allele = ref
		self.alt_allele = alt
		self.third_allele = third
		self.forth_allele = forth

		self.avail_alleles = avail_alleles

		if info is None:
			self.info = {}
			self.info['NS'] = -1
			self.info['DP'] = -1
			self.info['AF'] = -1
		else:
			self.info = info

		if fmt is None:
			self.format = {}
			self.format['GP1'] = ""
			self.format['GP2'] = ""
			self.format['GP3'] = ""
			self.format['GP4'] = ""
		else:
			self.format = fmt

		if sample_ids is None:
			self.sample_ids = []
		else:
			self.sample_ids = sample_ids

def format_header(sample_ids, cmdl):
	import time
	header = ""
	header += """##fileformat=VCFv4.1\n"""
	header += """##fileDate=%s\n""" % time.strftime("%Y-%m-%d %H:%M")
	header += """##source=https://github.com/snayfach/snpMLST\n"""
	header += """##command='%s'\n""" % cmdl
	header += """##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n"""
	header += """##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n"""
	header += """##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate Allele Frequency">\n"""
	header += """##FORMAT=<ID=GP1,Number=1,Type=Float,Description="Genotype 1 Probability">\n"""
	header += """##FORMAT=<ID=GP2,Number=1,Type=Float,Description="Genotype 2 Probability">\n"""
	header += """##FORMAT=<ID=GP3,Number=1,Type=Float,Description="Genotype 3 Probability">\n"""
	header += """##FORMAT=<ID=GP4,Number=1,Type=Float,Description="Genotype 4 Probability">\n"""

	col_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_ids
	header += """#%s\n""" % "\t".join(col_names)

	return header

def format_snp(snp):
	record = ""
	record += str(snp.chrom) + "\t" # CHROM
	record += str(snp.pos) + "\t" # POS
	record += str(snp.var_id) + "\t" # ID
	record += (snp.ref_allele + b"\t").decode() # REF
	record += (snp.avail_alleles + b"\t").decode() # ALT
	record += ".\t" # QUAL
	record += "PASS\t" # FILTER
	record += "%s\t" % format_info(snp) # INFO
	record += "%s\t" % ":".join(sorted(snp.format.keys())) # FORMAT
	record += "%s\n" % format_samples(snp) # GENOTYPES
	return record

def format_info(snp):
	return ";".join([key + "=" + str(value) for key, value in snp.info.items()])

def format_samples(snp):
	formats = sorted(snp.format.keys())
	indexes = range(len(snp.sample_ids))
	return "\t".join([":".join([str(snp.format[f][i]) for f in formats]) for i in indexes])

def write_vcf_header(snps, outdir, cmdl='unspecified'):
	import sys

	path = outdir+'/core_snps.vcf'
	if len(snps) > 0:
		with open(path, 'w') as file:
			file.write(format_header(snps[0].sample_ids, cmdl))

def write_vcf(snps, outdir, single_chrom_rep=False):
	import sys

	path = outdir+'/core_snps.vcf'
	if len(snps) > 0:
		with open(path, 'a') as file:
			for snp in snps:
				if single_chrom_rep is True:
					t_snp = snp
					t_snp.pos = t_snp.var_id
					file.write(format_snp(t_snp))
				else:
					file.write(format_snp(snp))
	else:
		print("Empty set of SNPs was found for the dataset, the file writing was skipped")

def write_coords_header(coords, out_dir):
	path = out_dir+'/coords.tsv'
	with open(path, 'w') as file:
		file.write('\t'.join(['chrom', 'start', 'end'])+'\n')

def write_coords(coords, out_dir):
	path = out_dir+'/coords.tsv'
	with open(path, 'a') as file:
		for d in coords:
			file.write('\t'.join([d['chrom'], str(d['start']), str(d['end'])])+'\n')

def merge_coords(coords, min_gap=1):
	if len(coords) > 1:
		merged_coords = []

		last_coord = coords[0]
		for i, coord in enumerate(coords[1:]):
			if coord['start'] - last_coord['end'] <= min_gap:
				if coord['chrom'] == last_coord['chrom']:
					new_coord = {'chrom':coord['chrom'], 'start':last_coord['start'], 'end':coord['end']}
					last_coord = new_coord
					continue

			merged_coords.append(last_coord)
			last_coord = coord

		merged_coords.append(last_coord)

		return merged_coords
	else:
		return coords

def write_genome(genome, out_dir):
	path = out_dir+'/consensus.fna'
	with open(path, 'w') as file:
		file.write('>consensus\n'+genome+'\n')

##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
