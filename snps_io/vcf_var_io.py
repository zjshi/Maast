class SNP:
	def __init__(self, chrom, variant_id, pos, ref="", alt="", info=None, fmt=None, sample_ids=None):
		self.chrom = chrom
		self.var_id = variant_id
		self.pos = pos

		self.ref_allele = ""
		self.alt_allele = ""

		if info is None:
			self.info = {}
			self.info['NS'] = -1
			self.info['DP'] = -1
			self.info['AF'] = -1

		if fmt is None:
			self.format = {}
			self.format['AF'] = ""

		if sample_ids is None:
			self.sample_ids = []


def format_header(sample_ids):
	import time
	header = ""
	header += """##fileformat=VCFv4.1\n"""
	header += """##fileDate=%s\n""" % time.strftime("%Y-%m-%d %H:%M")
	header += """##source=https://github.com/snayfach/snpMLST\n"""
	header += """##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n"""
	header += """##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n"""
	header += """##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate Allele Frequency">\n"""
	header += """##FORMAT=<ID=AF,Number=1,Type=Float,Description="Alternate Allele Frequency">\n"""
	header += """##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n"""
	header += """#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT %s\n""" % " ".join(sample_ids)
	return header

def format_snp(snp):
	record = ""
	record += str(snp.chrom) + "\t" # CHROM
	record += str(snp.pos) + "\t" # POS
	record += str(snp.var_id) + "\t" # ID
	record += snp.ref_allele + "\t" # REF
	record += snp.alt_allele + "\t" # ALT
	record += ".\t" # QUAL
	record += "PASS\t" # FILTER
	record += "%s\t" % format_info(snp) # INFO
	record += "%s\t" % ":".join(snp.format.keys()) # FORMAT
	record += "%s\n" % format_samples(snp) # GENOTYPES
	return record

def format_info(snp):
	return ";".join([key + "=" + str(value) for key, value in snp.info.items()])

def format_samples(snp):
	formats = snp.format.keys()
	indexes = range(len(snp.sample_ids))
	return "\t".join([":".join([str(snp.format[f][i]) for f in formats]) for i in indexes])

def write_vcf(snps, outdir):
	path = outdir+'/core_snps.vcf'
	if len(snps) > 0:
		with open(path, 'w') as file:
			file.write(format_header(snps[0].sample_ids))
			for snp in snps:
				file.write(format_snp(snp))
	else:
		print "Empty set of SNPs was found for the dataset, the file writing was skipped"
