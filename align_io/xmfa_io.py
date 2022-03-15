import sys
from operator import itemgetter

class Alignment:
	def __init__(self):
		self.desc = ''
		self.seqs = []

	def fetch_columns(self):
		# ncols is the length of first seq in Alignment
		# Can be inferred from this block of code that all seq in Alignment have the same length
		# seq.id = genome_info[seq.index]['SequenceFile']
		
		for col_index in range(self.ncols):
			pos = col_index + 1
			chars = [seq.seq[col_index] for seq in self.seqs]
			sample_ids = [seq.id for seq in self.seqs]
			yield AlignColumn(pos, chars, sample_ids)

class AlignColumn:
	def __init__(self, pos, chars, sample_ids):
		self.chrom = None
		self.pos = pos
		self.chars = chars
		self.sample_ids = sample_ids
		self.total_samples = len(chars)
		self.pooled_counts = {'A':0, 'T':0, 'G':0, 'C':0}
		self.pool_counts()
		self.pooled_depth = sum(self.pooled_counts.values())
		self.present_samples = sum(self.pooled_counts.values())
		self.call_snp()
		self.prev = self.present_samples / float(self.total_samples)
		self.allele_freqs = self.genotype()

	def pool_counts(self):
		""" Pool read-counts to 4 alleles """
		for char in self.chars:
			if char in self.pooled_counts:
				self.pooled_counts[char] += 1

	def call_snp(self):
		""" Identify major and minor alleles """
		counts = sorted(list(self.pooled_counts.items()), key=itemgetter(1), reverse=True)
		if self.present_samples == 0:
			self.cons_allele = 'N'
			self.cons_count = 0
			self.cons_freq = 0.0
			self.alt_freq = 0.0
		else:
			self.cons_allele, self.cons_count = counts.pop(0)
			self.cons_freq = self.cons_count/float(self.present_samples)
			if len(counts) > 0:
				self.alt_allele, self.alt_count = counts[0]
				self.alt_freq = self.alt_count/float(self.present_samples)
			else:
				self.alt_freq = 0.0

	def consensus(self):
		from collections import Counter
		return Counter(self.chars).most_common(1)[0][0]

	def percent_aligned(self, nseqs):
		""" Compute the % of genomes with observed data """
		gaps = self.chars.count('-')
		missing = self.chars.count('N')
		return (nseqs - gaps - missing) / float(nseqs)

	def genotype(self):
		genotypes = []
		for char in self.chars:
			if char == self.cons_allele:
				genotypes.append(0)
			elif char == self.alt_allele:
				genotypes.append(1)
			else:
				genotypes.append(None)
		return genotypes


class Sequence:
	def __init__(self, line):
		values = line.rstrip().lstrip('>').lstrip().split()
		self.index = values[0].split(':')[0]
		self.start = int(values[0].split(':')[1].split('-')[0])
		self.end = int(values[0].split(':')[1].split('-')[1])
		self.length = self.end - self.start + 1
		self.strand = values[1]
		self.seq = ''


def parse(fpath):

	# check version
	with open(fpath) as file:
		version = next(file).rstrip().split(' ', 1)[-1]
		if version != 'Parsnp v1.1':
			sys.exit("\nError: expected XMFA version 'Parsnp v1.1' but got '%s'\n" % version)

	# map genome index to info
	# keys = ['SequenceLength', 'SequenceFile', 'SequenceHeader']
	with open(fpath) as file:
		genome_info = {}
		header = ''
		for line in file:
			if line.startswith('##'): header += line.lstrip('#')
			elif line.startswith('#'): continue
			else: break
		for h in header.split('SequenceIndex ')[1:]:
			genome_index, info_string = h.rstrip('\n').split('\n', 1)
			genome_info[genome_index] = {}
			for info_record in info_string.split('\n'):
				field, value = info_record.split(' ', 1)
				genome_info[genome_index][field] = value

	# yield alignment blocks
	with open(fpath) as file:
		last = None
		current = Alignment()
		for line in file:
			if line.startswith('#'):
				continue
			elif line.startswith('='):
				current.nseqs = len(current.seqs)
				current.ncols = len(current.seqs[0].seq)
				last = current
				current = Alignment()
				yield last
			elif line.startswith('>'):
				seq = Sequence(line)
				seq.id = genome_info[seq.index]['SequenceFile']
				current.seqs.append(seq)
			else:
				current.seqs[-1].seq += line.rstrip().upper()
