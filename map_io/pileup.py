from __future__ import division

import csv, gzip, os, sys
import numpy as np
import pandas as pd

from time import time
from operator import itemgetter

class PileupAlignment:
	def __init__(self, chrom, chrom_len, chrom_pileups, min_depth):
		self.desc = '' #not used for now
		# chrom has identifier like 'cluster123', which mark a division on the core-genome, which is a conserved region found cross all reference genomes
		self.chrom = chrom
		self.pileups = chrom_pileups
		self.n_samples = len(chrom_pileups.keys())
		self.sample_ids = sorted(chrom_pileups.keys())
		self.ncols = chrom_len # number of sites

		self.gentle_check()
		"""
		attributes below were described in function update()
		"""
		self.local_pos = np.arange(self.ncols)

		self.counts_list = self.make_counts_list()
		self.pooled_counts = self.make_pooled_counts()
		self.pooled_depth = self.pooled_counts[0:4,:].sum(0)

		self.freq_mat = []

		self.ref_alleles = []
		self.alt_alleles = []
		self.third_alleles = []
		self.forth_alleles = []

		self.ref_prob_mat = []
		self.alt_prob_mat = []
		self.third_prob_mat = []
		self.forth_prob_mat = []

		self.sample_presence = []
		self.ref_freqs = []
		self.alt_freqs = []
		self.third_freqs = []
		self.forth_freqs = []

		self.prevalence = []
		self.aligned_pctg = []

		self.update()

	def gentle_check(self):
		gentle_msg = "something might have gone horribly wrong"

		if self.chrom is None:
			sys.stderr.write("chrom name is None, {}".format(gentle_msg))

		if self.n_samples <= 1:
			sys.stderr.write("one or zero sequences for encapsulation, {}".format(gentle_msg))

		if any(len(sample_id) == 0 for sample_id in self.sample_ids):
			sys.stderr.write("weird sample ids were found, {}".format(gentle_msg))

		if any([any(self.pileups[sample_id][:,1] > self.ncols) for sample_id in self.sample_ids]):
			sys.stderr.write("sequences were aligned out of the scope of chrom, {}".format(gentle_msg))

	def _make_pooled_counts(self):
		templ = np.repeat(0, 5 * self.ncols).reshape((5, self.ncols))

		for samp_id in self.sample_ids:
			plps = self.pileups[samp_id]
			plp_all_pos = np.array([plp.pos for plp in plps])
			plp_all_counts = np.transpose(np.array([plp.allele_list() for plp in plps]))

#			print plp_all_counts
#			print templ
#
#			print len(plp_all_pos)
#			print plp_all_counts.shape
#			print templ.shape
#			print templ[:,plp_all_pos-1].shape

			templ[:,plp_all_pos-1] += plp_all_counts

		return templ

	def make_allele_probs(self):
		n = len(self.sample_ids)

	def make_pooled_counts(self):
		templ = np.repeat(np.int32(0), 5 * self.ncols).reshape((5, self.ncols))

		for counts in self.counts_list:
			templ += counts

		return templ

	def make_counts_list(self):
		templs = []

		for samp_id in self.sample_ids:
			sub_templ = np.repeat(np.int32(0), 5 * self.ncols).reshape((5, self.ncols))

			plps = self.pileups[samp_id]

			plp_all_pos = plps[:,1].astype(int)
			plp_all_counts = np.transpose(plps[:,4:]).astype(np.int32)

			sub_templ[:,plp_all_pos-1] += plp_all_counts

			templs.append(sub_templ)

		return templs

	def update(self, min_depth=1):

		"""
		char_template: complete set of chars for each site on the sequences, from which the ref allele and alt allele will be selected
		"""
		char_template = np.array([
			np.repeat('A', self.pooled_counts.shape[1]),
			np.repeat('T', self.pooled_counts.shape[1]),
			np.repeat('G', self.pooled_counts.shape[1]),
			np.repeat('C', self.pooled_counts.shape[1])
			# np.repeat('N', self.pooled_counts.shape[1])
		])

		"""
		sorting the counts of different chars at each site,
		then select the indices of chars whose counts are in top 2,
		then using the indices of chars to select the chars,
		then the top 1 char (the char with highest count) at each site is considered as ref allele for now,
		then the char with the second highest count at each site is considered as alt allele for now,
		for those sites that have only one char, the counts of other chars are zeroes, so theorectically any of them can be selected by program, but in reality '-' will be selected, no exception was found.
		"""
		count_inds_mat = self.pooled_counts[0:4,:].argsort(axis=0)
		top2_inds = count_inds_mat[-4:,]
		top2_char_mat = np.choose(top2_inds, char_template)

		self.ref_alleles = top2_char_mat[3,:]
		self.alt_alleles = top2_char_mat[2,:]
		self.third_alleles = top2_char_mat[1,:]
		self.forth_alleles = top2_char_mat[0,:]

		# print self.ref_alleles
		# print self.alt_alleles

		"""
		frequency matrix has the same shape as the char matrix
		it is initialize to have None only,
		in the end, it is used to store the presence/absence of ref/alt allele for each site cross all samples with the following rules:
		- presence of ref allele: 1
		- absence of ref allele: 0
		- presence of N or -: None
		"""

		self.freq_mat = np.repeat(None, self.n_samples*self.ncols).reshape((self.n_samples, self.ncols))

		self.ref_freq_mat = np.repeat(None, self.n_samples*self.ncols).reshape((self.n_samples, self.ncols))
		self.alt_freq_mat = np.repeat(None, self.n_samples*self.ncols).reshape((self.n_samples, self.ncols))
		self.third_freq_mat = np.repeat(None, self.n_samples*self.ncols).reshape((self.n_samples, self.ncols))
		self.forth_freq_mat = np.repeat(None, self.n_samples*self.ncols).reshape((self.n_samples, self.ncols))

		self.ref_prob_mat = np.repeat(np.float16(0), self.n_samples * self.ncols).reshape((self.n_samples, self.ncols))
		self.alt_prob_mat = np.repeat(np.float16(0), self.n_samples * self.ncols).reshape((self.n_samples, self.ncols))
		self.third_prob_mat = np.repeat(np.float16(0), self.n_samples * self.ncols).reshape((self.n_samples, self.ncols))
		self.forth_prob_mat = np.repeat(np.float16(0), self.n_samples * self.ncols).reshape((self.n_samples, self.ncols))

		# print self.counts_list

		for i, count_mat in enumerate(self.counts_list):
			top2_count_mat = np.choose(top2_inds, count_mat[0:4,:])
			self.ref_freq_mat[i,:] = top2_count_mat[3,:]
			self.alt_freq_mat[i,:] = top2_count_mat[2,:]
			self.third_freq_mat[i,:] = top2_count_mat[1,:]
			self.forth_freq_mat[i,:] = top2_count_mat[0,:]

			# freq_sum = self.alt_freq_mat[i,:] + self.ref_freq_mat[i,:] + self.third_freq_mat[i,:] + self.forth_freq_mat[i,:]
			freq_sum = count_mat[0:4,:].sum(0)



			inval_mask = (freq_sum < min_depth)
			# print freq_sum
			# print inval_mask

			zero_mask = (freq_sum == 0)
			freq_sum[zero_mask] = 1

			self.freq_mat[i,:] = self.alt_freq_mat[i,:] / freq_sum
			self.freq_mat[i,inval_mask] = -2

			self.ref_prob_mat[i,:] = np.true_divide(self.ref_freq_mat[i,:], freq_sum)
			self.alt_prob_mat[i,:] = np.true_divide(self.alt_freq_mat[i,:], freq_sum)
			self.third_prob_mat[i,:] = np.true_divide(self.third_freq_mat[i,:], freq_sum)
			self.forth_prob_mat[i,:] = np.true_divide(self.forth_freq_mat[i,:], freq_sum)

			self.ref_prob_mat[i,inval_mask] = -2
			self.alt_prob_mat[i,inval_mask] = -2
			self.third_prob_mat[i,inval_mask] = -2
			self.forth_prob_mat[i,inval_mask] = -2

		# print self.freq_mat

		# sys.exit()

		"""
		the sample presence here only sum up the counts of A, T, G, C for each site, leave out the N and -,
		it facilitate the calculation of prevalence of the site on certain sequence alignment position
		ref and alt allele frequencies were calculated with denominator of sample_presence rather than the number of sample.
		naturally, there is another route to calculate them.
		"""
		top2_pool_count_mat = np.choose(top2_inds, self.pooled_counts)
		zero_mask = (self.pooled_depth == 0)
		temp_depth = self.pooled_depth
		temp_depth[zero_mask] = 1

		self.ref_freqs = top2_pool_count_mat[3,:] / temp_depth
		self.alt_freqs = top2_pool_count_mat[2,:] / temp_depth
		self.third_freqs = top2_pool_count_mat[1,:] / temp_depth
		self.forth_freqs = top2_pool_count_mat[0,:] / temp_depth

		self.ref_freqs[zero_mask] = 0
		self.alt_freqs[zero_mask] = 0
		self.third_freqs[zero_mask] = 0
		self.forth_freqs[zero_mask] = 0

		# print self.ref_freqs
		# print self.alt_freqs

		self.sample_presence = ((self.alt_freq_mat > 0) | (self.ref_freq_mat > 0) | (self.third_freq_mat > 0) | (self.forth_freq_mat > 0)).sum(0)
		self.prevalence = self.sample_presence/self.n_samples

		unaligned_masks = (self.pooled_counts[4,:] != 0)
		self.aligned_pctg = 1 - (unaligned_masks.sum(0) / self.n_samples)

def pileup_parse(plp_csv_path, uniq_chroms):
	#pileup_file = open(plp_csv_path)
	#pileup_reader = csv.reader(pileup_file, delimiter='\t')

	pileup_reader = pd.read_csv(plp_csv_path, sep=',', names=None).as_matrix()

	pileup = {}
	for chrom in uniq_chroms:
		chrom_mask = (pileup_reader[:,0] == chrom)
		pileup[chrom] = pileup_reader[chrom_mask,:]

	return pileup

def monolithic_parse(plp_csv_dir, map_genome):
	timeit = time()
	sys.stderr.write("\t[fetch contigs] start fetch contigs. \n")
	# fetch contig lengths
	contig_ids = []
	contig_lengths = []
	records = open(map_genome).read().split('>')[1:]
	for record in records:
		contig_id = record.split('\n', 1)[0].split()[0]
		contig_seq = record.split('\n', 1)[1].replace('\n', '')
		contig_ids.append(contig_id)
		contig_lengths.append([contig_id, len(contig_seq)])

	contig_ids = list(set(contig_ids))
	contig_lengths = sorted(contig_lengths)
	sys.stderr.write("\t[fetch contigs] elapsed seconds: %s\n" % round(time()-timeit,2) )

	# open pileups
	pileups = {}
	ext = '.pileup.csv'

	timeit = time()
	sys.stderr.write("\t[pileup parse] start to parse pileup. \n")
	for pileup_file in os.listdir(plp_csv_dir):
		sample_id = pileup_file[0:-len(ext)]
		pileup_path = os.path.join(plp_csv_dir, pileup_file)
		pileups[sample_id] = pileup_parse(pileup_path, contig_ids)
		sys.stderr.write("\t[pileup parse] finish loading {}\n".format(pileup_file))
	sys.stderr.write("\t[pileup parse] elapsed seconds: %s\n" % round(time()-timeit,2) )

	# sys.exit()

	timeit = time()
	sys.stderr.write("\t[identify sites] start identify sites. \n")
	min_depth = 1
	alns = []
	# loop over contigs
	for chrom, length in contig_lengths:
		chrom_pileups = {}
		for samp_id in pileups.keys():
			if chrom not in pileups[samp_id]:
				sys.stderr.write("{} not in the pileup set: {}\n".format(chrom, samp_id))
			else:
				chrom_pileups[samp_id] = pileups[samp_id][chrom]

		if len(chrom_pileups.keys()) != len(pileups.keys()):
			sys.stderr.write("chrom {} has not been aligned across all samples, so skipped\n".format(chrom))
		else:
			aln = PileupAlignment(chrom, length, chrom_pileups, min_depth)
			aln.update()
			alns.append(aln)

	sys.stderr.write("\t[identify sites] elapsed seconds: %s\n" % round(time()-timeit,2) )

	return alns
