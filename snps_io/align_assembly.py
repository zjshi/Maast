from __future__ import division

import numpy as np
from snps_io import vcf_io

class AlignAssembly:
	def __init__(self, algns=[], max_sites=float('inf'), gpos_offset=0):
		self.alignments = algns
		self.sample_ids = []

		self.chroms = np.array([])
		# self.char_mat = np.array([])
		self.global_pos = np.array([])
		self.local_pos = np.array([])
		self.ref_alleles = np.array([])
		self.alt_alleles = np.array([])
		self.third_alleles = np.array([])
		self.forth_alleles = np.array([])

		self.ref_prob_mat = np.array([])
		self.alt_prob_mat = np.array([])
		self.third_prob_mat = np.array([])
		self.forth_prob_mat = np.array([])

		self.freq_mat = np.array([])
		self.sample_presence = np.array([])
		self.prevalence = np.array([])
		self.ref_freqs = np.array([])
		self.alt_freqs = np.array([])
		self.third_freqs = np.array([])
		self.forth_freqs = np.array([])

		self.gpos_offset = gpos_offset
		self.is_spliced = False

		if self._check_():
			self.splice(max_sites)

		self.snps = []
		self.coords = []
		self.consensus_genome = ""

	def _check_(self):
		if len(self.alignments) > 2:
			if all([len(align.sample_ids) == len(self.alignments[0].sample_ids) for align in self.alignments]):
				return True
			else:
				return False
		else:
			return False

	def splice(self, max_sites=float('inf')):
		if (len(self.alignments) == 1):
			print("warning: cannot splice less than 2 alignments")

			first_align = self.alignments[0]
			self.sample_ids = first_align.sample_ids

			# self.char_mat = np.concatenate([algn.char_mat for algn in self.alignments], axis=1)
			self.chroms = np.array(np.repeat(first_align.chrom, first_align.ncols))
			self.local_pos = np.array(first_align.local_pos)

			self.ref_alleles = np.array(first_align.ref_alleles)
			self.alt_alleles = np.array(first_align.alt_alleles)
			self.third_alleles = np.array(first_align.third_alleles)
			self.forth_alleles = np.array(first_align.forth_alleles)

			self.ref_prob_mat = np.array(first_align.ref_prob_mat)
			self.alt_prob_mat = np.array(first_align.alt_prob_mat)
			self.third_prob_mat = np.array(first_align.third_prob_mat)
			self.forth_prob_mat = np.array(first_align.forth_prob_mat)

			self.freq_mat = np.array(first_align.freq_mat) # attention here
			self.sample_presence = np.array(first_align.sample_presence)
			self.prevalence = np.array(first_align.prevalence)
			self.ref_freqs = np.array(first_align.ref_freqs)
			self.alt_freqs = np.array(first_align.alt_freqs)
			self.third_freqs = np.array(first_align.third_freqs)
			self.forth_freqs = np.array(first_align.forth_freqs)

			self.global_pos = np.arange(len(self.chroms)) + self.gpos_offset

			self.is_spliced = True
			if max_sites < len(self.chroms):
				self.cut_short(max_sites)
		elif len(self.alignments) > 1 and all([len(align.sample_ids) == len(self.alignments[0].sample_ids) for align in self.alignments]):
			first_align = self.alignments[0]
			self.sample_ids = first_align.sample_ids

			# self.char_mat = np.concatenate([algn.char_mat for algn in self.alignments], axis=1)
			self.chroms = np.concatenate([np.repeat(algn.chrom, algn.ncols) for algn in self.alignments], axis=0)
			self.local_pos = np.concatenate([algn.local_pos for algn in self.alignments], axis=0)

			self.ref_alleles = np.concatenate([algn.ref_alleles for algn in self.alignments], axis=0)
			self.alt_alleles = np.concatenate([algn.alt_alleles for algn in self.alignments], axis=0)
			self.third_alleles = np.concatenate([algn.third_alleles for algn in self.alignments], axis=0)
			self.forth_alleles = np.concatenate([algn.forth_alleles for algn in self.alignments], axis=0)


			self.ref_prob_mat = np.concatenate([algn.ref_prob_mat for algn in self.alignments], axis=1)
			self.alt_prob_mat = np.concatenate([algn.alt_prob_mat for algn in self.alignments], axis=1)
			self.third_prob_mat = np.concatenate([algn.third_prob_mat for algn in self.alignments], axis=1)
			self.forth_prob_mat = np.concatenate([algn.forth_prob_mat for algn in self.alignments], axis=1)

			self.freq_mat = np.concatenate([algn.freq_mat for algn in self.alignments], axis=1)

			self.sample_presence = np.concatenate([algn.sample_presence for algn in self.alignments], axis=0)
			self.prevalence = np.concatenate([algn.prevalence for algn in self.alignments], axis=0)
			self.ref_freqs = np.concatenate([algn.ref_freqs for algn in self.alignments], axis=0)
			self.alt_freqs = np.concatenate([algn.alt_freqs for algn in self.alignments], axis=0)
			self.third_freqs = np.concatenate([algn.third_freqs for algn in self.alignments], axis=0)
			self.forth_freqs = np.concatenate([algn.forth_freqs for algn in self.alignments], axis=0)

			self.global_pos = np.arange(len(self.chroms)) + self.gpos_offset

			self.is_spliced = True

			if max_sites < len(self.chroms):
				self.cut_short(max_sites)
		else:
			print("errors: zero or uneven alignments")

		print("total number of sites: {}".format(len(self.chroms)))

		return self.is_spliced

	def cut_short(self, _max_sites):
		if self.is_spliced:
			if _max_sites < len(self.chroms):
				max_sites = int(_max_sites)

				self.chroms = self.chroms[:max_sites]
				self.global_pos = self.global_pos[:max_sites]
				self.local_pos = self.local_pos[:max_sites]
				self.ref_alleles = self.ref_alleles[:max_sites]
				self.alt_alleles = self.alt_alleles[:max_sites]
				self.freq_mat = self.freq_mat[:,:max_sites]

				self.ref_prob_mat = self.ref_prob_mat[:,:max_sites]
				self.alt_prob_mat = self.alt_prob_mat[:,:max_sites]
				self.third_prob_mat = self.third_prob_mat[:,:max_sites]
				self.forth_prob_mat = self.forth_prob_mat[:,:max_sites]

				self.sample_presence = self.sample_presence[:max_sites]
				self.prevalence = self.prevalence[:max_sites]
				self.ref_freqs = self.ref_freqs[:max_sites]
				self.alt_freqs = self.alt_freqs[:max_sites]
		else:
			print("warnings: impossible to cut short unspliced alignments, no changes has been done!")

	def id_core_genome(self, min_prev, min_alt_freq):
		print("min. prevalence: {}".format(min_prev))
		print("min. alt. frequency: {}".format(min_alt_freq))

		if self.is_spliced:
			prev_mask = (self.prevalence >= min_prev)
			snp_mask = (self.alt_freqs >= min_alt_freq) & (self.ref_alleles != b'N') & (self.ref_alleles != b'-')
			wildcard_mask = (self.ref_alleles != b'N') & (self.ref_alleles != b'-')

			# alt_freq_mask = ((1 - self.ref_freqs - self.alt_freqs) <= (min_alt_freq+0.000000001))


			# fake_mask = np.logical_not(alt_freq_mask)
			# print alt_freq_mask[fake_mask]
			# print (1 - self.ref_freqs - self.alt_freqs)[fake_mask]
			# print self.ref_freqs[fake_mask]
			# print self.alt_freqs[fake_mask]
			# print self.third_freqs[fake_mask]
			# print self.forth_freqs[fake_mask]

			self.consensus_genome = self.id_consensus_genome()

			shift_chroms = np.append(self.chroms[1:], self.chroms[-1])
			boundary_mask = np.logical_not((shift_chroms == self.chroms))
			#goodness_mask = (prev_mask & alt_freq_mask & wildcard_mask)
			goodness_mask = (prev_mask & wildcard_mask)

			self.coords = self.id_coordinates(boundary_mask, goodness_mask)

			print("masked by prev_mask: {}".format(np.sum(prev_mask)))
			print("masked by snp_mask: {}".format(np.sum(snp_mask)))
			# print "masked by alt_freq_mask: {}".format(np.sum(alt_freq_mask))
			print("masked by wildcard_mask: {}".format(np.sum(wildcard_mask)))

			calling_mask = goodness_mask & snp_mask
			self.snps = self.id_snps(calling_mask)
		else:
			sys.exit("premature call of id_core_genome, the multiple alignments were sliced yet.")

	def id_consensus_genome(self):
		if self.is_spliced:
			if len(self.ref_alleles) > 0:
				return b''.join([ref_allele for ref_allele in self.ref_alleles])
			else:
				return b''
		else:
			return b''

	def id_snps(self, calling_mask):
		if self.is_spliced:
			snps = []

			snp_chroms = self.chroms[calling_mask]
			snp_gb_pos = self.global_pos[calling_mask]
			snp_lc_pos = self.local_pos[calling_mask]
			snp_refs = self.ref_alleles[calling_mask]
			snp_alts = self.alt_alleles[calling_mask]
			snp_third = self.third_alleles[calling_mask]
			snp_forth = self.forth_alleles[calling_mask]

			snp_ref_prob_mat = self.ref_prob_mat[:,calling_mask]
			snp_alt_prob_mat = self.alt_prob_mat[:,calling_mask]
			snp_third_prob_mat = self.third_prob_mat[:,calling_mask]
			snp_forth_prob_mat = self.forth_prob_mat[:,calling_mask]

			snp_freqs = self.freq_mat[:,calling_mask]
			snp_presence = self.sample_presence[calling_mask]
			snp_prevs = self.prevalence[calling_mask]
			snp_ref_freqs = self.ref_freqs[calling_mask]
			snp_alt_freqs = self.alt_freqs[calling_mask]
			snp_third_freqs = self.third_freqs[calling_mask]
			snp_forth_freqs = self.forth_freqs[calling_mask]

			for i, chrom in enumerate(snp_chroms):
				var_id = str(snp_gb_pos[i])

				freq_row = snp_freqs[:,i]
				freq_row[snp_freqs[:,i] == None] = -1

				snp_ref_prob_row = snp_ref_prob_mat[:,i]
				snp_ref_prob_row[snp_ref_prob_mat[:,i] == None] = -1

				snp_alt_prob_row = snp_alt_prob_mat[:,i]
				snp_alt_prob_row[snp_alt_prob_mat[:,i] == None] = -1

				snp_third_prob_row = snp_third_prob_mat[:,i]
				snp_third_prob_row[snp_third_prob_mat[:,i] == None] = -1

				snp_forth_prob_row = snp_forth_prob_mat[:,i]
				snp_forth_prob_row[snp_forth_prob_mat[:,i] == None] = -1

				allele_mask = (np.array([snp_ref_prob_row.sum(), snp_alt_prob_row.sum(), snp_third_prob_row.sum(), snp_forth_prob_row.sum()]) > 0)

				alleles = np.array([snp_alts[i], snp_third[i], snp_forth[i]])
				alleles = alleles[allele_mask[1:]]

				if len(alleles) == 0:
					avail_alleles = b'.'
				else:
					avail_alleles = b','.join(alleles)

				snp = self._make_snp_(
					chrom, var_id, snp_lc_pos[i],
					snp_refs[i], snp_alts[i], snp_third[i], snp_forth[i], avail_alleles,
					len(self.sample_ids), snp_presence[i], round(snp_alt_freqs[i], 3),
					self.sample_ids, snp_ref_prob_row, snp_alt_prob_row, snp_third_prob_row, snp_forth_prob_row
				)
				snps.append(snp)

			self.snps = snps
			return snps
		else:
			return []

	def id_coordinates(self, boundary_mask, goodness_mask):
		if self.is_spliced:
			end_pos = np.array([])
			start_pos = np.array([])

			if len(self.alignments) > 1:
				end_pos = self.global_pos[boundary_mask]
				end_pos = np.append(end_pos, self.global_pos[-1])

				start_pos = np.array([self.global_pos[0]])
				shift_ends = end_pos[:-1] + 1
				start_pos = np.concatenate((start_pos, shift_ends))
			else:
				end_pos = np.array([self.global_pos[-1]])
				start_pos = np.array([self.global_pos[0]])

			bad_pos = self.global_pos[np.logical_not(goodness_mask)]

			rshift_bad_pos = bad_pos + 1
			lshift_bad_pos = bad_pos - 1

			start_pos = np.concatenate((start_pos, rshift_bad_pos))
			end_pos = np.concatenate((end_pos, lshift_bad_pos))

			start_pos = np.sort(start_pos)
			end_pos = np.sort(end_pos)

			good_region_mask = (start_pos <= end_pos)
			start_pos = start_pos[good_region_mask]
			end_pos = end_pos[good_region_mask]

			end_pos = np.sort(end_pos)

			coords = []
			for i, sp in enumerate(start_pos):
				coords.append({'chrom':self.chroms[sp-self.gpos_offset], 'start':sp, 'end':end_pos[i]})

			self.coords = coords
			return coords
		else:
			return []

	def _make_snp_(self, chrom, var_id, pos, ref, alt, third, forth, avail_alleles, NS, DP, AF, samp_ids, gp1, gp2, gp3, gp4):
		""" Format SNP for VCF """
		info = {}
		info['NS'] = NS
		info['DP'] = DP
		info['AF'] = AF

		dat_fmt = {}
		dat_fmt['GP1'] = gp1
		dat_fmt['GP2'] = gp2
		dat_fmt['GP3'] = gp3
		dat_fmt['GP4'] = gp4

		snp = vcf_io.SNP(chrom, var_id, pos, ref, alt, third, forth, avail_alleles, info, dat_fmt, samp_ids)

		return snp

def call_snps(aligns, max_sites, min_prev, snp_freq):
	"""
	Loop over each genomic site in each contig.
	For each site, fetch per-sample info from pileup files.
	Initialize GenomicSite object
	Determine site prevalence and allele frequency.
	Keep track of core-genome coordinates and SNPs in those regions.

	Args:
		max_sites:		int; max number of sites to process
		min_prev:		float; minimum prevalence for calling core sites
		snp_freq:		float; minimum minor allele frequency for snp calling
	"""

	aa = AlignAssembly(aligns, max_sites)

	if not aa.is_spliced:
		aa.splice()

	if max_sites < len(aa.chroms):
		aa.cut_short(max_sites)

	aa.id_core_genome(min_prev, snp_freq)

	return aa

def call_snps_iter(align_iterator, max_sites, min_prev, snp_freq):
	"""
	Loop over each genomic site in each contig.
	For each site, fetch per-sample info from pileup files.
	Initialize GenomicSite object
	Determine site prevalence and allele frequency.
	Keep track of core-genome coordinates and SNPs in those regions.

	Args:
		max_sites:		int; max number of sites to process
		min_prev:		float; minimum prevalence for calling core sites
		snp_freq:		float; minimum minor allele frequency for snp calling
	"""

	block_size = 100*1000
	counter = 0
	gb_pos = 0

	aligns = []
	for align in align_iterator:
		aligns.append(align)
		counter = counter + align.ncols

		if counter > block_size:
			aa = AlignAssembly(aligns, max_sites, gb_pos)

			if not aa.is_spliced:
				aa.splice()

			if max_sites < len(aa.chroms):
				aa.cut_short(max_sites)

			aa.id_core_genome(min_prev, snp_freq)

			for ali in aligns:
				gb_pos = gb_pos + ali.ncols

			aligns = []
			counter = 0

			yield aa


	if len(aligns) > 0:
		aa = AlignAssembly(aligns, max_sites, gb_pos)

		if not aa.is_spliced:
			aa.splice()

		if max_sites < len(aa.chroms):
			aa.cut_short(max_sites)

		aa.id_core_genome(min_prev, snp_freq)

		yield aa
