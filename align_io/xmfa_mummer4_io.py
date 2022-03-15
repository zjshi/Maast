from __future__ import division
import numpy as np
import copy
from align_io import seq_ali

class SimpleSequence:
	def __init__(self, indx, id="", seq=""):
		self.index = indx
		self.id = id
		self.seq = seq
		self.chrom = ""

def parse(fpath, max_sample=float('inf')):
		"""
		this section reads the body of the parsnp xmfa file,
		the '=' marked the end of a multiple sequence alignment,
		a multiple sequence alignment contains sequences in fasta format for each sample
		each sequence fasta header was splitted to different atttributes, and store in Sequence class for description.
		each sequence was added to the Sequence class
		"""
		alns = []
		total_len = 0
		with open(fpath) as file:
			last_aln = None
			cur_aln = seq_ali.Alignment()
			indx = 0
			for line in file:
				if line.startswith('='):
					cur_aln.nseqs = len(cur_aln.seqs)
					cur_aln.ncols = len(cur_aln.seqs[0].seq)
					cur_aln.chrom = cur_aln.seqs[0].chrom
					print(cur_aln.nseqs)
					cur_aln.update()
					last_aln = cur_aln
					cur_aln = seq_ali.Alignment()
					alns.append(last_aln)
					total_len = total_len + last_aln.ncols
				else:
					if len(cur_aln.seqs) <= max_sample:
						if line.startswith('>'):
							seq = SimpleSequence(indx)
							temp_items = line.rstrip().split(' ')
							seq.id = temp_items[0][1:]
							seq.chrom = temp_items[1]
							cur_aln.seqs.append(seq)
							indx = indx + 1
						else:
							cur_aln.seqs[-1].seq += line.rstrip().upper()
					else:
						if len(cur_aln.seqs) > max_sample:
							cur_aln.seqs.pop()
						pass

		print("total length of alignments: {}".format(total_len))
		return alns

def iter_parse(fpath, max_sample=float('inf')):
		"""
		this section reads the body of the parsnp xmfa file,
		the '=' marked the end of a multiple sequence alignment,
		a multiple sequence alignment contains sequences in fasta format for each sample
		each sequence fasta header was splitted to different atttributes, and store in Sequence class for description.
		each sequence was added to the Sequence class
		"""
		lines = []
		n_align = 0
		with open(fpath) as file:
			for line in file:
				lines.append(line)
				if line.startswith('='):
					n_align = n_align + 1

		if n_align > 1:
			last_aln = None
			cur_aln = seq_ali.Alignment()
			indx = 0
			for line in lines:
				if line.startswith('='):
					cur_aln.nseqs = len(cur_aln.seqs)
					cur_aln.ncols = len(cur_aln.seqs[0].seq)
					cur_aln.chrom = cur_aln.seqs[0].chrom

					cur_aln.update()
					last_aln = cur_aln
					cur_aln = seq_ali.Alignment()
					yield last_aln
				else:
					if len(cur_aln.seqs) <= max_sample:
						if line.startswith('>'):
							seq = SimpleSequence(indx)
							temp_items = line.rstrip().split(' ')
							seq.id = temp_items[0][1:]
							seq.chrom = temp_items[1]
							cur_aln.seqs.append(seq)
							indx = indx + 1
						else:
							cur_aln.seqs[-1].seq += line.rstrip().upper()
					else:
						if len(cur_aln.seqs) > max_sample:
							cur_aln.seqs.pop()
						pass
		else:
			max_iter_stride = 200*1000

			last_aln = None
			cur_aln = seq_ali.Alignment()
			indx = 0

			for line in lines:
				if line.startswith('='):
					cur_aln.nseqs = len(cur_aln.seqs)
					cur_aln.ncols = len(cur_aln.seqs[0].seq)
					cur_aln.chrom = cur_aln.seqs[0].chrom

					if cur_aln.ncols <= max_iter_stride:
						# print cur_aln.ncols
						cur_aln.update()
						yield cur_aln
					else:
						for sp in range(0, cur_aln.ncols, max_iter_stride):
							next_aln = seq_ali.Alignment()
							for seq in cur_aln.seqs:
								temp_seq = copy.deepcopy(seq)
								if sp+max_iter_stride <= cur_aln.ncols+1000:
									temp_seq.seq = temp_seq.seq[sp:(sp+max_iter_stride)]
								else:
									temp_seq.seq = temp_seq.seq[sp:]
								next_aln.seqs.append(temp_seq)
								#print next_aln.seqs[-1].chrom

							next_aln.nseqs = len(next_aln.seqs)
							next_aln.ncols = len(next_aln.seqs[0].seq)
							next_aln.chrom = next_aln.seqs[0].chrom
							next_aln.update()
							yield next_aln
				else:
					if len(cur_aln.seqs) <= max_sample:
						if line.startswith('>'):
							seq = SimpleSequence(indx)
							temp_items = line.rstrip().split(' ')
							seq.id = temp_items[0][1:]
							seq.chrom = temp_items[1]
							cur_aln.seqs.append(seq)
							indx = indx + 1
						else:
							cur_aln.seqs[-1].seq += line.rstrip().upper()
					else:
						if len(cur_aln.seqs) > max_sample:
							cur_aln.seqs.pop()
						pass
