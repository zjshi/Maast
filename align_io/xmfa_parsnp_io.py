from __future__ import division
import numpy as np
from snp_mlst.align_io import seq_ali

class Sequence:
	def __init__(self, line):
		values = line.rstrip().lstrip('>').lstrip().split()
		self.index = values[0].split(':')[0]
		self.start = int(values[0].split(':')[1].split('-')[0])
		self.end = int(values[0].split(':')[1].split('-')[1])
		self.length = self.end - self.start + 1
		self.strand = values[1]
		self.chrom = values[2]
		self.seq = ''

def extract_genome_info(fpath):
	# check version
	with open(fpath) as file:
		version = next(file).rstrip().split(' ', 1)[-1]
		if version != 'Parsnp v1.1':
			sys.exit("\nError: expected XMFA version 'Parsnp v1.1' but got '%s'\n" % version)

	"""
	this section reads the header of the parsnp xmfa file,
	header are the lines beginning with ## and #,
	map genome index to info,
	keys (field in code) = ['SequenceLength', 'SequenceFile', 'SequenceHeader']
	"""
	genome_info = {}
	with open(fpath) as file:
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

	return genome_info

def parse(fpath, max_sample=float('inf')):

	genome_info = extract_genome_info(fpath)

	"""
	this section reads the body of the parsnp xmfa file,
	the '=' marked the end of a multiple sequence alignment,
	a multiple sequence alignment contains sequences in fasta format for each sample
	each sequence fasta header was splitted to different atttributes, and store in Sequence class for description.
	each sequence was added to the Sequence class
	"""
	alns = []
	with open(fpath) as file:
		last_aln = None
		cur_aln = seq_ali.Alignment()
		for line in file:
			if line.startswith('#'):
				continue
			elif line.startswith('='):
				cur_aln.nseqs = len(cur_aln.seqs)
				cur_aln.ncols = len(cur_aln.seqs[0].seq)
				cur_aln.chrom = cur_aln.seqs[0].chrom
				cur_aln.update()
				last_aln = cur_aln
				cur_aln = seq_ali.Alignment()
				alns.append(last_aln)
			elif len(cur_aln.seqs) <= max_sample:
				if line.startswith('>'):
					seq = Sequence(line)
					seq.id = genome_info[seq.index]['SequenceFile']
					cur_aln.seqs.append(seq)
				else:
					cur_aln.seqs[-1].seq += line.rstrip().upper()
			else:
				if len(cur_aln.seqs) > max_sample:
					cur_aln.seqs.pop()
				pass
	return alns

def iter_parse(fpath, max_sample=float('inf')):
	genome_info = extract_genome_info(fpath)

	"""
	this section reads the body of the parsnp xmfa file,
	the '=' marked the end of a multiple sequence alignment,
	a multiple sequence alignment contains sequences in fasta format for each sample
	each sequence fasta header was splitted to different atttributes, and store in Sequence class for description.
	each sequence was added to the Sequence class
	"""

	with open(fpath) as file:
		last_aln = None
		cur_aln = seq_ali.Alignment()
		for line in file:
			if line.startswith('#'):
				continue
			elif line.startswith('='):
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
						seq = Sequence(line)
						seq.id = genome_info[seq.index]['SequenceFile']
						cur_aln.seqs.append(seq)
					else:
						cur_aln.seqs[-1].seq += line.rstrip().upper()
				else:
					if len(cur_aln.seqs) > max_sample:
						cur_aln.seqs.pop()
					pass
