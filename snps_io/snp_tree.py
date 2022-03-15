from __future__ import division

import subprocess as sp, os, time, sys, argparse
import numpy as np
import scipy as sp

from skbio import DistanceMatrix
from skbio.tree import nj

def parse_args():
	""" Return dictionary of command line arguments
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS)

	parser.add_argument('program', help=argparse.SUPPRESS)

	parser.add_argument('--msa', type=str, dest='msa', required=True,
		help="""Path to multiple sequnece alignment file as input""")
	parser.add_argument('--out-dir', type=str, dest='out_dir',required=True,
		help="""Path to output dir (default same as input), two output files (snp_dist.tsv, snp_dist.tree) will be written""")
	parser.add_argument('--max-samples', type=int, metavar='INT', default=float('inf'),
		help="""Only use a subset of genomes or samples for building snp tree""")

	return vars(parser.parse_args())

def parse_msa(fpath, max_sample=float('inf')):
	"""
	this section reads the body of the parsnp xmfa file,
	the '=' marked the end of a multiple sequence alignment,
	a multiple sequence alignment contains sequences in fasta format for each sample
	each sequence fasta header was splitted to different atttributes, and store in Sequence class for description.
	each sequence was added to the Sequence class
	"""
	seqs = {}

	with open(fpath) as file:
		cur_id = ""
		for line in file:
			if line.startswith('='):
				pass
			elif line.startswith('>'):
				temp_items = line.rstrip().split(' ')
				cur_id = temp_items[0][1:]
				if cur_id not in seqs:
					if len(seqs.keys()) < max_sample:
						seqs[cur_id] = ""
			else:
				if cur_id in seqs:
					seqs[cur_id] += line.rstrip().upper()

	for id in seqs.keys():
		seqs[id] = np.array([np.fromstring(seqs[id], dtype='c')])

	return seqs

def main_vec():
	args = parse_args()

	genomes = parse_msa(args['msa'], args['max_samples'])

	try: os.makedirs(args['out_dir'])
	except: pass

	print("Count SNPs")
	dist_path = os.path.join(args['out_dir'], 'snp_dist.tsv')
	print("   path: %s" % dist_path)
	dist_file = open(dist_path, 'w')
	matrix = []

	occurs = []
	for id in genomes:
		occurs.append(genomes[id][0] != '-')
	occurs = np.array(occurs)

	for i, id in enumerate(genomes.keys()):
		occ_row = occurs[i]
		cooccs = occ_row & occurs

		diffs = []
		for sid in genomes:
			diffs.append(genomes[id][0] != genomes[sid][0])
		diffs = np.array(diffs)

		raw_counts = np.sum(diffs & cooccs, axis=1)
		norm_counts = raw_counts / np.sum(cooccs, axis=1)

		for j, sid in enumerate(genomes.keys()):
			dist_file.write('\t'.join([id, sid, str(raw_counts[j]), str(norm_counts[j])])+'\n')

		matrix.append(norm_counts)

	print("Build SNP tree")
	tree_path = os.path.join(args['out_dir'], 'snp_dist.tree')
	print("   path: %s" % tree_path)
	dm = DistanceMatrix(matrix, genomes.keys())
	tree = nj(dm, result_constructor=str)
	open(tree_path,'w').write(tree)

	print("\nDone!")

def main():
	args = parse_args()

	genomes = parse_msa(args['msa'], args['max_samples'])

	try: os.makedirs(args['out_dir'])
	except: pass

	print("Count SNPs")
	dist_path = os.path.join(args['out_dir'], 'snp_dist.tsv')
	print("   path: %s" % dist_path)
	dist_file = open(dist_path, 'w')
	matrix = []

	for id1 in genomes:
		array = []
		is_present1 = genomes[id1] != '-'
		for id2 in genomes:
			is_present2 = genomes[id2] != '-'
			is_diff = genomes[id1] != genomes[id2]
			co_occur = is_present1 & is_present2
			raw_count = (is_diff & co_occur).sum()

			norm_count = 0
			co_sum = co_occur.sum()

			if raw_count != 0 and co_sum != 0:
				norm_count = float(raw_count) / co_sum

			array.append(norm_count)
			dist_file.write('\t'.join([id1, id2, str(raw_count), str(norm_count)])+'\n')
		matrix.append(array)

	print("Build SNP tree")
	tree_path = os.path.join(args['out_dir'], 'snp_dist.tree')
	print("   path: %s" % tree_path)
	dm = DistanceMatrix(matrix, genomes.keys())
	tree = nj(dm, result_constructor=str)
	open(tree_path,'w').write(tree)

	print("\nDone!")

if __name__ == "__main__":
	main()
