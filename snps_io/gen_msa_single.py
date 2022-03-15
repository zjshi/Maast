import sys, argparse

from Bio import SeqIO

def parse_args():
	""" Return dictionary of command line arguments
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS)

	parser.add_argument('--delta', type=str, dest='delta_path', required=True,
		help="""Path to delta file indicating the whole genome alignment information. This delta file is also the output of mummer4""")
	parser.add_argument('--ref-seq', type=str, dest='ref_path', required=True,
		help="""Path to reference genome file as input in multiple fasta format""")
	parser.add_argument('--qry-seq', type=str, dest='qry_path', required=True,
		help="""Path to query genome file as input in multiple fasta format""")
	parser.add_argument('--ref-name', type=str, dest='ref_name', required=True,
		help="""Specify reference genome name""")
	parser.add_argument('--qry-name', type=str, dest='qry_name', required=True,
		help="""Specify query genome name""")
	parser.add_argument('--out', type=str, dest='out', default="/dev/out",
		help="""Path to output (double genome alignment) will be written in MSA formart""")

	return vars(parser.parse_args())


def rc(seq):
	base_map = {
		'A': 'T', 'a': 'T',
		'C': 'G', 'c': 'G',
		'G': 'C', 'g': 'C',
		'T': 'A', 't': 'A',
		'N': 'N', 'n': 'N'
	} 

	return ''.join([base_map[c] for c in seq[::-1]])

def read_genome(genome_path):
	ordered_chroms = []
	genome_seqs = dict()
	
	for seq in SeqIO.parse(genome_path, "fasta"):
		ordered_chroms.append(seq.id)
		genome_seqs[seq.id] = seq.seq

	return genome_seqs, ordered_chroms

def parse_delta(delta_path):
	align_blocs = []

	with open(delta_path, 'r') as fh:
		fh.readline()
		fh.readline()

		r_tag = ""
		q_tag = ""
		
		r_len = ""
		q_len = ""

		bloc = []

		for line in fh:
			if line[0] == '>':
				items = line[1:].rstrip().split(' ')
				r_tag = items[0]
				q_tag = items[1]

				r_len = int(items[2])
				q_len = int(items[3])
			else:
				if ' ' in line:
					items = line.rstrip().split(' ')
					bloc = [ int(item) for item in items[:4] ]
				else:
					diff = line.rstrip()
					
					if diff == '0':
						align_blocs.append([r_tag, q_tag, r_len, q_len] + bloc)
						bloc = []
					else:
						bloc.append(int(diff))
	
	return align_blocs

def main():
	args = parse_args()

	ref_genome, ref_chroms = read_genome(args['ref_path'])
	qry_genome, qry_chroms  = read_genome(args['qry_path'])

	align_blocs = parse_delta(args['delta_path'])

	aligned_qry = dict()

	for chrom_id in ref_genome.keys():
		aligned_qry[chrom_id] = '-' * len(ref_genome[chrom_id])

	for bloc in align_blocs:
		r_tag = bloc[0]
		q_tag = bloc[1]

		r_len = bloc[2]
		q_len = bloc[3]

		r_start = bloc[4]
		r_end = bloc[5]

		assert r_end > r_start

		q_start = bloc[6]
		q_end = bloc[7]

		q_seq = ""
		if q_end < q_start:
			q_start = q_len - q_start + 1
			q_end = q_len - q_end + 1

			q_seq = rc(qry_genome[q_tag])[q_start-1:q_end]
		else:
			q_seq = qry_genome[q_tag][q_start-1:q_end]


		pos = 0
		for diff in bloc[8:]:
			if diff < 0:
				pos = pos+abs(diff)
				q_seq = q_seq[:pos-1] + q_seq[pos:]
				pos = pos - 1
			else:
				pos = pos + diff
				q_seq = q_seq[:pos-1] + '-' + q_seq[pos-1:]

		# print [r_start, r_end, r_end - r_start, q_start, q_end, q_end - q_start,  len(q_seq)]

if __name__ == "__main__":
	main()
