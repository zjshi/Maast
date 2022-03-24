#!/usr/bin/env python3

from __future__ import division

import sys, os, time, argparse, shutil, hashlib, math
import numpy as np
from operator import itemgetter

from Bio import SeqIO

from snps_io import id_genome_clusters, id_centroid
from snps_io import vcf_io, concat_alleles, gen_msa, align_assembly

from db_io import build_db

def get_data_type():
	""" Get program specified by user (species, genes, or snps) """
	import sys
	if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
		cmd = '%s ' % (os.path.basename(sys.argv[0]))
		print('usage: %s <module> [options]' % cmd)
		print('')
		print('description: identify and genotype core-genome snps from <module>')
		print('')
		print('modules:')
		print('	end_to_end         Run full Maast pipeline from begining to end')
		print('	genomes            perform multiple alignment of genomes to call core-genome SNPs')
		print('	db                 build kmer database targeting snps')
		print('	genotype           call core-genome SNPs for single genomes and isolate sequencing data')
		print('	tree               build SNP tree using identified genotypes')
		print('')
		print("use '%s <module> -h' for usage on a specific command" % cmd)
		print('')
		quit()
	elif sys.argv[1] not in ['end_to_end', 'genomes', 'db', 'genotype', 'tree']:
		sys.exit("\nError: invalid subcommand\n\nSupported subcommand: genomes, db, genotype, end_to_end, tree\n")
	else:
		return sys.argv[1]

def parse_args():

	data_type = get_data_type()

	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		add_help=False,
		usage=argparse.SUPPRESS
		)

	parser.add_argument('data_type', help=argparse.SUPPRESS)

	if data_type == 'end_to_end':
		end2end_input = parser.add_argument_group('end2end_input')
		end2end_input.add_argument('--in-dir', type=str, metavar='PATH',required=True,
			help = """Path to directory of raw-read-files in FASTQ format (.fastq or .fq; gzipped or not) or whole-genome sequences in FASTA format (.fna, .fa, .fsa or .fasta). (Required)""")

	io = parser.add_argument_group('input/output')
	io.add_argument('--out-dir', type=str, metavar='PATH', required=True,
		help="""Directory to store output (required)""")

	if data_type in ['genomes']: 
		io.add_argument('--fna-dir', type=str, metavar='PATH', required=True,
			help = """Path to directory of genomes in FASTA format (required)""")
	
	if data_type in ['genomes', 'end_to_end']:
		io.add_argument('--rep-fna', type=str, metavar='PATH', default=None,
			help = """Path to the reference genome serving as the template for whole genome alignment. If provided, Maast will not identify and use centroid genome for reference (default None)""")
		io.add_argument('--skip-align', action='store_true', default=False,
			help = """skip whole genome sequence or short read alignment, only applicable when alignment has already been done (default False)""")
		io.add_argument('--has-completeness', action='store_true', default=False,
			help = """Toggle for specifying completeness for supplied genomes sequenes. If toggled on, it requries to supply either --completeness or --completeness-list (default False)""")
		io.add_argument('--completeness', type=float, metavar='FLOAT', default=None,
			help = """Single completeness value for all genomes sequenes (i.e. all genomes have the same completeness) (default False)""")
		io.add_argument('--completeness-list', type=str, metavar='PATH', default=None,
			help = """Path to list of pairs of genome file name and completeness value, separated by tab character. (note: genome file names should have no duplicates, and should cover all files specified in --fna-dir) (default None)""")
		io.add_argument('--missing-ratio', type=float, metavar='FLOAT', default=0.05,
			help = """Parameter defining the missing ratio of core sites even when completeness is 1 (default 0.05)""")
		io.add_argument('--min-pid', type=float, metavar='FLOAT', default=0,
			help = """Parameter defining the minimal identity for including each aligned block, [0, 100] (default 0)""")
		io.add_argument('--min-aln-len', type=int, metavar='INT', default=10,
			help = """Parameter defining the minimal length for including each aligned block (default 10)""")
		io.add_argument('--max-pid-delta', type=float, metavar='FLOAT', default=0.1,
			help = """Parameter defining the maximum identity gap between identity of each aligned block and whole-genome ANI, all alignments with identity less than ANI * (1 - delta) will be purged, [0, 1] (default 10)""")
		io.add_argument('--mem', action='store_true', default=False,
			help = """calling SNPs by genomic segment, option for memory saving (default False)""")

	if data_type in ['genomes', 'end_to_end']:
		prep = parser.add_argument_group('preprocessing')
		prep.add_argument('--keep-redundancy', action='store_true', default=False,
			help="""If toggled on, Maast will skip redundancy removal and move on with all input genomes (default=False)""")
		prep.add_argument('--skip-centroid', action='store_true', default=False,
			help="""If toggled on, Maast will not attempt to identify and use centroid genome for reference (default=False)""")
		prep.add_argument('--sketch-k', type=int, metavar='INT', default=21,
			help="""k-mer size for building Mash sketch (default=21)""")
		prep.add_argument('--sketch-size', type=int, metavar='INT', default=5000,
			help="""The number of k-mers per Mash sketch (default=5000)""")
		prep.add_argument('--precut', type=float, metavar='FLOAT', default=0.05,
			help="""Limit searches among pair of genomes with distance smaller than the provided value (default=0.05)""")
		prep.add_argument('--start-cutoff', type=float, metavar='FLOAT', default=0.02,
			help="""The cutoff from which Maast will start to search a distance cutoff, which generate the good number of genome clusters and tag genomes based on a given MAF (default=0.02)""")
		prep.add_argument('--end-cutoff', type=float, metavar='FLOAT', default=0.0001,
			help="""Similiar to --start-cutoff, the cutoff at which Maast will end the search for a distance cutoff. This value should be smaller than --start-cutoff (default=0.0001)""")
		prep.add_argument('--range-factor', type=float, metavar='FLOAT', default=1.2,
			help="""This factor times the minimum number of genomes needed for a given MAF will create the upper bound of a range satisfying the search. It should be larger than 1 (default=1.2)""")

	if data_type in ['genomes', 'end_to_end']:
		snps = parser.add_argument_group('snp-calling')
		snps.add_argument('--max-sites', type=int, metavar='INT', default=float('inf'),
			help="""Maximum genomic sites to parse (use all); useful for testing (default=inf)""")
		snps.add_argument('--min-prev', type=float, metavar='FLOAT', default=1.0,
			help="""Minimum prevalence (default=1.0)""")
		snps.add_argument('--snp-freq', type=float, metavar='FLOAT', default=0.01,
			help="""Minimum minor allele frequency for SNP calling (default=0.01)""")
		snps.add_argument('--max-samples', type=int, metavar='INT', default=float('inf'),
			help="""Only use a subset of genomes or metagenomes for snp calling (default=inf)""")

	if data_type in ['db', 'end_to_end']:
		db = parser.add_argument_group('db-building')
	if data_type in ['db']:
		db.add_argument('--ref-genome', type=str, dest='ref_genome', required=True,
			help="""Path to reference genome sequence file (required)""")
		db.add_argument('--vcf', type=str, dest='vcf', required=True,
			help="""Path to a vcf file describing core snps/genetic variants called based on multiple sequence alignments (required)""")
		db.add_argument('--msa', type=str, dest='msa', required=True,
			help="""Path to multiple sequence alignment file (required)""")
		db.add_argument('--tag-fna-list', type=str, dest='tag_list', required=True,
			help="""Path to a list of paths to the tag genomes (FASTA format) which are included in multiple sequence alignment file (required)""")
		db.add_argument('--fna-dir', type=str, dest='fna_dir', default=None,
			help="""Path to a list of paths to the tag genomes (FASTA format) which are included in multiple sequence alignment file (default=None)""")
		db.add_argument('--coords', type=str, dest='coords', default=None,
			help="""Path to core genome block coordinate file (default=None)""")

	if data_type in ['db', 'end_to_end']:
		db.add_argument('--genome-name', dest='genome_name', type=str, default='100000',
			help="""Name of the core-genome corresponding to INPUT. Should be six digits with the first digit in [1, 9] (default=100000)""")
		db.add_argument('--overwrite', dest='overwrite', action='store_true', help="""Overwrite existing output files""")
		db.add_argument('--kmer-type', dest='kmer_type', default='all',
			choices=['all', 'center'],
			help="""
                    Choose type of kmers to be fetched
                    all: all elligible kmers 
                        1) covered snp at any position
                        and 2) do not cover any bad sites (e.g. N or -)
                        and 3) were well contained on its coordinate division (default)
                    center: all kmers whose target snps was at their centers.""")
		db.add_argument('--snp-cover', dest='snp_type', default='all',
			choices=['all', 'l1-tags', 'l2-tags'],
			help="""
                    Choose object to kmerize
                    all: all snps from the cluster will be attempted for kmer search; most kmers (default)
                    l1-tags: only representative snps from all snp blocks will be attempted
                    l2-tags: only representative snps from representative snp blocks will be attempted; fewest kmers
                    * note: all kmers must uniquely match an allele and intersect >= 1 SNP""")

	if data_type in ['genotype', 'end_to_end']:
		genotype_input = parser.add_argument_group('genotype_input')

	if data_type in ['genotype']:
		genotype_input.add_argument('--in-dir', type=str, metavar='PATH',required=True,
			help = """Path to directory of raw-read-files in FASTQ format (.fastq or .fq; gzipped or not) or whole-genome sequences in FASTA format (.fna, .fa, .fsa or .fasta) (required)""")
		genotype_input.add_argument('--ref-genome', type=str, dest='ref_genome', required=True,
			help="""Path to reference genome sequence file (required)""")
		genotype_input.add_argument('--db', type=str, metavar='PATH', dest='kmer_db_path', required=True,
			help = """Path to directory of raw-read-files in FASTQ format (.fastq or .fq; gzipped or not) or whole-genome sequences in FASTA format (.fna, .fa, .fsa or .fasta) (required)""")
		genotype_input.add_argument('--vcf', type=str, dest='vcf', required=True,
			help="""Path to a vcf file describing core snps/genetic variants called based on multiple sequence alignments (required)""")
		single_genome = parser.add_argument_group('genome-genotyping')
		single_genome.add_argument('--min-pid', type=float, metavar='FLOAT', default=0,
			help = """Parameter defining the minimal identity for including each aligned block, [0, 100] (default=0)""")
		single_genome.add_argument('--min-aln-len', type=int, metavar='INT', default=10,
			help = """Parameter defining the minimal length for including each aligned block (default=10)""")
		single_genome.add_argument('--max-pid-delta', type=float, metavar='FLOAT', default=0.1,
			help = """Parameter defining the maximum identity gap between identity of each aligned block and whole-genome ANI, all alignments with identity less than ANI * (1 - delta) will be purged, [0, 1] (default=0.1)""")

	if data_type in ['genotype', 'end_to_end']:
		genotype_input.add_argument('--merge-pairs', action='store_true', default=False,
			help = """Flag to merge paired raw reads files in <in-dir>; indicated by ext '_1*' and '_2*'""")

		align = parser.add_argument_group('reads-genotyping')
		align.add_argument('--mode', default='very-sensitive',
			choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
			help = """Alignment speed/sensitivity (default=very-sensitive)""")
		align.add_argument('--max-reads', type=int, metavar='INT',
			help = """Maximum # reads to use from each FASTQ file (default=None; use all)""")

	if data_type in ['genomes', 'genotype', 'end_to_end']:
		io.add_argument('--subset-list', type=str, metavar='PATH', default=None, 
			help = """Path to file contains the names of the fullset or subset of the files in the input directory. Files not in the list will not be included for snp calling (default=None; use all)""")

	if data_type in ['tree']:	
		tree_io = parser.add_argument_group('tree_io')
		tree_io.add_argument('--input-list', type=str, dest='input_list', required=True,
			help="""A list of input pairs. Each pair per row contains a path to a genotype result file generated from Maast genotype command and a unique name of the file. (required)
                    The path and name must be separated by a tab.
                    Example
                    /file/path/1	name1
                    /file/path/2	name2
                    /file/path/3    name3
                    ...""")
		tree_io.add_argument('--min-sites', type=int, dest='min_sites_per_sample', default=1000,
			help="""Minimum SNP sites. Any allele sequence with a number of non-empty sites lower than this value will not be included (default=1000)""")
		tree_io.add_argument('--max-gap-ratio', type=float, dest='max_gap_ratio', default=0.5,
			help="""Maximum ratio of gaps. Any allele sequence with a ratio of gap higher than this value will not be included (default=0.5)""")
		tree_io.add_argument('--min-site-prev', type=float, dest='min_site_prev', default=0.9,
			help="""Minimum site prevalence. Any site with an actual allele presents in a fraction of sequences lower than this value will not be included (default=0.9)""")
		tree_io.add_argument('--min-MAF', type=float, dest='min_maf', default=0.01,
			help="""Minimum allele frequency. Any site with MAF lower than this value will not be included (default=0.01)""")
		tree_io.add_argument('--min-MAC', type=float, dest='min_mac', default=1,
			help="""Minimum allele count. Any site with MAC lower than this value will not be included (default=1)""")
		tree_io.add_argument('--min-depth', type=float, dest='min_depth', default=1,
			help="""Minimum read depth. Any site supported by a number of reads lower than this value will not be included. This option is only for genotypes identified from sequencing reads. Default value is 1 and any value >1 will effectively exclude all whole genome assemblies from analysis. Caution is advised (default=1)""")

	misc = parser.add_argument_group('misc')
	misc.add_argument("-h", "--help", action="help",
		help="""Show this help message and exit""")
	misc.add_argument('--threads', type=int, metavar='INT', default=1,
		help="""Number of CPUs to use (default=1)""")

	args = vars(parser.parse_args())

	args['data_type'] = data_type

	return args

def run_command(cmd, env=None):
	import subprocess as sp
	if env:
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, env=env)
	else:
		p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = p.communicate()
	if p.returncode != 0:
		err_msg =  "\nError: the following returned non-zero status: '%s':\n" % cmd
		err_msg += "\n%s" % err
		sys.exit(err_msg)
	else:
		return out, err

def parallel(function, argument_list, threads):
	""" Based on: https://gist.github.com/admackin/003dd646e5fadee8b8d6 """
	import multiprocessing as mp
	import signal
	import time

	def init_worker():
		signal.signal(signal.SIGINT, signal.SIG_IGN)

	pool = mp.Pool(int(threads), init_worker)

	try:
		results = []
		for arguments in argument_list:
			p = pool.apply_async(function, args=arguments)
			results.append(p)
		pool.close()

		while True:
			if all(r.ready() for r in results):
				return [r.get() for r in results]
			time.sleep(1)

	except KeyboardInterrupt:
		pool.terminate()
		pool.join()
		sys.exit("\nKeyboardInterrupt")

def reformat_sequence_headers(args):
	"""
	Reformat sequence headers in input genomes to prevent parsnp from crashing
	"""
	import Bio.SeqIO
	if 'fna_dir' in args:
		try: os.makedirs(args['out_dir']+'/temp/genomes')
		except: pass
		for file in os.listdir(args['fna_dir']):
			infile = open(args['fna_dir']+'/'+file)
			outfile = open(args['out_dir']+'/temp/genomes/'+file, 'w')
			for seq in Bio.SeqIO.parse(infile, 'fasta'):
				seq.id = seq.id.replace('-', '_')
				seq.seq = str(seq.seq).upper()
				outfile.write('>'+seq.id+'\n'+seq.seq+'\n')
		infile.close()
		outfile.close()
		args['fna_dir'] = args['out_dir']+'/temp/genomes'

	if 'rep_fna' in args and args['rep_fna'] is not None:
		infile = open(args['rep_fna'])
		outfile = open(args['out_dir']+'/temp/'+os.path.basename(args['rep_fna']), 'w')
		for seq in Bio.SeqIO.parse(infile, 'fasta'):
			seq.id = seq.id.replace('-', '_')
			seq.seq = str(seq.seq).upper()
			outfile.write('>'+seq.id+'\n'+seq.seq+'\n')
		infile.close()
		outfile.close()
		args['rep_fna'] = args['out_dir']+'/temp/'+os.path.basename(args['rep_fna'])

def locate_fpaths(args, in_dir, rep_fna=None, subset_list=None):
	subset_map = dict()

	for f in os.listdir(in_dir):
		subset_map[f] = 1

	if subset_list is not None:
		subset_map = dict()
		with open(subset_list, 'r') as fh:
			for ln in fh:
				subset_map[ln.rstrip()] = 1

	args["subset_map"] = subset_map

	ref_path = ""
	fpaths = []

	# Using the largest genome file in direcory for reference intead of randomly selecting anyone
	lg_fpath = ""
	cur_size = 0
	for f in os.listdir(in_dir):
		if f in subset_map:
			fpath = in_dir.rstrip('/')+'/'+f
			ftype = id_input_type(fpath)

			if os.path.isfile(fpath) and ftype == "fasta":
				fstats = os.stat(fpath)
				fpaths.append(fpath)
				if fstats.st_size >= cur_size:
					cur_size = fstats.st_size
					lg_fpath = fpath
			else:
				sys.stderr.write("skip {}: not fasta format\n".format(fpath))	

		else:
			sys.stderr.write("skip {}\n".format(f))

	if rep_fna is not None: # Using speficied reference genome
		ref_path = rep_fna
	else:
		ref_path = lg_fpath

	args['rep_fna'] = ref_path
	args['fna_paths'] = fpaths

def detect_single_chrom(ref_path):
	single_chrom = True
	chrom_cnt = 0
	with open(ref_path, 'r') as fh:
		for line in fh:
			if line[0] == '>':
				chrom_cnt = chrom_cnt + 1

				if chrom_cnt == 1:
					pass
				else:
					single_chrom = False
					break

	return single_chrom

def register_run_id(args, in_dir, single=False):
	args['run_id'] = in_dir.rstrip('/').split('/')[-1]

	if single is True:
		args['run_id'] = args['run_id'] + "_single"

	return args['run_id']

def register_msa_id(args, ref_path, fpaths):
	order_names = []

	for fpath in fpaths:
		order_names.append(fpath.rstrip('/').split('/')[-1])

	order_names.append(ref_path.rstrip('/').split('/')[-1])

	in_string = "".join(order_names)
	args['msa_id'] = hashlib.md5(in_string.encode()).hexdigest()

	return args['msa_id']

def auto_min_pid_by_delta(coords_path, idt_delta):
	min_pid_by_delta = 0

	# fields = [('s1',int),('e1',int),
	#		  ('s2',int),('e2',int),
	#		  ('len1',int),('len2',int),
	#		  ('pid',float),
	#		  ('c1',str),('c2',str)]

	pids = []
	with open(coords_path) as f:
		for i in range(5):
			next(f)
		for l in f:
			values = l.replace(' | ', ' ').split()
			pid = float(values[6])
			pids.append(pid)

	avg_pid = sum(pids)/len(pids)

	min_pid_by_delta = avg_pid * (1 - idt_delta)

	return min_pid_by_delta

def run_mummer4_single(fpath, genome_id, ref_fpath, rep_id, out_dir, skip_align, min_pid, min_aln_len, max_pid_delta, internal_thread_num):
	print("	%s - %s" % (rep_id, genome_id))

	try: os.makedirs(out_dir)
	except: pass

	log = open(out_dir+'/log','w')

	if skip_align is True and os.path.isfile("%s/%s.delta" % (out_dir, genome_id)):
		log.write('nucmer alignment was skipped\n')
		print('	nucmer alignment skipped\n')
	else:
		command = "nucmer "
		command += "-t %s " % internal_thread_num
		command += "%s " % ref_fpath
		command += "%s " % fpath
		command += "--prefix %s/%s " % (out_dir, genome_id)
		out, err = run_command(command)
		log.write(str(out)+'\n'+str(err))

	command = "delta-filter -q -r "
	command += "-i %s " % str(min_pid)
	command += "-l %s " % str(min_aln_len)
	command += "%s/%s.delta " % (out_dir, genome_id)
	command += "> %s/%s.filter.delta.1" % (out_dir, genome_id)
	out, err = run_command(command)
	log.write(str(out)+'\n'+str(err))

	command = "show-coords "
	command += "%s/%s.filter.delta.1 " % (out_dir, genome_id)
	command += "> %s/%s" % (out_dir, 'coords.tmp')
	out, err = run_command(command)
	log.write(str(out)+'\n'+str(err))

	coords_path = "{}/{}".format(out_dir, 'coords.tmp')
	min_pid_by_delta = auto_min_pid_by_delta(coords_path, max_pid_delta)

	command = "delta-filter -q -r "
	command += "-i %s " % str(min_pid_by_delta)
	command += "-l %s " % str(min_aln_len)
	command += "%s/%s.delta " % (out_dir, genome_id)
	command += "> %s/%s.filter.delta" % (out_dir, genome_id)
	out, err = run_command(command)

	for utility in ['coords', 'snps', 'diff']:
		command = "show-%s " % utility
		command += "%s/%s.filter.delta " % (out_dir, genome_id)
		command += "> %s/%s" % (out_dir, utility)
		out, err = run_command(command)
		log.write(str(out)+'\n'+str(err))


def run_mummer4(args):
	fpaths = args['fna_paths']
	if 'tag_genome_paths' in args:
		fpaths = args['tag_genome_paths']
		
	ref_fpath = args['rep_fna']
	if 'tag_ref' in args:
		ref_fpath = args['tag_ref']

	register_run_id(args, args['fna_dir'])
	register_msa_id(args, ref_fpath, fpaths)

	print("reference genome path: %s" % ref_fpath)

	args['mummer4_dir'] = args['out_dir']+'/temp/mummer4/'+args['run_id']
	try: os.makedirs(args['mummer4_dir'])
	except: pass

	shutil.copy(ref_fpath, os.path.join(args['mummer4_dir'], 'reference.fna'))

	arg_list = []
	rep_id = ref_fpath.split('/')[-1].replace('.fna', '')

	print("[paired alignment]: start")
	for fpath in fpaths:
		genome_id = fpath.split('/')[-1].replace('.fna', '')
		out_dir = '%s/aln/%s' % (args['mummer4_dir'], genome_id)

		arg_list.append([fpath, genome_id, ref_fpath, rep_id, out_dir, args['skip_align'], args['min_pid'], args['min_aln_len'], args['max_pid_delta'], 1])

	print("[paired alignment]: done")

	parallel(run_mummer4_single, arg_list, args['threads'])

	msa_path = gen_msa.build_msa(indir=args['mummer4_dir'], overwrite=True, subset=args["subset_map"])

	shutil.copy(os.path.join(args['mummer4_dir'], 'reference.fna'), args['out_dir'])

	args['msa_path'] = args['out_dir'] + '/tag_msa.fna'
	shutil.move(msa_path, args['msa_path'])

	args['msa_type'] = 'xmfa-mummer4'
	
	args['tag_list_path'] = args['out_dir'] + '/tag_paths.list'

	with open(args['tag_list_path'], 'w') as fh:
		for fpath in fpaths:
			fh.write("{}\n".format(fpath.rstrip()))
	

def run_mash_scketch(args):
	ref_fpath = args['rep_fna']
	fpaths = args['fna_paths']

	register_run_id(args, args['fna_dir'])
	register_msa_id(args, ref_fpath, fpaths)

	print("reference genome path: %s" % ref_fpath)

	args['mash_dir'] = args['out_dir']+'/temp/mash/'+args['run_id']

	try: os.makedirs(args['mash_dir'])
	except: pass

	args['fna_list_path'] = args['mash_dir'] + '/in_fna.list'

	with open(args['fna_list_path'], 'w') as fh:
		for fpath in fpaths:
			fh.write("{}\n".format(fpath))
	
	print("[building mash sketch]: start")

	command = "mash sketch "
	command += "-k %s " % str(args['sketch_k'])
	command += "-s %s " % str(args['sketch_size'])
	command += "-p %s " % str(args['threads'])
	command += "-o %s " % (args['mash_dir']+'/mash_sketch')
	command += "-l %s " % args['fna_list_path']

	out, err = run_command(command)
	sys.stderr.write(str(out)+'\n'+str(err))

	args['mash_sketch_path'] = args['mash_dir']+'/mash_sketch.msh'

def run_mash_dist(args):
	sketch_path = args['mash_sketch_path']

	assert os.path.exists(sketch_path)

	args['mash_dist_path'] = args['mash_dir'] + '/mash_dist.tsv'

	print("[calculating mash distance]: start")

	command = "mash dist "
	command += "-p %s " % str(args['threads'])
	command += "%s %s " % (sketch_path, sketch_path)
	command += "> %s " % args['mash_dist_path'] 

	out, err = run_command(command)
	sys.stderr.write(str(out)+'\n'+str(err))

def do_precut(args):
	dist_path = args['mash_dist_path']

	assert os.path.exists(dist_path)

	args['cut_dist_path'] = args['mash_dir'] + '/mash_dist.cut.tsv'

	print("[cut mash distance: {}]: start".format(str(args['precut'])))

	command = "awk '$3 < %s' " % str(args['precut'])
	command += "%s " % dist_path
	command += "> %s " % args['cut_dist_path']
	
	out, err = run_command(command)
	sys.stderr.write(str(out)+'\n'+str(err))

def id_clusters(args):
	run_mash_scketch(args)

	run_mash_dist(args)

	s_cut = args['start_cutoff']
	e_cut = args['end_cutoff']
	r_fac = args['range_factor']

	total_n = len(args['fna_paths'])

	maf = args['snp_freq']

	critical_n = math.ceil(1 / maf)
	
	do_precut(args)
	dist_path = args['cut_dist_path']
	assert os.path.exists(dist_path)

	optimal_clusters, optimal_d, optimal_n = [], None, None
	while s_cut <= args['precut']:
		optimal_clusters, optimal_d, optimal_n, firstcut_exit = id_genome_clusters.build_genome_blocks(dist_path, total_n, critical_n, s_cut, e_cut, r_fac)
		if firstcut_exit is True:
			s_cut = s_cut + 0.01
		else:
			break
		

	clust_genomes = dict()
	tag_genomes = []
	for cluster in optimal_clusters:
		tag_genomes.append(cluster.tag_genome)
		for genome in cluster.genomes: 
			clust_genomes[genome] = 1

	
	for fpath in args['fna_paths']:
		if fpath not in clust_genomes:
			tag_genomes.append(fpath)
	
	args['tag_genome_paths'] = tag_genomes

def id_tag_ref(args):
	if 'mash_dist_path' not in args or not os.path.exists(args['mash_dist_path']):
		run_mash_scketch(args)
		run_mash_dist(args)

	dist_path = args['mash_dist_path']

	tag_paths = args['fna_paths']
	if 'tag_genome_paths' in args and len(args['tag_genome_paths']) > 1:
		tag_paths = args['tag_genome_paths']

	centroid = id_centroid.identify(tag_paths, dist_path)
	
	print(centroid)

	args['tag_ref'] = centroid
	args['rep_fna'] = centroid

def run_kmerset_validate(args):
	assert os.path.exists(args['kmer_set'])
	assert os.path.exists(args['tag_list'])

	args['kmer_prof_path'] = args['out_dir']+'/kmer_prof.tsv'

	args['check_fna_paths'] = args['out_dir']+'/check_fna_paths.list'
	if 'fna_paths' in args:
		with open(args['check_fna_paths'], 'w') as fh:
			for fpath in args['fna_paths']:
				fh.write("{}\n".format(fpath))

	print("[validating kmer set]: start")

	command = "callm_db_val "
	command += "-d %s " % args['kmer_set']
	command += "-n %s " % args['genome_name']
	command += "-t %s " % args['threads']
	#command += "-L %s " % args['tag_list'] 
	command += "-L %s " % args['check_fna_paths']
	command += "-o %s " % args['kmer_prof_path']

	out, err = run_command(command)
	sys.stderr.write(str(out)+'\n'+str(err))

def filter_kmers(args):
	assert os.path.exists(args['kmer_prof_path'])

	args['filtered_kmer_path'] = args['out_dir']+'/selected_kmers.tsv'

	with open(args['filtered_kmer_path'], 'w') as fw:
		with open(args['kmer_prof_path'], 'r') as fh:
			for line in fh:
				items = line.rstrip().split('\t')

				nonsingle_hit = int(items[8])
				
				null_hit = int(items[6])
				single_hit = int(items[7])
				
				ref_hit = int(items[10])
				alt_hit = int(items[11])

				if nonsingle_hit > 0:
					continue

				if single_hit / (single_hit + null_hit) < 0.5:
					continue

				if ref_hit == 0 or alt_hit == 0:
					continue

				rec1 = "{}\t{}0{}".format(items[2], items[9], items[0])
				rec2 = "{}\t{}1{}".format(items[3], items[9], items[0])
				rec3 = "{}\t{}0{}".format(items[4], items[9], items[0])
				rec4 = "{}\t{}1{}".format(items[5], items[9], items[0])

				fw.write("{}\n{}\n{}\n{}\n".format(rec1, rec2, rec3, rec4))

def run_build_db(args):
	assert args['filtered_kmer_path']
	
	args['kmer_db_path'] = args['out_dir']+'/kmer_db.bin'

	command = "callm_db_build "
	command += "%s " % args['filtered_kmer_path']
	command += "> %s " % args['kmer_db_path']

	out, err = run_command(command)
	sys.stderr.write(str(out)+'\n'+str(err))

def read_input_dir(args, in_dir, subset_list=None):
	subset_map = dict()

	for f in os.listdir(in_dir):
		subset_map[f] = 1

	if subset_list is not None:
		subset_map = dict()
		with open(subset_list, 'r') as fh:
			for ln in fh:
				subset_map[ln.rstrip()] = 1

	args["subset_map"] = subset_map

	fna_paths = []
	fq_paths = []

	for f in os.listdir(in_dir):
		if f in subset_map:
			fpath = in_dir.rstrip('/')+'/'+f
			print(fpath)
			
			if os.path.isdir(fpath):
				continue

			assert os.path.isfile(fpath)
			ftype = id_input_type(fpath)

			if ftype == "unknown":
				sys.stderr.write("skip {}: unknown input type\n".format(fpath))	
			elif ftype == "not_supported":
				sys.stderr.write("skip {}: compressed fasta is not supported yet\n".format(fpath))
			elif ftype == "fasta":
				fna_paths.append(fpath)
			elif ftype in ["fastq", "fastq.gz", "fastq.lz4", "fastq.bz2"]:
				fq_paths.append(fpath)
			else:
				assert False
		else:
			sys.stderr.write("skip {}\n".format(f))

	fq_pairs = []
	if len(fq_paths) > 1:
		fq_pairs = pair_inputs(fq_paths)

	args['fna_paths'] = fna_paths
	args['fq_paths'] = fq_paths
	args['fq_pairs'] = fq_pairs


def id_input_type(fpath):
	in_type = "fastq" #default

	fn_its = fpath.split("/")[-1].split(".")

	fn_end = ""
	if fn_its[-1] in ['gz', 'lz4', 'bz2']:
		fn_end = fn_its[-2]
	else:
		fn_end = fn_its[-1]

	if fn_end in ['fa', 'fsa', 'fna', 'fasta']:
		in_type = "fasta"
	elif fn_end in ['fq', 'fastq']:
		in_type = "fastq"
	else:
		in_type = "unknown"

	if fn_its[-1] in ['gz', 'lz4', 'bz2']:
		if fn_end in ['fa', 'fsa', 'fna', 'fasta']:
			in_type = "not_supported"
		else:
			in_type = in_type + '.' + fn_its[-1]
			
	return in_type

def pair_inputs(fq_paths):
	pairs = dict()

	for fqpath in fq_paths:
		fn_its = fqpath.split("/")[-1].split(".")
		fq_name_parts = fn_its[0].split("_")
		
		if len(fq_name_parts) != 2:
			continue

		if fq_name_parts[1] not in ["1", "2"]:
			continue

		if fq_name_parts[0] not in pairs:
			pairs[fq_name_parts[0]] = dict()

		pairs[fq_name_parts[0]][fq_name_parts[1]] = fqpath
	
	real_pairs = []
	for name in pairs.keys():
		if "1" in pairs[name] and "2" in pairs[name]:
			real_pairs.append([pairs[name]["1"], pairs[name]["2"], name])
			
	return real_pairs

def genotype_single_genomes(args):
	ref_fpath = args['ref_genome']
	fpaths = args['fna_paths']

	print("reference genome path: %s" % ref_fpath)

	args['genotype_dir'] = args['out_dir']+'/temp/genotype'
	try: os.makedirs(args['genotype_dir'])
	except: pass

	args['gt_results_dir'] = args['out_dir']+'/gt_results'
	try: os.makedirs(args['gt_results_dir'])
	except: pass

	arg_list = []
	arg_list_gt = []
	rep_id = ref_fpath.split('/')[-1].replace('.fna', '')

	global ref 
	ref = read_ref(ref_fpath)

	global genos 
	genos = extract_genotypes(args['vcf'])

	print("[paired alignment]: start")
	for fpath in fpaths:
		genome_id = fpath.split('/')[-1]
		out_dir = '%s/aln/%s' % (args['genotype_dir'], genome_id)
		arg_list.append([fpath, genome_id, ref_fpath, rep_id, out_dir, False, args['min_pid'], args['min_aln_len'], args['max_pid_delta'], 1])

		coord_path = out_dir + '/coords'
		snp_path = out_dir + '/snps'
		output = args['gt_results_dir'] + '/' + genome_id + ".tsv"
		arg_list_gt.append([genos, ref, coord_path, snp_path, output])

	print("[paired alignment]: done")

	parallel(run_mummer4_single, arg_list, args['threads'])
	parallel(run_single_fasta_gt, arg_list_gt, args['threads'])

def read_ref(fpath):
	seq_recs = list(SeqIO.parse(fpath, "fasta"))

	rec_table = dict()
	for rec in seq_recs:
		rec_table[rec.id] = str(rec.seq).upper()

	return rec_table

def extract_genotypes(vcf_path):
	genos = []
	with open(vcf_path, 'r') as fh:
		for l in fh:
			if l[0] == "#":
				continue
			else:
				values = l.rstrip().split('\t')[:5]

				chrom = values[0]
				pos_r = int(values[1])
				gid = values[2]
				allele_ma = values[3]
				allele_mi = values[4]

				if len(allele_mi) > 1:
					continue

				genos.append([chrom, str(pos_r), gid, allele_ma, allele_mi])

	return genos 

def run_single_fasta_gt(genos, ref, coord_path, snp_path, output):
	coord_map = dict()

	with open(coord_path, 'r') as fh:
		for i in range(5):
			next(fh)
		for l in fh:
			values = l.replace(' | ', ' ').split()
			# position in coords file is 1 indexed compared to 0 indexed in vcf
			start = int(values[0]) - 1
			end = int(values[1]) - 1
			chrom = values[7]

			assert end > start

			if chrom not in coord_map:
				coord_map[chrom] = []

			coord_map[chrom].append([start, end])


	snp_map = dict()
	with open(snp_path) as fh:
		for i in range(5):
			next(fh)
		for l in fh:
			values = l.replace(' | ', ' ').split()
			# position in snps file is 1 indexed compared to 0 indexed in vcf
			pos_r = int(values[0]) - 1 
			allele_r = values[1]
			allele_a = values[2]
			chrom = values[10]
			
			if allele_r == "." or allele_a == ".":
				continue

			if chrom not in snp_map:
				snp_map[chrom] = dict()

			snp_map[chrom][pos_r] = [allele_r, allele_a]
	
	gtypes = []
	for geno in genos:
		chrom = geno[0]
		pos_r = int(geno[1])
		gid = geno[2]
		allele_ma = geno[3]
		allele_mi = geno[4]
		
		if chrom not in coord_map:
			continue

		for g_range in coord_map[chrom]:
			if pos_r >= g_range[0] and pos_r <= g_range[1]:
				if chrom in snp_map and pos_r in snp_map[chrom]:
					if allele_mi == snp_map[chrom][pos_r][1]:
						gtypes.append([chrom, str(pos_r), gid, allele_ma, allele_mi, '0', '1'])
					else:
						gtypes.append([chrom, str(pos_r), gid, allele_ma, allele_mi, '1', '0'])
				else:
					assert chrom in ref	
					allele_r = ref[chrom][pos_r]
					if allele_mi == allele_r:
						gtypes.append([chrom, str(pos_r), gid, allele_ma, allele_mi, '0', '1'])
					else:
						gtypes.append([chrom, str(pos_r), gid, allele_ma, allele_mi, '1', '0'])

	with open(output, 'w') as fw:
		for gtype in gtypes:
			fw.write("{}\n".format("\t".join(gtype)))

def genotype_reads(args):
	fpaths = args['fq_paths']

	args['genotype_dir'] = args['out_dir']+'/temp/genotype'
	try: os.makedirs(args['genotype_dir'])
	except: pass

	args['gt_results_dir'] = args['out_dir']+'/gt_results'
	try: os.makedirs(args['gt_results_dir'])
	except: pass

	gt_paths = []
	outname = '%s/iso_gt' % args['genotype_dir']
	try: os.makedirs(outname)
	except: pass

	mode = 2
	if args['mode'] == "very-fast":
		mode = 10
	elif args['mode'] == "fast":
		mode = 5
	elif args['mode'] == 'sensitive':
		mode = 2
	elif args['mode'] == 'very-sensitive':
		mode = 1
	else:
		assert False

	command = "iso_gt_mtar "
	command += "-d %s " % args['kmer_db_path']
	command += "-t %s " % args['threads']
	command += "-j %s " % mode
	command += "-o %s/" % outname
	command += "%{in} "
	command += "-f "

	for fpath in fpaths:
		command += "%s " % fpath
		gt_paths.append(outname + '/' + extract_fastq_path_name(fpath) + ".tsv")
	
	out, err = run_command(command)
	sys.stderr.write(str(out)+'\n'+str(err))

	merge_paths = []
	if args["merge_pairs"]:
		assert "fq_pairs" in args

		for fq_pair in args["fq_pairs"]:
			fq_1 = fq_pair[0]
			fq_2 = fq_pair[1]
			fq_name = fq_pair[2]

			fq_gt_1 = extract_fastq_path_name(fq_1) + ".tsv"
			fq_gt_2 = extract_fastq_path_name(fq_2) + ".tsv"
			
			fq_merge = dict()
			for fq_gt in [fq_gt_1, fq_gt_2]:
				with open(fq_gt, 'r') as fh:
					for line in fh:
						items = line.rstrip().split('\t')
						if items[0] not in fq_merge:
							fq_merge[items[0]] = int(items[0])
						else:
							fq_merge[items[0]] += int(items[0])

			merge_output = outname + "/" + fq_name + ".merged.tsv"
			with open(merge_output, 'w') as fw:
				for snp in fq_merge.keys():
					fw.write("{}\t{}\n".format(snp, str(fq_merge[snp])))

			merge_paths.append(merge_output)
	
	arg_list = []
	for gt_path in gt_paths + merge_paths:
		fq_id = '.'.join(gt_path.split('/')[-1].split('.')[:-1])
		output = args['gt_results_dir'] + '/' + fq_id + '.reads.tsv'
		arg_list.append([args['vcf'], gt_path, output])

	parallel(run_parse_single, arg_list, args['threads'])			

def extract_fastq_path_name(fpath):
	# chop off all leading '.' and '/'
	pparts = []
	real_idx = 0
	for i, ppart in enumerate(fpath.split('/')):
		if ppart == '.' or ppart == "..":
			continue
		else:
			real_idx = i
			break
	
	vpath = '/'.join(fpath.split('/')[real_idx:])

	path_parts = vpath.split('.')
	real_parts = []
	if path_parts[-1] in ['gz', 'lz4', 'bz2']:
		real_parts = path_parts[:-2]
	elif path_parts[-1] in ['fq', 'fastq']:
		real_parts = path_parts[:-1]
	else:
		assert False
	
	return ".".join(real_parts).replace('/', '_').replace('.','_')


def run_parse_single(vcf_path, gt_path, output):
	snp_map = dict()

	with open(gt_path, 'r') as fh:
		for line in fh:
			values = line.rstrip().split('\t')
			snp = values[0]
			count = values[1]

			allele_type = int(snp[6])
			assert allele_type in [0, 1]

			gid = snp[7:]
			
			if gid not in snp_map:
				snp_map[gid] = [0, 0]
			
			snp_map[gid][allele_type] = snp_map[gid][allele_type] + int(count)

	gtypes = []
	with open(vcf_path, 'r') as fh:
		for l in fh:
			if l[0] == "#":
				continue
			else:
				values = l.rstrip().split('\t')[:5]

				chrom = values[0]
				pos_r = int(values[1])
				gid = values[2]
				allele_ma = values[3]
				allele_mi = values[4]

				if len(allele_mi) > 1:
					continue

				if gid in snp_map: 	
					gtypes.append([chrom, str(pos_r), gid, allele_ma, allele_mi, str(snp_map[gid][0]), str(snp_map[gid][1])])
	
	with open(output, 'w') as fw:
		for gtype in gtypes:
			fw.write("{}\n".format("\t".join(gtype)))

def call_snps_main(args):
	cmdl_str = ' '.join(sys.argv[1:])

	if args['data_type'] in ['genomes', 'end_to_end']:
		locate_fpaths(args, args['fna_dir'], args['rep_fna'], args['subset_list'])


	if args['data_type'] in ['genomes', 'end_to_end']:
		if args["has_completeness"]:
			if args["completeness"]:
				args["min_prev"] = (1 - float(args["missing_ratio"])) * float(args["completeness"])
			elif args["completeness_list"]:
				completeness_map = {}
				with open(args["completeness_list"], 'w') as fh:
					for line in fh:
						items = line.rstrip('').split('\t')
						completeness_map[items[0]] = float(items[1])
				
				ref_fpath = args['rep_fna']
				fpaths = args['fna_paths']
				
				completenesses = []

				for fpath in fpaths:
					fname = fpath.rstrip('/').split('/')[-1]
					if fname in completeness_map:
						completenesses.append(completeness_map[fname])
					else:
						sys.exit("missing completeness: {}".format(fpath))

				avg_completeness = sum(completenesses)/len(completenesses)
				args["min_prev"] = (1 - float(args["missing_ratio"])) * avg_completeness
			else:
				print("useless option --has-completeness")
	
	if len(args['fna_paths']) <= math.ceil(1 / args['snp_freq']):
		print("[Warning] Total number of genomes ({}) < min. number of genomes required for effective SNP calling with MAF {} ({})".format(len(args['fna_paths']), args['snp_freq'], math.ceil(1 / args['snp_freq'])))
		print("[Warning] Skip tag genome selection, all genomes will be used")
		args['keep_redundancy'] = True

	if args['data_type'] in ['genomes', 'end_to_end']:
		if not args['keep_redundancy']:
			id_clusters(args)

		if args['skip_centroid']:
			assert args['rep_fna'] is not None
			assert os.path.exists(args['rep_fna'])
		else:
			id_tag_ref(args)


	# >>> 1. Generate multiple-genome-alignment or pileups

	# data type is genomes: use parsnp to perform multiple genome alignment
	start = time.time()
	if args['data_type'] in ['genomes', 'end_to_end']:
		print("Running mummer4; start")
		run_mummer4(args)
		#args['mummer4_dir'] = '/Users/jasonshi/Documents/zjshi_github/snpMLST/unit_test_raw/snps_from_genomes/Borrelia_burgdorferi_56121/temp/mummer4/54d64396-732c-42b0-8e88-3de63e8a665e/msa.fna'
		# msa_path = gen_msa.build_msa(indir=args['mummer4_dir'], max_genomes=1280)
		# args['msa_path'] = '/Users/jasonshi/Documents/zjshi_github/snpMLST/unit_test_raw/snps_from_genomes/Borrelia_burgdorferi_56121/temp/mummer4/54d64396-732c-42b0-8e88-3de63e8a665e/msa.fa'
		# args['msa_type'] = 'xmfa-mummer4'
		print("Running mummer4; done!")
	print("Elapsed time: {}".format(time.time()-start))


	# >>> 2. Parse multiple-genome-alignment or pileup and call SNPs

	# fetch generator to parse msa columns or mpileup sites
	start = time.time()
	print("Fetching file-type-specific parser; start")
	if args['data_type'] in ['genomes', 'end_to_end', 'msa']:
		from align_io import msa
		if args['mem']:
			site_assembly = msa.iter_parse(args['msa_path'], args['msa_type'], args['max_samples'])
		else:
			site_assembly = msa.monolithic_parse(args['msa_path'], args['msa_type'], args['max_samples'])

	print("Fetching file-type-specific parser; done")
	print("Elapsed time: {}".format(time.time()-start))


	# id core-genome coords and snps
	start = time.time()
	print("Identifying core-snps; start")
	print("max sites: {}".format(args['max_sites']))
	print("min prevalence: {}".format(args['min_prev']))
	print("min MAF: {}".format(args['snp_freq']))

	if args['mem']:
		align_assembs = align_assembly.call_snps_iter(site_assembly, args['max_sites'], args['min_prev'], args['snp_freq'])
	else:
		align_assembs = align_assembly.call_snps(site_assembly, args['max_sites'], args['min_prev'], args['snp_freq'])
	print("Identifying core-snps; done")
	print("Elapsed time: {}".format(time.time()-start))

	# sys.exit()

	single_chrom_rep = False

	if args['mem'] is True and args['rep_fna'] is not None:
		single_chrom_rep = detect_single_chrom(args['rep_fna'])

	# write output files
	start = time.time()
	print("Writing snps to VCF; start")
	if args['mem']:
		header_ready = False
		coords_buffer = []
		for align_assemb in align_assembs:
			if len(align_assemb.snps) > 0:
				if not header_ready:
					vcf_io.write_coords_header(coords_buffer, args['out_dir'])
					vcf_io.write_vcf_header(align_assemb.snps, args['out_dir'], cmdl_str)
					header_ready = True

				# vcf_io.write_genome(core_genome.consensus_genome, args['out_dir'])
				coords_buffer = coords_buffer + align_assemb.coords
				vcf_io.write_vcf(align_assemb.snps, args['out_dir'], single_chrom_rep)

		vcf_io.write_coords(vcf_io.merge_coords(coords_buffer), args['out_dir'])
		# vcf_io.write_coords(coords_buffer, args['out_dir'])
	else:
		vcf_io.write_coords_header(align_assembs.coords, args['out_dir'])
		vcf_io.write_vcf_header(align_assembs.snps, args['out_dir'], cmdl_str)
		vcf_io.write_coords(align_assembs.coords, args['out_dir'])
		# vcf_io.write_genome(core_genome.consensus_genome, args['out_dir'])
		vcf_io.write_vcf(align_assembs.snps, args['out_dir'])
	print("Writing snps to VCF; done!")
	print("Elapsed time: {}".format(time.time()-start))


def build_db_main(args):
	args['kmer_size'] = 31

	genome_path, vcf_path, coords_path, tag_list_path = args['ref_genome'], args['vcf'], args['coords'], args['tag_list']
	k_size, k_type = args['kmer_size'], args['kmer_type']

	if args['fna_dir'] is not None:
		locate_fpaths(args, args['fna_dir'])

	genome_seq = build_db.open_genome_seq(genome_path)
	#snps = build_db.open_vcf_file(vcf_path)

	coords = None
	if coords_path is not None:
		coords = build_db.read_coords(coords_path)

	snp_gb_pos, snp_alleles = build_db.open_vcf_file_local(vcf_path)
	#snp_gb_pos = [int(snp.ID) for snp in snps]
	#snp_alleles = [[str(snp.REF), str(snp.ALT[0])] for snp in snps]
	#snp_kmers = fetch_snp_kmers(genome_seq, snp_gb_pos, snp_alleles, k_size, k_type, coords)

	genome_seqs = build_db.load_msa(args['msa'])
	snp_kmers = build_db.fetch_all_from_msa(genome_seqs, genome_seq, snp_gb_pos, snp_alleles, k_size, coords)

	args['kmer_set'] = args['out_dir'] + '/nr_kmer_set.tsv'

	build_db.dump_tsv(snp_kmers, args['kmer_set'])

	run_kmerset_validate(args)

	filter_kmers(args)

	run_build_db(args)


def genotype_main(args):
	read_input_dir(args, args['in_dir'], args['subset_list'])
	
	try: os.makedirs(args['out_dir'])
	except: pass

	if len(args["fna_paths"]) > 0:
		genotype_single_genomes(args)
	
	if len(args["fq_paths"]) > 0:
		genotype_reads(args)

def tree_main(args):
	concat_alleles.concat_allele_tree(args)

def end2end_main(args):
	try: os.makedirs(args['out_dir'])
	except: pass

	args['fna_dir'] = args['in_dir']
	locate_fpaths(args, args['in_dir'], args['rep_fna'], args['subset_list'])
	call_snps_main(args)

	args['kmer_size'] = 31
	args['ref_genome'] = args['rep_fna']
	args['vcf'] = args['out_dir'].rstrip('/') + '/core_snps.vcf' 
	args['coords'] = args['out_dir'].rstrip('/') + '/coords.tsv'
	args['tag_list'] = args['out_dir'].rstrip('/') + '/tag_paths.list'
	args['msa'] = args['out_dir'].rstrip('/') + '/tag_msa.fna'

	build_db_main(args)

	read_input_dir(args, args['in_dir'], args['subset_list'])
	if len(args["fna_paths"]) > 0:
		genotype_single_genomes(args)
	
	if len(args["fq_paths"]) > 0:
		genotype_reads(args)

def main():
	args = parse_args()

	try: os.makedirs(args['out_dir'])
	except: pass

	if args['data_type'] == 'genomes':
		call_snps_main(args)
	elif args['data_type'] == 'db':
		build_db_main(args)
	elif args['data_type'] == 'genotype':
		genotype_main(args)
	elif args['data_type'] == 'tree':
		tree_main(args)
	elif args['data_type'] == 'end_to_end':
		end2end_main(args)
	else:
		sys.exit("\nError: invalid subcommand\nSupported subcommand: genomes, db, genotype, tree, end_to_end\n")

if __name__ == "__main__":
	main()
