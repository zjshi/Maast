import os, time, sys
import numpy as np

def parse_seqs(path):
	with open(path) as file:
		try: id = next(file).split()[0].lstrip('>')
		except: return
		seq = ''
		for line in file:
			if line[0]=='>':
				yield id, seq
				try: id = line.split()[0].lstrip('>')
				except: return
				seq = ''
			else:
				seq += line.rstrip()
		yield id, seq

def parse_coords(fpath):
	fields = [('s1',int),('e1',int),
			  ('s2',int),('e2',int),
			  ('len1',int),('len2',int),
			  ('pid',float),
			  ('c1',str),('c2',str)]
	with open(fpath) as f:
		for i in range(5):
			next(f)
		for l in f:
			values = l.replace(' | ', ' ').split()
			yield dict([(f[0],f[1](v)) for f,v in zip(fields, values)])

def parse_snps(fpath):
	fields = [('p1',int),('b1',str),('b2',str),('p2',int),
			  ('buf',int),('dist',int),
			  ('r',int),('q',int),
			  ('s1',int),('s2',int),
			  ('c1',str),('c2',str)]
	with open(fpath) as f:
		for i in range(5):
			next(f)
		for l in f:
			values = l.replace(' | ', ' ').split()
			yield dict([(f[0],f[1](v)) for f,v in zip(fields, values)])

def build_msa(indir, overwrite=True, max_genomes=None, max_sites=None, msa_id=None, subset=None):
	start = time.time()

	aln_dir = os.path.join(indir, 'aln')

	if not os.path.exists(indir):
		sys.exit("Error: dir does not exist: %s" % indir)

	print("Reading reference genome")
	ref = {}
	chroms = []
	local_pos = np.array([])
	for id, seq in parse_seqs(os.path.join(indir, 'reference.fna')):
		chroms.append(id)
		ref[id] = np.array(list(seq.upper()))
		if len(local_pos) == 0:
			local_pos = np.arange(len(seq))
		else:
			local_pos = np.concatenate([local_pos, np.arange(len(seq))])
	print("   count contigs: %s" % len(ref))
	print("   count sites: %s" % sum([len(_) for _ in ref.values()]))

	print("Initializing alignments")
	genome_ids = os.listdir(aln_dir)

	if max_genomes is not None:
		genome_ids = genome_ids[:max_genomes]

	if len(subset) == 0:
		pass
	else:
		genome_ids = []
		for genome_id in os.listdir(aln_dir):
			if "{}.fna".format(genome_id) in subset:
				genome_ids.append(genome_id)

	print("   count genomes: %s" % len(genome_ids))
	genomes = {}
	for genome_id in genome_ids:
		genomes[genome_id] = {}
		for id, seq in ref.items():
			genomes[genome_id][id] = np.array(['-']*len(seq))

	print("Reading alignment blocks")
	for genome_id in genome_ids:
		fpath = '%s/%s/coords' % (aln_dir, genome_id)
		aln_length = 0
		for r in parse_coords(fpath):
			aln_length += (r['e1'] - r['s1'])
			genomes[genome_id][r['c1']][r['s1']-1:r['e1']] = ref[r['c1']][r['s1']-1:r['e1']]

	print("Reading SNPs")
	for genome_id in genome_ids:
		fpath = '%s/%s/snps' % (aln_dir, genome_id)
		for r in parse_snps(fpath):
			if r['b1'] == '.':
				continue
			elif r['b2'] == '.':
				genomes[genome_id][r['c1']][r['p1']-1] = '-'
			else:
				genomes[genome_id][r['c1']][r['p1']-1] = r['b2']

	chrom_aligns = {}
	for chrom in chroms:
		chrom_aligns[chrom] = ''
		for genome_id in genomes:
			chrom_aligns[chrom] = chrom_aligns[chrom] + '>{} {}\n{}\n'.format(genome_id, chrom, ''.join(genomes[genome_id][chrom]))

		chrom_aligns[chrom] = chrom_aligns[chrom] + '=\n'

	print("Writing fasta")

	fname = "msa.fa"

	if msa_id is not None:
		fname = "{}.fa".format(msa_id)

	msa_path = os.path.join(indir, fname)

	if overwrite is True:
		pass
	else:
		indx = 1
		while (os.path.isfile(msa_path)):
			fname = "msa.{}.fa".format(indx)
			msa_path = os.path.join(indir, fname)

	print("   path: %s" % msa_path)
	with open(msa_path, 'w') as f:
		for chrom in chroms:
			f.write(chrom_aligns[chrom])

	print("\nDone!")
	print("Time (s):", round(time.time()-start,2))

	return msa_path
"""
	print("Writing fasta")
	msa_path = os.path.join(indir, 'msa.fa')
	print("   path: %s" % msa_path)
	with open(msa_path, 'w') as f:
		for chrom in chroms:
			for genome_id in genomes:
				f.write('>%s %s\n' % (genome_id, chrom))
				f.write(''.join(genomes[genome_id][chrom])+'\n')
			f.write("=\n")

	print("\nDone!")
	print("Time (s):", round(time.time()-start,2))
"""
