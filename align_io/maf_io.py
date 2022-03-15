

class Alignment:
	def __init__(self, line):
		self.desc = line
		self.seqs = []

class Sequence:
	def __init__(self, line):
		values = line.rstrip().split()
		self.chrom = values[1]
		self.start = values[2]
		self.end = int(values[5])
		self.length = int(values[3])
		self.strand = values[4]

		self.seq = values[6].upper()

def parse(fpath):
	with open(fpath) as file:
		for line in file:
			if line[0] == '#': continue
			else: break
		alignment = Alignment(line)
		for line in file:
			if line[0] == 'a':
				yield alignment
				alignment = Alignment(line)
			elif line[0] == 's':
				sequence = Sequence(line)
				alignment.seqs.append(sequence)
	yield alignment

def iter_parse(fpath):
	with open(fpath) as file:
		for line in file:
			if line[0] == '#': continue
			else: break
		alignment = Alignment(line)
		for line in file:
			if line[0] == 'a':
				yield alignment
				alignment = Alignment(line)
			elif line[0] == 's':
				sequence = Sequence(line)
				alignment.seqs.append(sequence)
	yield alignment
