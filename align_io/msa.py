
def parse_control(msa_path, msa_type, max_sample=float('inf')):
	if msa_type == 'xmfa-parsnp':
		from align_io.xmfa_parsnp_io import parse
	elif msa_type == 'xmfa-mummer4':
		from align_io.xmfa_mummer4_io import parse
	elif msa_type == 'xmfa-mauve':
		from align_io.xmfa_mauve_io import parse
	elif msa_type == 'maf-mugsy':
		from align_io.maf_io import parse
	else:
		import sys
		sys.exit("\nError: invalid value for --msa-format: %s\n" % msa_type)
	return parse(msa_path, max_sample)

def monolithic_parse(msa_path, msa_type, max_sample=float('inf')):
	return parse_control(msa_path, msa_type, max_sample)

def iter_parse_control(msa_path, msa_type, max_sample=float('inf')):
	if msa_type == 'xmfa-parsnp':
		from align_io.xmfa_parsnp_io import iter_parse
	elif msa_type == 'xmfa-mummer4':
		from align_io.xmfa_mummer4_io import iter_parse
	elif msa_type == 'xmfa-mauve':
		from align_io.xmfa_mauve_io import iter_parse
	elif msa_type == 'maf-mugsy':
		from align_io.maf_io import iter_parse
	else:
		import sys
		sys.exit("\nError: invalid value for --msa-format: %s\n" % msa_type)
	return iter_parse(msa_path, max_sample)

def iter_parse(msa_path, msa_type, max_sample=float('inf')):
	for align in iter_parse_control(msa_path, msa_type, max_sample):
		yield align

def iterate_cols(msa_path, msa_type):
	for align in parse(msa_path, msa_type):
		for column in align.fetch_columns():
			yield column
