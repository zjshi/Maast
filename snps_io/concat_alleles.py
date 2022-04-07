import sys, os, argparse
import operator
from time import time

def read_input(in_dir, subset_list=None):
    fpaths = []
    fnames = []
    
    subset_map = dict()

    for f in os.listdir(in_dir):
        subset_map[f] = 1

    if subset_list is not None:
        subset_map = dict()
        with open(subset_list, 'r') as fh:
            for ln in fh:
                items = ln.rstrip().split('\t')
                assert len(items) == 2
                fname = items[0].split('/')[-1]
                subset_map[fname] = 1
    
    for f in os.listdir(in_dir):
        if f in subset_map:
            fpath = in_dir.rstrip('/')+'/'+f

            if os.path.isfile(fpath):
                fstats = os.stat(fpath)
                if fstats.st_size >= 0:
                    fpaths.append(fpath)
                    fnames.append(f)
                else:
                    sys.stderr.write("skip {}: empty file\n".format(fpath))
            else:
                sys.stderr.write("skip {}: not exist\n".format(fpath))

        else:
            sys.stderr.write("skip {}\n".format(f))

    return fpaths, fnames

def read_msa(msa_in):
    valid_chars = dict()
    valid_chars['A'] = 'A'
    valid_chars['a'] = 'A'
    valid_chars['C'] = 'C'
    valid_chars['c'] = 'C'
    valid_chars['G'] = 'G'
    valid_chars['g'] = 'G'
    valid_chars['T'] = 'T'
    valid_chars['t'] = 'T'
    valid_chars['-'] = '-'

    alns = dict()
    cur_seq = ""
    cur_sample = ""
    cur_contig = ""
    with open(msa_in, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                items = line.rstrip().split(' ')
                if items[0] not in alns:
                    alns[items[0]] = dict()
                if items[1] not in alns[items[0]]:
                    alns[items[0]][items[1]] = ""
                cur_sample = items[0]
                cur_contig = items[1]
            elif line[0] == '=':
                pass
            else:
                elems = []
                for char in line.rstrip().split():
                    if char not in valid_chars:
                        elems.append(char)
                    else:
                        elems.append(valid_chars[char])
                alns[cur_sample][cur_contig] = "".join(elems)

    aln_recs = []
    for sample in alns.keys():
        for contig in alns.keys():
            aln_recs.append([sample, contig, alns[cur_sample][cur_contig]])

    sorted_alns = sorted(aln_recs, key = lambda x: (x[1], x[0]))
    concat_alns = dict()
    for aln_rec in sorted_alns:
        if aln_rec[0] not in concat_alns:
            concat_alns[aln_rec[0]] = ""
        else:
            concat_alns[aln_rec[0]] += aln_rec[2]

    return concat_alns 

def read_aln(aln_in):
    valid_chars = dict()
    valid_chars['A'] = 'A'
    valid_chars['a'] = 'A'
    valid_chars['C'] = 'C'
    valid_chars['c'] = 'C'
    valid_chars['G'] = 'G'
    valid_chars['g'] = 'G'
    valid_chars['T'] = 'T'
    valid_chars['t'] = 'T'
    valid_chars['-'] = '-'

    alns = dict()
    cur_seq = ""
    cur_sample = ""
    with open(aln_in, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                cur_sample = line.rstrip()
                if cur_sample not in alns:
                    alns[cur_sample] = ""
            else:
                elems = []
                for char in line.rstrip().split():
                    if char not in valid_chars:
                        elems.append(char)
                    else:
                        elems.append(valid_chars[char])
                alns[cur_sample] += "".join(elems)

    return alns 

def write_aln(alns, out_path, max_gap=0.2):
    with open(out_path, 'w') as fh:
        for aln_key in alns.keys():
            n_gaps = 0
            total_len = len(alns[aln_key])
            for base in alns[aln_key]:
                if base == '-':
                    n_gaps += 1

            if n_gaps/total_len > max_gap:
                print("{}: skip {}".format(n_gaps/total_len, aln_key))
            else:
                fh.write("{}\n{}\n".format(aln_key, alns[aln_key]))

def read_gtp(input_path, min_depth):
    input_recs = dict()

    with open(input_path, 'r') as fh:
        for line in fh:
            items = line.rstrip().split("\t")

            contig_id = items[0]
            contig_pos = items[1]
            snp_key = contig_id + "__" + contig_pos

            cnt_allele_1 = int(items[5])
            cnt_allele_2 = int(items[6])
            if cnt_allele_1 + cnt_allele_2 < min_depth:
                continue

            allele = ""
            if cnt_allele_1 > cnt_allele_2:
                allele = items[3]
            else:
                allele = items[4]
            input_recs[snp_key] = allele 

    return input_recs

def union_inputs(inputs, names):
    first_union_in = dict()

    for input_recs in inputs:
        for snp_key in input_recs:
            if snp_key not in first_union_in:
                first_union_in[snp_key] = 1
    print("first_union_in: {}".format(len(first_union_in.keys())))

    union_in = dict()
    n_samples = len(names)
    for snp_key in first_union_in.keys():
        n_prev = 0
        allele_col = dict()
        for input_recs in inputs:
            if snp_key in input_recs:
                n_prev += 1
                allele_col[input_recs[snp_key]] = 1
        #if n_prev / n_samples >= 10 and len(allele_col.keys()) > 1:
        if len(allele_col.keys()) > 1:
            union_in[snp_key] = 1
        
    print("union_in: {}".format(len(union_in.keys())))
    all_keys = [[key.split('__')[0], int(key.split('__')[1])] for key in union_in.keys()]
    sorted_keys = sorted(all_keys, key = lambda x: (x[0], x[1]))
    #sorted_keys = sorted(all_keys, key = operator.itemgetter(0, 1))

    allele_aln = dict()
    for i, input_recs in enumerate(inputs):
        alleles = []
        for elem in sorted_keys:
            snp_key = elem[0] + "__" + str(elem[1])
            if snp_key in input_recs:
                alleles.append(input_recs[snp_key])
            else:
                alleles.append('-')
        allele_aln[names[i]] = alleles

    return allele_aln 

def concat_snps(allele_aln, allele_aln_fasta, max_gap, min_prev, min_maf, min_mac):
    with open(allele_aln_fasta, 'w') as fh:
        good_names = []
        for name in allele_aln.keys():
            n_gaps = 0
            total_len = len(allele_aln[name])
            for base in allele_aln[name]:
                if base == '-':
                    n_gaps += 1

            if n_gaps/total_len > max_gap:
                print("{}: skip {}".format(n_gaps/total_len, name))
            else:
                good_names.append(name)

        print(good_names)
        comm_aln = dict()
        comm_inds = []
        n_samples = len(good_names)

        for i, allele in enumerate(allele_aln[good_names[0]]):
            n_gaps = 0
            alleles = []
            allele_track = dict()
            for name in good_names:
                alleles.append(allele_aln[name][i])
                if allele_aln[name][i] == '-':
                    n_gaps += 1
                else:
                    if allele_aln[name][i] not in allele_track:
                        allele_track[allele_aln[name][i]] = 1
                    else:
                        allele_track[allele_aln[name][i]] += 1

            if (1 - n_gaps/n_samples) < min_prev:
                print("low prevalence: {}: skip {}".format(1 - n_gaps/n_samples, i))
                continue
            elif len(allele_track.keys()) <= 1:
                print("not a SNP site: skip {}".format(i))
                continue
            else:
                sorted_allele_track = sorted(allele_track.items(), key=lambda item: item[1], reverse=True)
                major_count = sorted_allele_track[0][1]
                minor_count = sorted_allele_track[1][1]
                if minor_count/n_samples < min_maf:
                    print("low min MAF: {}: skip {}".format(minor_count/n_samples, i))
                    continue
                elif minor_count < min_mac:
                    print("low min MAC: {}: skip {}".format(minor_count, i))
                    continue
                else:
                    comm_inds.append(i)

        print("number of good sites: {}".format(len(comm_inds)))
        print("number of good samples: {}".format(len(good_names)))
        
        for name in good_names:
            if name not in comm_aln:
                comm_aln[name] = []
            for i in comm_inds:
                comm_aln[name].append(allele_aln[name][i])
            fh.write(">{}\n{}\n".format(name, "".join(comm_aln[name])))

    return allele_aln_fasta


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

def run_fasttree(snp_str_fasta, out_dir):
    sys.stderr.write("[start] inferring max. likelihood tree\n")
    sys.stderr.write("\tsnp string fasta path: {}\n".format(snp_str_fasta))

    o_mat_path = out_dir + "/concat_allele.aln.mat"

    command = "FastTreeMP -makematrix -nt -gtr < "
    command += snp_str_fasta
    command += " > "
    command += o_mat_path

    environ = os.environ.copy()
    run_command(command, environ)
    sys.stderr.write("\tfinishing up, distance matrix is writtedn to {}\n".format(o_mat_path))

    o_tre_path = out_dir + "/concat_allele.aln.tre"

    command = "FastTreeMP -nt -gtr < "
    command += snp_str_fasta
    command += " > "
    command += o_tre_path

    environ = os.environ.copy()
    run_command(command, environ)

    sys.stderr.write("\tfinishing up, tree is writtedn to {}\n".format(o_tre_path))
    sys.stderr.write("[done] inferring max. likelihood tree\n")

def concat_allele_tree(args):
    in_dir = args['input_dir']
    in_path = args['input_list']

    out_dir = args['out_dir'].rstrip('/')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    min_sites_per_sample = args['min_sites_per_sample']
    max_gap_ratio = args['max_gap_ratio']
    min_site_prev = args['min_site_prev']
    min_maf = args['min_maf']
    min_mac = args['min_mac']

    paths, names = read_input(in_dir, in_path) 
    if len(names) != len(set(names)):
        sys.stderr.write("\n[error] names of input files are not unqiue.\n")
        sys.exit()
    input_recs = []
    aln_fasta = out_dir + '/concat_allele.aln.fasta'
    nonempty_names = []
    for i, path in enumerate(paths):
        if not os.path.exists(path):
            print("Skip input: {} does not exists.".format(path))
            continue
        if args["min_depth"] is not None:
            min_depth = args["min_depth"]
        input_rec = read_gtp(path, min_depth)
        if len(input_rec.keys()) < min_sites_per_sample:
            print("{}: skipped {}".format(len(input_rec.keys()), path))
        else:
            input_recs.append(input_rec)
            nonempty_names.append(names[i])
    allele_aln = union_inputs(input_recs, nonempty_names)
    concat_snps(allele_aln, aln_fasta, max_gap_ratio, min_site_prev, min_maf, min_mac)
    run_fasttree(aln_fasta, out_dir)
