"""
Program to classify the snps based on the following:
1. whether are at stop codons
2. where are located on a functional gene sequence
4. whether are synonymous (S) or non-synonymous (N)
5. whether are fourfold degenerate sites (S) or non-degenerate sites (N)
and count the number the number for each of the class,
and will ulimately estimate the following:
1. Pi: Nucleotide diversity: 2*p*(1-p) where p is the alternative allele frequency
2. Pi(N)/Pi(S): Pi(non-degenerate sites) / Pi(fourfold degenerate sites) ???
3. pN/pS: the ratio of non-synonymous to synonymous polymorphism rates
"""

from __future__ import division

import re, vcf, copy, argparse
import numpy as np

from Bio import SeqIO
from gff3 import Gff3

import pnps_dnds

def parse_args():
    """ Return dictionary of command line arguments
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS)

    parser.add_argument('program', help=argparse.SUPPRESS)

    parser.add_argument('--name', type=str, dest='name', required=True,
        help="""Name or identifier of the project""")
    parser.add_argument('--fna', type=str, dest='fna', required=True,
        help="""Path to consensus genome sequence""")
    parser.add_argument('--vcf', type=str, dest='vcf', required=True,
        help="""Path to core-genome SNPs""")
    parser.add_argument('--coords', type=str, dest='coords', default=None,
        help="""Path to core-genome coordinates""")
    parser.add_argument('--gff', type=str, dest='gff', default=None,
        help="""Path to gff file defining the CDS coordinates""")
    parser.add_argument('--with-mags', action='store_true', default=False,
        help="""The genome sequences or samples can be separated into MAGs and others""")
    parser.add_argument('--min-prev', type=float, dest='min_prev', default=None,
        help="""Minimal prevalence for SNP calling""")
    parser.add_argument('--with-header', action='store_true', default=False,
        help="""Output the header for the parameters""")
    parser.add_argument('--out', type=str, dest='out', default='/dev/stdout',
        help="""Path to output file (/dev/stdout)""")

    return vars(parser.parse_args())

def open_vcf_file(fpath, min_prev=None):
    """
    * ``Record.CHROM``; string
    * ``Record.POS``; int
    * ``Record.ID``; None
    * ``Record.REF``; string
    * ``Record.ALT``; list
    * ``Record.QUAL``; None
    * ``Record.FILTER``; list
    * ``Record.INFO``; dictionary

    additional attributes:
    * ``Record.FORMAT``; string
    * ``Record.samples``; list
    * ``Record.genotype``; object
    """
    vcf_reader = vcf.Reader(open(fpath, 'r'))

    raw_snps = [snp for snp in vcf_reader]

    sample_names = [sample.sample for sample in raw_snps[0].samples]
    n_sample = len(sample_names)

    vcf_snps = []
    for snp in raw_snps:
        if min_prev is not None:
            prev = snp.INFO['NS']/n_sample
            if prev > min_prev:
                vcf_snps.append(snp)
        else:
            vcf_snps.append(snp)

    return vcf_snps

def snp_stats(snps):
    n_snps = len(snps)

    bi_snps, tri_snps, quad_snps = classify_snps(snps)
    n_snps = len(snps)
    n_bi_snps = len(bi_snps)
    n_tri_snps = len(tri_snps)
    n_quad_snps = len(quad_snps)

    sample_names = [sample.sample for sample in snps[0].samples]
    n_sample = len(sample_names)

    values = [n_sample, n_snps, n_bi_snps, n_tri_snps, n_quad_snps]
    header = ['n_sample', 'n_snps', 'n_bi_snps', 'n_tri_snps', 'n_quad_snps']

    return values, header

def mag_snp_stats(snps):
    n_mag_snps = 0
    are_mags = ['_' in sample.sample for sample in snps[0].samples]
    are_refs = ['_' not in sample.sample for sample in snps[0].samples]

    if np.all(are_mags):
        all_mags = 1
        n_mag_snps = len(snps)
    else:
        all_mags = 0

        ref_inds = []
        for i, is_ref in enumerate(are_refs):
            if is_ref:
                ref_inds.append(i)

        static_sample = snps[0].samples[ref_inds[0]]

        for snp in snps:
            flgs = []
            for ind in ref_inds[1:]:
                sample = snp.samples[ind]
                flgs.append(static_sample.data.GP1 == sample.data.GP1)
                flgs.append(static_sample.data.GP2 == sample.data.GP2)
                flgs.append(static_sample.data.GP3 == sample.data.GP3)
                flgs.append(static_sample.data.GP4 == sample.data.GP4)

                if not np.all(flgs):
                    break

            if np.all(flgs):
                n_mag_snps = n_mag_snps + 1
            else:
                pass

    values = [n_mag_snps, all_mags]
    header = ['n_mag_snps', 'all_mags']

    return values, header

def core_genome_stats(coords):
    n_cg_bloc = 0

    cg_bloc_sizes = []
    avg_cg_bloc_size = 0
    mid_cg_bloc_size = 0

    cg_bloc_dists = []
    avg_cg_bloc_dist = 0
    mid_cg_bloc_dist = 0

    last_end = 0
    with open(coords, "r") as fh:
        next(fh)
        for line in fh:
            items = line.rstrip().split("\t")
            cg_bloc_size = int(items[2]) - int(items[1])+1
            cg_bloc_sizes.append(cg_bloc_size)

            if last_end > 0:
                cg_bloc_dist = int(items[1]) - last_end
                cg_bloc_dists.append(cg_bloc_dist)

            last_end = int(items[2])

    n_cg_bloc = len(cg_bloc_sizes)
    cg_size = sum(cg_bloc_sizes)
    avg_cg_bloc_size = cg_size / n_cg_bloc
    mid_cg_bloc_size = sorted(cg_bloc_sizes)[int((n_cg_bloc-1)/2)]

    n_cg_dist = len(cg_bloc_dists)
    avg_cg_bloc_dist = sum(cg_bloc_dists)/n_cg_dist
    mid_cg_bloc_dist = sorted(cg_bloc_dists)[int((n_cg_dist-1)/2)]

    values = [cg_size, n_cg_bloc, avg_cg_bloc_size, mid_cg_bloc_size, avg_cg_bloc_dist, mid_cg_bloc_dist]
    header = ['cg_size', 'n_cg_bloc', 'avg_cg_bloc_size', 'mid_cg_bloc_size', 'avg_cg_bloc_dist', 'mid_cg_bloc_dist']

    return values, header

def classify_snps(snps):
    bi_snps = []
    tri_snps = []
    quad_snps = []

    for snp in snps:
        if len(snp.ALT) == 1:
            bi_snps.append(snp)
        elif len(snp.ALT) == 2:
            tri_snps.append(snp)
        elif len(snp.ALT) == 3:
            quad_snps.append(snp)
        else:
            assert False

    return bi_snps, tri_snps, quad_snps

def genic_stats(seq_recs, rec_table, cds_recs, snps):
    genic_region, genic_masks = get_genic_region(rec_table, cds_recs)

    genic_region_size = len(genic_region)

    bi_snps, tri_snps, quad_snps = classify_snps(snps)

    alt_rec_tb = get_alt_seq_recs(rec_table, bi_snps)

    alt_genic_region, alt_genic_masks = get_genic_region(alt_rec_tb, cds_recs)

    sample_names = [sample.sample for sample in snps[0].samples]
    n_sample = len(sample_names)

    genic_array = [ genic_masks[snp.CHROM][snp.POS] for snp in snps]
    genic_snp_counts = sum(genic_array)

    bi_genic_array = [ genic_masks[snp.CHROM][snp.POS] for snp in bi_snps]
    bi_genic_snp_counts = sum(bi_genic_array)

    # print genic_region
    # print len(genic_region)
    # print alt_genic_region
    # print len(alt_genic_region)

    pn, ps, dn, ds, pn2ps, dn2ds = pnps_dnds.pnps_dnds(genic_region, alt_genic_region)

    snps2codons = map_codons(rec_table, cds_recs, bi_snps)

    n_syn_snp = 0
    n_nonsyn_snp = 0
    n_premature_stop_snp = 0
    n_stop_disrupt_snp = 0

    pi_all = 0
    pi_syn = 0
    pi_nonsyn = 0
    pi_4f_deg = 0
    pi_non_deg = 0

    for snp in bi_snps:
        pi_val = 2*float(snp.INFO['AF'][0])*(1-float(snp.INFO['AF'][0]))
        if snp.ID in snps2codons:
            codon_pair = snps2codons[snp.ID]
            if len(codon_pair) > 0:
                if not pnps_dnds.is_synonymous(codon_pair[0], codon_pair[1]):
                    pi_nonsyn = pi_nonsyn + pi_val
                    n_nonsyn_snp = n_nonsyn_snp + 1
                    if pnps_dnds.is_stop_codon(codon_pair[0]):
                        n_stop_disrupt_snp = n_stop_disrupt_snp + 1
                    elif pnps_dnds.is_stop_codon(codon_pair[1]):
                        n_premature_stop_snp = n_premature_stop_snp + 1
                    else:
                        pass
                else:
                    pi_syn = pi_syn + pi_val
                    n_syn_snp = n_syn_snp + 1

                if pnps_dnds.is_4f_deg(codon_pair[0], codon_pair[2]):
                    pi_4f_deg = pi_4f_deg + pi_val

                if pnps_dnds.is_non_deg(codon_pair[0], codon_pair[2]):
                    pi_non_deg = pi_non_deg + pi_val

        pi_all = pi_all + pi_val

    values = [genic_region_size, genic_snp_counts, bi_genic_snp_counts, pn2ps, pn, ps, dn2ds, dn, ds, n_syn_snp, n_nonsyn_snp, n_premature_stop_snp, n_stop_disrupt_snp, pi_all, pi_syn, pi_nonsyn, pi_4f_deg, pi_non_deg]
    header = ['genic_region_size', 'genic_snp_counts', 'bi_genic_snp_counts', 'pn2ps', 'pn', 'ps',  'dn2ds', 'dn', 'ds', 'n_syn_snp', 'n_nonsyn_snp', 'n_premature_stop_snp', 'n_stop_disrupt_snp', 'pi', 'pi_syn', 'pi_nonsyn', 'pi_4f_deg', 'pi_non_deg']

    return values, header

def get_alt_seq_recs(seq_rec_tb, snps):
    alt_rec_tb = copy.deepcopy(seq_rec_tb)

    for snp in snps:
        if seq_rec_tb[snp.CHROM][snp.POS] != snp.REF:
            alt_rec_tb[snp.CHROM] = str(alt_rec_tb[snp.CHROM][:snp.POS]) + str(snp.REF) + str(alt_rec_tb[snp.CHROM][snp.POS+1:])
        elif seq_rec_tb[snp.CHROM][snp.POS] != snp.ALT[0]:
            alt_rec_tb[snp.CHROM] = str(alt_rec_tb[snp.CHROM][:snp.POS]) + str(snp.ALT[0]) + str(alt_rec_tb[snp.CHROM][snp.POS+1:])
        else:
            print "{}: {} - {}: {}, {}, {}, {}, {}".format(
                    snp.ID, snp.REF, snp.ALT[0], seq_rec_tb[snp.CHROM][snp.POS-2], seq_rec_tb[snp.CHROM][snp.POS-1],
                    seq_rec_tb[snp.CHROM][snp.POS], seq_rec_tb[snp.CHROM][snp.POS+1], seq_rec_tb[snp.CHROM][snp.POS+2])

    return alt_rec_tb

### map codons to snps
def map_codons(seq_rec_tb, cds_recs, snps):
    codons = dict()
    for cds_rec in cds_recs:
        sid = cds_rec["seqid"]
        if sid not in codons:
            codons[sid] = dict()

        genic_region = seq_rec_tb[sid][cds_rec["start"]-1:cds_rec["end"]]

        for i in range(cds_rec["start"]-1, cds_rec["end"], 3):
            start_pos = i - cds_rec["start"] + 1
            codon = genic_region[start_pos:start_pos+3]
            codons[sid][i] = "0.{}".format(codon)
            codons[sid][i+1] = "1.{}".format(codon)
            codons[sid][i+2] = "2.{}".format(codon)

    snps2codons = dict()

    n_nonsyn = 0
    for snp in snps:
        sid = snp.CHROM
        snps2codons[snp.ID] = []
        if sid in codons:
            pos = snp.POS
            if pos in codons[sid]:
                codon_code = codons[sid][pos]
                i, codon = codon_code.split(".")
                i = int(i)
                alt_codon = codon[:i]+str(snp.ALT[0])+codon[i+1:]

                if codon == alt_codon:
                    alt_codon = codon[:i]+str(snp.REF)+codon[i+1:]

                snps2codons[snp.ID].append(codon)
                snps2codons[snp.ID].append(alt_codon)
                snps2codons[snp.ID].append(i)
                n_nonsyn = n_nonsyn + 1
                # print [i, codon, alt_codon]
    return snps2codons


### get genic masks
def get_genic_region(seq_rec_tb, cds_recs):
    genic_region = ""
    genic_masks = dict()

    for sid, seq in seq_rec_tb.iteritems():
        genic_mask = np.repeat(False, len(seq_rec_tb[sid]))
        genic_masks[sid] = genic_mask

    for cds_rec in cds_recs:
        sid = cds_rec["seqid"]
        genic_masks[sid][cds_rec["start"]-1:cds_rec["end"]] = True
        genic_region = genic_region + seq_rec_tb[sid][cds_rec["start"]-1:cds_rec["end"]]

    return genic_region, genic_masks

### get pi genic diversity
def get_pi(snps):
    all_pis = []

    for snp in snps:
        pi = 2*snps.INFO['AF']*(1-snps.INFO['AF'])
        all_pis.append(pi)

    return all_pis

def open_gff_file(gff):
    gff_hd = Gff3(gff)

    cds_recs = []
    for line in gff_hd.lines[4:]:
        if 'seqid' in line and 'type' in line and 'strand' in line:
            if line['type'] == 'CDS' and line['strand'] == '+':
                cds_recs.append(line)

    print "number of cds: {}".format(len(cds_recs))
    return cds_recs

def get_gff_summary(cds_recs):
    n_cds = 0

    cds_sizes = []
    avg_cds_size = 0
    mid_cds_size = 0

    n_cds = len(cds_recs)

    for rec in cds_recs:
        cds_size = rec["end"] - rec["start"] + 1
        cds_sizes.append(cds_size)
        if (cds_size % 3) != 0:
            print "odd cds with starting pos at: {}".format(rec["start"])

    total_cds_size = sum(cds_sizes)
    avg_cds_size = total_cds_size / n_cds
    mid_cds_size = sorted(cds_sizes)[int((n_cds-1)/2)]

    return [n_cds, avg_cds_size, mid_cds_size]

def get_seq_recs(fna):
    seq_recs = list(SeqIO.parse(fna, "fasta"))

    rec_table = dict()
    for rec in seq_recs:
        rec_table[rec.id] = str(rec.seq).upper()

    return seq_recs, rec_table

def get_ref_genome_size(seq_recs):
    rg_size = 0

    for rec in seq_recs:
        rg_size = rg_size + len(rec.seq)

    return rg_size

def get_stats(args):
    seq_recs, rec_table = get_seq_recs(args['fna'])
    rg_size = get_ref_genome_size(seq_recs)

    header = ['name', 'rg_size']
    params = [args['name'], rg_size]

    snps = open_vcf_file(args['vcf'], args['min_prev'])
    snp_values, snp_header = snp_stats(snps)
    params = params + snp_values
    header = header + snp_header

    if args['with_mags']:
        mag_snp_values, mag_snp_header = mag_snp_stats(snps)
        params = params + mag_snp_values
        header = header + mag_snp_header

    if args['coords'] is not None:
        core_genome_values, core_genome_header = core_genome_stats(args['coords'])
        params = params + core_genome_values
        header = header + core_genome_header

    if args['gff'] is not None:
        cds_recs = open_gff_file(args['gff'])
        genic_values, genic_header = genic_stats(seq_recs, rec_table, cds_recs, snps)
        params = params + genic_values
        header = header + genic_header

    return params, header

def main():
    args = parse_args()
    summary, header = get_stats(args)

    with open(args['out'], 'w') as fh:
        if args['with_header']:
            fh.write("\t".join([str(item) for item in header])+"\n")

        fh.write("\t".join([str(item) for item in summary])+"\n")

if __name__ == "__main__":
    main()
