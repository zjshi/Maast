"""dnds
This module is a reference implementation of estimating nucleotide substitution
neutrality by estimating the percent of synonymous and nonsynonymous mutations.
"""
from __future__ import print_function, division
from math import log
from fractions import Fraction
import logging

BASES = {'A', 'G', 'T', 'C'}

codons = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "STOP",
    "TAG": "STOP",
    "TGT": "C",
    "TGC": "C",
    "TGA": "STOP",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}

ff_deg_sites = {
    'GGG.2': 1,
    'ACT.2': 1,
    'GCA.2': 1,
    'ACG.2': 1,
    'CTC.2': 1,
    'TCC.2': 1,
    'ACC.2': 1,
    'GTA.2': 1,
    'GTG.2': 1,
    'GCG.2': 1,
    'GGA.2': 1,
    'CCG.2': 1,
    'GGC.2': 1,
    'CCT.2': 1,
    'ACA.2': 1,
    'CGC.2': 1,
    'GCT.2': 1,
    'GTT.2': 1,
    'TCG.2': 1,
    'CGA.2': 1,
    'CTT.2': 1,
    'CGT.2': 1,
    'TCT.2': 1,
    'CGG.2': 1,
    'CCC.2': 1,
    'CTG.2': 1,
    'CCA.2': 1,
    'GGT.2': 1,
    'GCC.2': 1,
    'CTA.2': 1,
    'TCA.2': 1,
    'GTC.2': 1
}

non_deg_sites = {
    'GCA.0': 1,
    'AAA.0': 1,
    'AAA.1': 1,
    'GCA.1': 1,
    'CGT.1': 1,
    'CGT.0': 1,
    'TTG.1': 1,
    'AGG.1': 1,
    'CCT.1': 1,
    'CCT.0': 1,
    'AGA.1': 1,
    'AAG.0': 1,
    'AAG.1': 1,
    'GAT.1': 1,
    'GAT.0': 1,
    'TTA.1': 1,
    'AGC.0': 1,
    'AGC.1': 1,
    'GAC.0': 1,
    'GTA.0': 1,
    'TGA.2': 1,
    'TTC.0': 1,
    'TTC.1': 1,
    'GTG.0': 1,
    'TGT.0': 1,
    'TGT.1': 1,
    'ACT.1': 1,
    'ACT.0': 1,
    'AGT.1': 1,
    'AGT.0': 1,
    'CCA.0': 1,
    'CCA.1': 1,
    'TGG.1': 1,
    'TGG.0': 1,
    'TGG.2': 1,
    'GGA.0': 1,
    'GGA.1': 1,
    'CAA.0': 1,
    'CAA.1': 1,
    'GTT.0': 1,
    'GTT.1': 1,
    'CGC.0': 1,
    'CGC.1': 1,
    'CAC.0': 1,
    'CAC.1': 1,
    'ATG.1': 1,
    'ATG.0': 1,
    'ATG.2': 1,
    'ATA.1': 1,
    'ATA.0': 1,
    'CGG.1': 1,
    'TAT.0': 1,
    'TAT.1': 1,
    'CAG.0': 1,
    'CAG.1': 1,
    'CTG.1': 1,
    'GAA.0': 1,
    'GAA.1': 1,
    'ATC.1': 1,
    'ATC.0': 1,
    'GGG.0': 1,
    'GGG.1': 1,
    'GAC.1': 1,
    'GTA.1': 1,
    'CTC.1': 1,
    'TGA.0': 1,
    'TAA.0': 1,
    'GTG.1': 1,
    'CTC.0': 1,
    'GCG.0': 1,
    'GCG.1': 1,
    'CCG.0': 1,
    'CCG.1': 1,
    'TGC.1': 1,
    'TGC.0': 1,
    'GGC.0': 1,
    'GGC.1': 1,
    'CAT.1': 1,
    'CAT.0': 1,
    'CGA.1': 1,
    'CTT.0': 1,
    'CTT.1': 1,
    'TCT.0': 1,
    'TCT.1': 1,
    'CCC.0': 1,
    'CCC.1': 1,
    'ATT.0': 1,
    'ATT.1': 1,
    'GAG.0': 1,
    'GAG.1': 1,
    'GGT.1': 1,
    'GGT.0': 1,
    'GCC.0': 1,
    'GCC.1': 1,
    'AAT.1': 1,
    'AAT.0': 1,
    'GTC.1': 1,
    'GTC.0': 1,
    'ACC.0': 1,
    'ACG.0': 1,
    'ACG.1': 1,
    'ACC.1': 1,
    'TCC.1': 1,
    'TCC.0': 1,
    'TAG.1': 1,
    'TAG.0': 1,
    'TAC.1': 1,
    'TAC.0': 1,
    'GCT.1': 1,
    'GCT.0': 1,
    'ACA.0': 1,
    'ACA.1': 1,
    'TCG.1': 1,
    'TCG.0': 1,
    'TTT.1': 1,
    'TTT.0': 1,
    'AAC.0': 1,
    'AAC.1': 1,
    'CTA.1': 1,
    'TCA.1': 1,
    'TCA.0': 1
}

def split_seq(seq, n=3):
    '''Returns sequence split into chunks of n characters, default is codons'''
    return [seq[i:i + n] for i in range(0, len(seq), n)]

def average_list(l1, l2):
    """Return the average of two lists"""
    return [(i1 + i2) / 2 for i1, i2 in zip(l1, l2)]

def dna_to_protein(codon):
    '''Returns single letter amino acid code for given codon'''
    if codon in codons:
        return codons[codon]
    else:
        return 'Z'

def translate(seq):
    """Translate a DNA sequence into the 1-letter amino acid sequence"""
    return "".join([dna_to_protein(codon) for codon in split_seq(seq)])

def is_ambiguous(codon):
    if codon not in codons:
        return True
    else:
        return False

def is_synonymous(codon1, codon2):
    '''Returns boolean whether given codons are synonymous'''
    return dna_to_protein(codon1) == dna_to_protein(codon2)

def is_stop_codon(codon):
    if dna_to_protein(codon) == 'STOP':
        return True
    else:
        return False

def is_4f_deg(codon, pos):
    if "{}.{}".format(codon, pos) in ff_deg_sites:
        return True
    else:
        return False

def is_non_deg(codon, pos):
    if "{}.{}".format(codon, pos) in non_deg_sites:
        return True
    else:
        return False

def dnds_codon(codon):
    '''Returns list of synonymous counts for a single codon.
    Calculations done per the methodology taught in class.
    http://sites.biology.duke.edu/rausher/DNDS.pdf
    '''
    syn_list = []
    for i in range(len(codon)):
        base = codon[i]
        other_bases = BASES - {base}
        syn = 0
        for new_base in other_bases:
            new_codon = codon[:i] + new_base + codon[i + 1:]
            syn += int(is_synonymous(codon, new_codon))
        syn_list.append(Fraction(syn, 3))
    return syn_list


def dnds_codon_pair(codon1, codon2):
    """Get the dN/dS for the given codon pair"""
    return average_list(dnds_codon(codon1), dnds_codon(codon2))


def syn_sum(seq1, seq2):
    """Get the sum of synonymous sites from two DNA sequences"""
    syn = 0
    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]
        dnds_list = dnds_codon_pair(codon1, codon2)
        syn += sum(dnds_list)
    return syn


def hamming(s1, s2):
    """Return the hamming distance between 2 DNA sequences"""
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)) + abs(len(s1) - len(s2))

def codon_subs(codon1, codon2):
    """Returns number of synonymous substitutions in provided codon pair
    Methodology for multiple substitutions from Dr. Swanson, UWashington
    https://faculty.washington.edu/wjs18/dnds.ppt
    """
    diff = hamming(codon1, codon2)
    if diff == 0:
        return 0
    elif diff == 1:
        return int(translate(codon1) == translate(codon2))
    elif diff == 2:
        syn = 0
        for i in range(len(codon1)):
            base1 = codon1[i]
            base2 = codon2[i]
            if base1 != base2:
                new_codon = codon1[:i] + base2 + codon1[i + 1:]
                syn += int(is_synonymous(codon1, new_codon))
                syn += int(is_synonymous(codon2, new_codon))
        return syn / diff
    elif diff == 3:
        syn = 0
        neighbors1 = []
        neighbors2 = []

        for i in range(len(codon1)):
            base1 = codon1[i]
            base2 = codon2[i]

            neighbor1 = codon1[:i] + base2 + codon1[i + 1:]
            neighbor2 = codon2[:i] + base1 + codon2[i + 1:]

            neighbors1.append(neighbor1)
            neighbors2.append(neighbor2)

            syn += 2*int(is_synonymous(codon1, neighbor1))
            syn += 2*int(is_synonymous(codon2, neighbor2))

        for neighbor1 in neighbors1:
            for neighbor2 in neighbors2:
                if hamming(neighbor1, neighbor2) == 1:
                    syn += int(is_synonymous(neighbor1, neighbor2))

        return syn/(diff*2)
    else:
        assert False


def substitutions(seq1, seq2):
    """Returns number of synonymous and nonsynonymous substitutions"""
    dna_changes = hamming(seq1, seq2)

    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)

    syn = 0
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]

        syn += codon_subs(codon1, codon2)

    return (syn, dna_changes - syn)

def simple_count(codon1, codon2):
    diff = hamming(codon1, codon2)

    syn = 0
    prm_stop = 0
    dsrpt_stop = 0
    for i in range(len(codon1)):
        base1 = codon1[i]
        base2 = codon2[i]
        if base1 != base2:
            new_codon = codon1[:i] + base2 + codon1[i + 1:]
            if is_synonymous(codon1, new_codon):
                syn = syn + 1
            else:
                if is_stop_codon(codon1):
                    dsrpt_stop = dsrpt_stop + 1
                elif is_stop_codon(new_codon):
                    prm_stop = prm_stop + 1

    return syn, diff-syn, prm_stop, dsrpt_stop

def snp_stats(seq1, seq2):
    n_snp = hamming(seq1, seq2)

    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)

    syns = 0
    non_syns = 0
    prm_stops = 0
    dsrpt_stops = 0

    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]

        if hamming(codon1, codon2) == 0:
            pass
        else:
            syn, non_syn, prm_stop, dsrpt_stop = simple_count(codon1, codon2)
            syns = syns + syn
            non_syns = non_syns + non_syn
            prm_stops = prm_stop + prm_stop
            dsrpt_stops = dsrpt_stops + dsrpt_stop

    return syns, non_syns, prm_stops, dsrpt_stops

def clean_sequences(seq1, seq2):
    """Clean up provided sequence by removing whitespace."""
    seq1 = seq1.replace(' ', '')
    seq2 = seq2.replace(' ', '')

    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)

    new_codon_list1 = []
    new_codon_list2 = []

    n_wild_codons = 0
    for i, codon in enumerate(codon_list1):
        if is_ambiguous(codon_list1[i]) or is_ambiguous(codon_list2[i]):
            n_wild_codons = n_wild_codons + 1
            pass
        else:
            new_codon_list1.append(codon_list1[i])
            new_codon_list2.append(codon_list2[i])

    seq1 = "".join(new_codon_list1)
    seq2 = "".join(new_codon_list2)

    print("codons containing wildcards (stripped): {}".format(n_wild_codons))

    return seq1, seq2

def pnps_dnds(seq1, seq2):
    """Main function to calculate pN/pS between two DNA sequences."""
    """Main function to calculate dN/dS between two DNA sequences per Nei &
    Gojobori 1986. This includes the per site conversion adapted from Jukes &
    Cantor 1967.
    """
    # Strip any whitespace from both strings and remove wildcard codons
    seq1, seq2 = clean_sequences(seq1, seq2)
    # Check that both sequences have the same length
    assert len(seq1) == len(seq2)
    # Check that sequences are codons
    assert len(seq1) % 3 == 0

    snp_stats(seq1, seq2)


    syn_sites = syn_sum(seq1, seq2)
    non_sites = len(seq1) - syn_sites
    logging.info('Sites (syn/nonsyn): {}, {}'.format(syn_sites, non_sites))
    syn_subs, non_subs = substitutions(seq1, seq2)
    logging.info('pN: {} / {}\t\tpS: {} / {}'.format(non_subs, round(non_sites), syn_subs, round(syn_sites)))

    pn = ps = dn = ds = pn2ps = dn2ds = -1

    try:
        pn = non_subs / non_sites
    except Exception as error:
        pass

    try:
        ps = syn_subs / syn_sites
    except Exception as error:
        pass

    try:
        dn = -(3 / 4) * log(1 - (4 * pn / 3))
    except Exception as error:
        pass

    try:
        ds = -(3 / 4) * log(1 - (4 * ps / 3))
    except Exception as error:
        pass

    logging.info('dN: {}\t\tdS: {}'.format(round(dn, 3), round(ds, 3)))

    if pn == -1 or ps == -1:
        pass
    else:
        try:
            pn2ps = pn / ps
        except Exception as error:
            pass

    if dn == -1 or ds == -1:
        pass
    else:
        try:
            dn2ds = dn / ds
        except Exception as error:
            pass

    return pn, ps, dn, ds, pn2ps, dn2ds
