import pysam
import re
import sys
import numpy as np 
import plotly.express as px
import statsmodels.api as sm
import plotly.graph_objects as go
import scipy
import math
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from pysam import FastxFile
from tqdm import tqdm
from utils import *
import seaborn as sns

def get_read_stats(cigar_splitted, strand):
    match, mismatch, insertion, deletion, n_clips = [], [], [], [], [0, 0]
    for index in range(0, len(cigar_splitted)-1, 2): # items[-1] is actually ""
        count, symbol = int(cigar_splitted[index]), cigar_splitted[index+1]
        if symbol == "=":
            match.append(int(count))
        elif symbol == "X":
            mismatch.append(int(count))
        elif symbol == "I":
            insertion.append(int(count))
        elif symbol == "D":
            deletion.append(int(count))
        elif symbol == "S":
            if strand == "+":
                if index == 0:
                    n_clips[0] += count
                else:
                    n_clips[1] += count
            if strand == "-":
                if index == 0:
                    n_clips[1] += count
                else:
                    n_clips[0] += count
    return match, mismatch, insertion, deletion, n_clips

def get_sample_stats(reference, sam, qscores=None, lengths=None, subsample = 10000000, return_counts = False):
    # print("Processing sample", sam)
    seqs = read_samfile(sam, subsample)
    read_stats, read_q, empirical_q = {}, [], []
    n_seqs = np.minimum(len(seqs), subsample)
    
    read_stats["metrics"] = np.empty((n_seqs, 3), dtype = np.float64)
    read_stats["acc"] = np.empty((n_seqs), dtype = np.float64)
    read_stats["qscores"] = np.empty((n_seqs, 3), dtype = np.float64)
    
    read_stats["align_ratio"] = np.empty(n_seqs, dtype = np.float64)
    read_stats["original_read_length"] = np.empty(n_seqs, dtype = np.int64)
    read_stats["aligned_read_length"] = np.empty(n_seqs, dtype = np.int64)
    read_stats["read_q"] = np.empty(n_seqs, dtype = np.float64)
    read_stats["empirical_q"] = np.empty(n_seqs, dtype = np.float64)
    read_stats["n_clips"] = np.empty((n_seqs, 2), dtype = np.float64)
    
    mismatches, insertions, deletions = [], [], []
    
    for seq_num in tqdm(range(n_seqs)):
        name, flag, start, cigar, sample_seq, mapq, chromo, qs_base = seqs[seq_num]
        # qscore, length = qscores[name], lengths[name]
        
        strand = '-' if int(flag) & 0x10 else '+'
        
        cigar_splitted = re.split('([^0-9])',cigar.strip())
        match, mismatch, insertion, deletion, n_clips = get_read_stats(cigar_splitted, strand = strand)
        
        if return_counts == True:
            mismatches += mismatch
            insertions += insertion
            deletions += deletion
        
        match, mismatch, insertion, deletion = np.sum(match), np.sum(mismatch), np.sum(insertion), np.sum(deletion)
        align_ratio = 1 - (n_clips[0] + n_clips[1]) / len(sample_seq)
        read_l = match + mismatch + insertion + deletion

        read_stats["metrics"][seq_num, 0] = mismatch/read_l
        read_stats["metrics"][seq_num, 1] = insertion/read_l
        read_stats["metrics"][seq_num, 2] = deletion/read_l
        read_stats["acc"][seq_num] = match / read_l
        
        read_stats["align_ratio"][seq_num] = align_ratio
        read_stats["n_clips"][seq_num, :] = n_clips
    
    if return_counts == True:
        return read_stats, np.array(mismatches), np.array(insertions), np.array(deletions)
    else:
        return read_stats

def read_samfile(samfile, subsample = 10000000, mapq_thres = 20, seqlen_thres = 100, return_ids = False):
    seqs, supp_count, ids = [], 0, set()
    with open(samfile, "r") as f:
        for line in f:
            if len(seqs) == subsample:
                break
            if line[0] != "@":
                records = line.strip().split("\t")
                name, flag, start, mapq, cigar, seq, chromo, qscore = records[0], records[1], records[3], records[4], records[5], records[9], records[2], records[10]
                if int(flag) & 2048:
                    supp_count += 1
                else:
                    if int(mapq) > mapq_thres and len(seq) > seqlen_thres:
                        if return_ids == True:
                            ids.add(name)
                        else:
                            seqs.append((name, flag, start, cigar, seq, mapq, chromo, qscore))
    # print("The number of supplementary alignments in this sample is", supp_count) 
    print("The number of sequences in this sample is", len(seqs))
    if return_ids == True:
        return ids
    return seqs

base2index = {"A":0, "C":1, "G":2, "U":3, "_":4}

def update_err_dist(err_dist, symbol, homo_flag):
    symbol2name = {"X": "mis", "D": "del", "I": "ins"}
    if symbol in symbol2name:
        if homo_flag:
            err_dist["ho_" + symbol2name[symbol]] += 1
        else:
            err_dist["he_" + symbol2name[symbol]] += 1
                        
def get_kmer_profiles(samfile, ref, subsample = 100000, kmer_size = 5):
    # print("Processing sample", samfile)
    seqs = read_samfile(samfile, subsample)
    
    kmer_dict, confusion = init_kmer_dict(kmer_size), np.zeros((5, 4), dtype = int)
    motif_acc = {kmer: {} for kmer in ['GA', 'AG', 'GU', 'UG', 'GC', 'CG', 'AC', 'CA', 'AU','UA', 'UC', 'CU',]}
    err_dist = {"ho_mis": 0, "ho_del": 0, "ho_ins": 0, "he_mis": 0, "he_del": 0, "he_ins": 0}
    
    for seq_num in tqdm(range(np.minimum(len(seqs), subsample))):
        name, flag, start, cigar, sample_seq, mapq, chromo, qs_base = seqs[seq_num]
        cigar_splitted, cigar_expanded, n_clips = re.split('([^0-9])', cigar.strip()), '', 0
        for index in range(0, len(cigar_splitted)-1, 2): # items[-1] is ""
            count, symbol = int(cigar_splitted[index]), cigar_splitted[index+1]
            if symbol == "S": # Is hard-clipping H dealt with?
                n_clips += count
            else:
                cigar_expanded = cigar_expanded + ''.join([symbol * count])
                
        ref_start = int(start) - 1 # samfiles are 1-based
        ref_end = ref_start + len(sample_seq) - cigar_expanded.count("I") + cigar_expanded.count("D") + cigar_expanded.count("N") - n_clips
        ref_seq = ref.fetch(reference = chromo, start = ref_start, end = ref_end)
        if "N" in ref_seq or ref_seq.upper() != ref_seq: # Filtering non-standard reference sequences?
            continue
        
        # to clip the soft-clipped bases
        sample_seq = sample_seq[int(cigar_splitted[0]):] if cigar_splitted[1] == "S" else sample_seq
        sample_seq = sample_seq[:-int(cigar_splitted[-3])] if cigar_splitted[-2] == "S" else sample_seq
            
        if int(flag) & 0x10 == '-':
            sample_seq, ref_seq, cigar_expanded = rev_compl(sample_seq), rev_compl(ref_seq), cigar_expanded[::-1]
        
        # to find the positions of homopolymers on the read, from left -1 to right + 1
        homo_pos = set()
        for kmer, homo_start, l in get_homo_pos(ref_seq, 1000000):
            for increment in range(-1, l+1, 1):
                homo_pos.add(homo_start + increment)
            
        # Expand sample sequence by adding deletions _ and removing insertions
        expanded_sample, p_sample, p_ref, insert_pos = '', 0, 0, {}
        for i in range(len(cigar_expanded)):
            
            if cigar_expanded[i] == "=" or cigar_expanded[i] == "X":
                expanded_sample += sample_seq[p_sample]
                update_err_dist(err_dist, cigar_expanded[i], p_ref in homo_pos)
                p_sample += 1
                p_ref += 1
                
            elif cigar_expanded[i] == "D":
                expanded_sample += "_"
                update_err_dist(err_dist, cigar_expanded[i], p_ref in homo_pos)
                p_ref += 1
                
            elif cigar_expanded[i] == "N":
                expanded_sample += "N"
                p_ref += 1
                
            elif cigar_expanded[i] == "I":
                insert_pos[p_ref] = sample_seq[p_sample] if sample_seq[p_sample] != "T" else "U" 
                update_err_dist(err_dist, cigar_expanded[i], p_ref in homo_pos)
                p_sample += 1
                
#         if len(expanded_sample) != len(ref_seq):
#             print("expansion didn't work!!")
        
        ref_seq = ref_seq.replace("T", "U")
        expanded_sample = expanded_sample.replace("T", "U")
                    
        for i in range(len(ref_seq)):
            if expanded_sample[i] != "N" and ref_seq[i] in base2index:
                    confusion[base2index[expanded_sample[i]], base2index[ref_seq[i]]] += 1
        
        step = 2
        for i in range(len(ref_seq)-step+1):
            kmer = ref_seq[i:i+step]
            sample_kmer = expanded_sample[i:i+step]
            
            ins_flag, inserts = False, []
            for j in range(i+1, i+step):
                if j in insert_pos:
                    ins_flag = True
                    inserts.append((j-i-1, insert_pos[j]))
            if ins_flag == True:
                for loc, insert in inserts:
                    inserted_kmer = sample_kmer[:loc+1] + insert + sample_kmer[loc+1:]
                sample_kmer = inserted_kmer

            if "N" not in sample_kmer and kmer in motif_acc:
                if sample_kmer not in motif_acc[kmer]:
                    motif_acc[kmer][sample_kmer] = 0
                motif_acc[kmer][sample_kmer] += 1
    
    return kmer_dict, motif_acc, confusion, err_dist
