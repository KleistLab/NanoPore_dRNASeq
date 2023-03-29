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

def get_sample_stats(rna, reference, qscores, lengths, subsample = 10000000, return_counts = False):
    infile = samfiles[rna]
    rodan_flag = True if "rodan" in rna else False
    print("Processing sample", infile)
    seqs = read_samfile(rna, infile, subsample)
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
        qscore, length = qscores[name], lengths[name]
        
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
        read_stats["aligned_read_length"][seq_num] = match + mismatch + insertion
        read_stats["original_read_length"][seq_num] = length
        read_stats["n_clips"][seq_num, :] = n_clips
        
        # qscores are only relevant to Guppy
        if rodan_flag == False:
            read_stats["read_q"][seq_num] = qscore
            read_stats["empirical_q"][seq_num] = symbol2qscore(qs_base)
            clip_start, clip_end = get_clipped_bases(cigar_splitted)
            qscore_clipped = symbol2qscore(qs_base[clip_start : -clip_end-1])
            read_stats["qscores"][seq_num, 0] = qscore
            read_stats["qscores"][seq_num, 1] = qscore_clipped
        
        # fit the linear model between mean_qscore_template and empirical mean qscore
    if rodan_flag == False:
        mod = sm.OLS(read_stats["empirical_q"], sm.add_constant(read_stats["read_q"]))
        const, x1 = mod.fit().params    
        read_stats["qscores"][:, 2] = np.divide(np.add(read_stats["qscores"][:, 1], -const), x1)
        read_stats["fitted_params"] = (const, x1)
    
    if return_counts == True:
        return read_stats, np.array(mismatches), np.array(insertions), np.array(deletions)
    else:
        return read_stats

def read_samfile(rna, samfile, subsample = 10000000, mapq_thres = 20, seqlen_thres = 100, return_ids = False):
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
    print("The number of supplementary alignments in this sample is ", supp_count) 
    print("The number of sequences in this sample is ", len(seqs))
    if return_ids == True:
        return ids
    return seqs


