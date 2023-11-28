import pysam
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import re
import statsmodels.api as sm

analysis_set = {}

analysis_set["main"] = ["Human", "Human_IVT",  "Mouse",  "Zebrafish", "C. elegans", "Arabidopsis",  "H. volcanii", "E. coli", "shortRNAs_IVT", "Yeast", "SARS2", "SARS2_IVT"]
 
    
refs, fastq_files, summary_files = {}, {}, {}

#refs["SARS2_IVT"] = pysam.Fastafile("sars2/sars2ivt/SARS-CoV-2.fa") # CTTT instead of CTTC in the paper...
#refs["SARS2"] = refs["SARS2_IVT"]

#refs["shortRNAs_IVT"] = pysam.Fastafile("nanorms/small_RNAs.fasta") 

#refs["Yeast"] = pysam.Fastafile("epinano/epinano_yeast_rep2/Yeast_sk1.fa")

#refs["Human"] = pysam.Fastafile("eligos/HS_H460/GCF_000001405.40_GRCh38.p14_genomic.fna")
#refs["Human_IVT"] = refs["Human"]

#refs["Mouse"] = pysam.Fastafile("mouse/GCF_000001635.27_GRCm39_genomic.fna")
#refs["Zebrafish"] = pysam.Fastafile("nano3pseq/GCF_000002035.6_GRCz11_genomic.fna")
#refs["H. volcanii"] = pysam.Fastafile("regensburg/hvo/GCF_000025685.1_ASM2568v1_genomic.fna")
#refs["E. coli"] = pysam.Fastafile("regensburg/ecoli/escherichia_000005845.2_ASM584v2_genomic.fna")
#refs["Arabidopsis"] = pysam.Fastafile("parkersimpson2020/arabidopsis_TAIR10.1_genomic.fna")
#refs["C. elegans"] = pysam.Fastafile("celegans/GCF_000002985.6_WBcel235_genomic.fna")
    
samfiles = {"Human_IVT": "eligos/human_ivt/human_ivt_guppy.sam",
             "Human": "eligos/HS_H460/h460_guppy.sam",
            
            "Mouse": "mouse/mouse_guppy.sam",

            "SARS2": "sars2/sars2_guppy.sam",
            "SARS2_IVT": "sars2/sars2ivt/sars2ivt_guppy_full.sam",

             "Yeast": "epinano/epinano_yeast_rep2/yeast_guppy.sam", 

             "Zebrafish":"nano3pseq/zebrafish_guppy.sam",
            "H. volcanii": "regensburg/hvo/hvo_guppy.sam",
            "E. coli": "regensburg/ecoli/ecoli_guppy.sam",
            "Arabidopsis": "parkersimpson2020/arabidopsis_guppy.sam",
            "C. elegans": "celegans/celegans_guppy.sam",
            "shortRNAs_IVT": "nanorms/small_rna_guppy.sam"
            }

rodan_trainingset = set(["Human", "E. coli", "Arabidopsis", "C. elegans"])
analysis_set["rodan"] = []
for sample in analysis_set["main"]:
    key =  sample + "_rodan"
    analysis_set["rodan"].append(key)
 #   refs[key] = refs[sample]
    samfiles[key] = samfiles[sample].replace("guppy", "rodan")
    
for sample in samfiles.keys():
    paths = samfiles[sample].split("/")
    fastq_dirs = ''
    for path in paths[:-1]:
        fastq_dirs += path
        fastq_dirs += '/'
    fastq_files[sample] = fastq_dirs + "guppy_6.1.7/merged.fastq"
    summary_files[sample] = fastq_dirs + "guppy_6.1.7/sequencing_summary.txt"

    
def symbol2qscore(symbols):
    quals = np.fromiter(
            (ord(x) - 33 for x in symbols),
            dtype=int, count=len(symbols))
    mean_p = np.mean(np.power(10, quals/-10))
    return -10*np.log10(mean_p)

def get_clipped_bases(cigar_splitted):
    clip_start = int(cigar_splitted[0]) if cigar_splitted[1] == "S" else 0
    clip_end = int(cigar_splitted[-3]) if cigar_splitted[-2] == "S" else 0
    return clip_start, clip_end

def init_homodict(min_size, max_size):
    homos = set()
    alphabet = ["A", "C", "G", "U"]
    for i in range(min_size, max_size + 1):
        for letter in alphabet:
            homos.add(letter * i)
    return homos

def get_homo_pos(seq, max_len = 9):
    homo2pos = []
    for homopolymer in re.finditer(r'([ACGU])\1{1,}', seq):
        homo_seq = homopolymer.group()
        homo_start = homopolymer.start()
        l = len(homo_seq)
        if l <= max_len:
            homo2pos.append((homo_seq, homo_start, l))
    return homo2pos

def is_homop(seq, letter):
    if "_" in seq:
        seq = "".join(seq.split("_"))
    if (len(set(seq)) == 1 and letter == seq[0]) or seq == '':
        return seq
    return False


def rev_compl(seq):
    base2comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(base2comp[base] for base in reversed(seq))


def get_summaries(summaryfile, subsample = None):
    qscores, lengths = {},{}
    with open(summaryfile, "r") as f:
        next(f)
        for line in f:
            if subsample is not None and len(qscores) < subsample:
                items = line.strip().split()
                qscores[items[1]] = float(items[14])
                lengths[items[1]] = float(items[13])
            if subsample is None:
                items = line.strip().split()
                qscores[items[1]] = float(items[14])
                lengths[items[1]] = float(items[13])
    return qscores, lengths

def init_kmer_dict(k = 5):
    kmer_dict, kmers = {}, ["A", "C", "G", "U"]
    alphabet = set("ACGU")
    
    for i in range(k-1):
        kmers_extended = []
        for kmer in kmers:
            for letter in alphabet:
                kmers_extended.append(kmer+letter)
        kmers = kmers_extended
        
    for kmer in kmers:
        kmer_dict[kmer] = dict()
        
    return kmer_dict


def get_kmer_profile(expanded_sample, ref_seq, insert_pos, kmer_size = 5):
    kmer_dict = init_kmer_dict(kmer_size)
    if len(expanded_sample) != len(ref_seq):
        print("expansion didn't work!!")
    for i in range(len(ref_seq)-kmer_size+1):
        kmer = ref_seq[i:i+kmer_size]
        sample_kmer = expanded_sample[i:i+kmer_size]
        ins_flag = False
        base2insert = []
        for j in range(i+1, i+kmer_size):
            if j in insert_pos:
                ins_flag = True
                base2insert.append((j-i, insert_pos[j]))
        if ins_flag == True:
            for loc, insert in base2insert:
                inserted_kmer = sample_kmer[:loc] + insert + sample_kmer[loc:]
            sample_kmer = inserted_kmer
            # print(i, kmer, sample_kmer, inserted_kmer, base2insert)
        if kmer in kmer_dict:
            if sample_kmer not in kmer_dict[kmer]:
                kmer_dict[kmer][sample_kmer] = 0
            kmer_dict[kmer][sample_kmer] += 1
    return kmer_dict

    
def plot_qual_length(qscores, lengths, sample, aligned_reads):
    x = np.arange(2, 15)
    counts, counts_aligned, lengths_sum = np.zeros(len(x)), np.zeros(len(x)), []
    for coord in x:
        counts[int(coord-2)] = 0
        counts_aligned[int(coord-2)] = 0
        lengths_sum.append([])
    for read_id in qscores:
        qs, l = qscores[read_id], lengths[read_id]
        index = int(np.floor(qs))-2
        if index < len(counts):
            counts[index] += 1
            lengths_sum[index].append(l)
            if read_id in aligned_reads:
                counts_aligned[index] += 1

    lengths_mean, lengths_median = np.zeros(len(x)), np.zeros(len(x))
    for i in range(len(lengths_sum)):
        lengths_mean[i] = np.mean(lengths_sum[i])
        lengths_median[i] = np.median(lengths_sum[i])
        
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.bar(x, counts, color = "#69b3a2", edgecolor = "black")
    ax1.bar(x, counts_aligned, color = "black", edgecolor = "black")
    # print(lengths_sum, counts)
    ax2.plot(x, lengths_mean, color = "#3399e6", lw=3, label = "mean read length")
    ax2.plot(x, lengths_median, color = 'red', lw=3, label = "median read length")
    ax2.legend(loc="upper left")
    ax1.set_xlabel("qscore"), ax1.set_ylabel("count"), ax2.set_ylabel("length")
    fig.suptitle(sample + " " + str(len(qscores)) + " reads")
    plt.show()
    
    
def get_cumSum(concent):
    cummulative = []
    partial = 0
    for concent in concent:
        partial += concent 
        cummulative.append(partial)
    return cummulative

def get_cumSum(concent):
    cummulative = []
    partial = 0
    for concent in concent:
        partial += concent 
        cummulative.append(partial)
    return cummulative

def get_content(kmer2acc):
    alphabet = ["G", "A", "C", "T"]
    one_content, one_cummu, two_content, two_cummu = {}, {}, {}, {}
    for letter in alphabet:
        one_content[letter] = []
    for kmer, acc in kmer2acc: 
        for letter in alphabet:
            one_content[letter].append(kmer.count(letter))

    for letter in alphabet:
        one_cummu[letter] = get_cumSum(one_content[letter])
    
    for i in range(4):
        for j in range(i, 4):
            two_content[alphabet[i] + alphabet[j]] = []

    for kmer, acc in kmer2acc:
        for twoMer in two_content:
            twoMer_flipped = twoMer[1] + twoMer[0]
            two_content[twoMer].append(kmer.count(twoMer) + kmer.count(twoMer_flipped))

    for twoMer in two_content:
        two_cummu[twoMer] = get_cumSum(two_content[twoMer])
    
    return one_cummu, two_cummu


def get_qscore(qualityfile):
    qscores = {}
    with open(qualityfile, "r") as f:
        next(f)
        for line in f:
            items = line.strip().split()
            qscores[items[1]] = float(items[14])
    print(qualityfile, len(qscores))
    return qscores

