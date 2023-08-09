from core_eval import *

ref = pysam.Fastafile(sys.argv[1])
sam = sys.argv[2]

print("Parsing the entire dataset to compute general stats...")
read_stats = get_sample_stats(ref, sam)
print("-" * 60)

mis, ins, delet = np.median(read_stats["metrics"], axis = 0)
print("The overall accuracy on this dataset is %.3f (mismatch %.3f insertion %.3f deletion %.3f)" % 
      (np.median(read_stats["acc"]), mis, ins, delet))
print("-" * 60)

subsample = 200000
print("Subsampling %s reads to compute systematic error profiles..." % (subsample))
kmer_profiles, motif_acc, confusion, err_dist = get_kmer_profiles(sam, ref, subsample = subsample, kmer_size=3)

print("-" * 60)

homo_total = err_dist['ho_mis'] + err_dist['ho_ins'] + err_dist['ho_del']
hetero_total = err_dist['he_mis'] + err_dist['he_ins'] + err_dist['he_del']
total = homo_total + hetero_total
print("The sequencing errors are %.1f%% heteropolymers and %.1f%% homopolymers" % 
      (hetero_total / total * 100, homo_total / total * 100))
print("-" * 60)

motif_stats = {}
for motif in motif_acc:
    total, n_del, n_ins = np.sum(list(motif_acc[motif].values())), 0, 0
    for basecall in motif_acc[motif]:
        if "_" in basecall:
            n_del += motif_acc[motif][basecall]
        if len(basecall) > 2:
            n_ins += motif_acc[motif][basecall]

    deletion = n_del / total
    insertion = n_ins / total 
    mismatch = (total - motif_acc[motif][motif] - n_del - n_ins) / total
    motif_stats[motif] = [mismatch, insertion, deletion]

print("The error profiles of 2-mer motifs")
print("\t Mis \t Ins\t Del\t")
for motif, errors in motif_stats.items():
    print("%s\t%.3f\t%.3f\t%.3f\t" % (motif, errors[0], errors[1], errors[2]))
print("-" * 60)

print("The confusion matrix of the 4 reference bases ACGU ")
confusion_matrix = np.flip(confusion / np.sum(confusion, axis = 0)[None, :], axis = 0)
basecalled = ["del", "U", "G", "C", "A"]
for i in range(len(basecalled)):
    print("%s\t%.3f\t%.3f\t%.3f\t%.3f\t" 
          % (basecalled[i], confusion_matrix[i, 0], confusion_matrix[i, 1], confusion_matrix[i, 2], confusion_matrix[i, 3]))
print("\tA\tC\tG\tU\t")