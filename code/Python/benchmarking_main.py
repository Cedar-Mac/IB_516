import os, re
import pipeline
import true_errors
import matplotlib.pyplot as plt

path_to_data = "../../data/test_data"

fasta_suffix = ".fasta"

errors = []
seq_count = []
error_count = []

c1 = [60, 61, 62]
c2 = [69, 70, 71]
c3 = [96, 97, 98]

if "freq_filtered" not in os.listdir(path_to_data):
    pipeline.merge_pairs(data_dir=path_to_data)

    pipeline.trim_primers(data_dir=path_to_data)

    pipeline.length_filter(data_dir=path_to_data, amplicon_length=142)

    pipeline.chimera_filter(data_dir=path_to_data)

    pipeline.quality_filter(data_dir=path_to_data)

    pipeline.frequency_filter(data_dir=path_to_data, min_seq_count=3, min_site_occurance=3)

# Run through combinations of denoising parameters
for option in ["ent", "no_ent"]:
    for alpha in [1, 3, 5, 7, 9, 11, 13]:
        args = [str(alpha), "3", "1", "1", "4"] if option == "no_ent" else [str(alpha), "3", "1", "1", "-y", "4"]
        denoise_dir = f"{path_to_data}/benchmarking/{option}_alpha_{alpha}_denoised/"

        # Denoise with each parameter
        pipeline.denoise(data_dir=path_to_data, output_dir=f"benchmarking/{option}_alpha_{alpha}_denoised", DnoisE_args=args)

        # Calculate errors for each fasta file
        method_stats = []
        for file in os.listdir(denoise_dir):
            if ".csv" not in file and "log" not in file:
                file_error_rate = true_errors.check_seqs(
                    data_dir=denoise_dir, fasta=file, codon_1=c1, codon_2=c2, codon_3=c3
                )
                method_stats.append(file_error_rate)

        # Ensure method_stats is not empty before calculations
        if not method_stats:
            print(f"No valid FASTA files found for {option}_alpha_{alpha}. Skipping calculations.")
            continue

        error_rate = [fasta_stats[2] for fasta_stats in method_stats]
        num_seqs = [fasta_stats[1] for fasta_stats in method_stats]
        total_errors = [fasta_stats[0] for fasta_stats in method_stats]

        # Handle division by zero safely
        method_error_rate = (f"{option}_alpha_{alpha}", sum(error_rate) / len(error_rate)) if error_rate else (f"{option}_alpha_{alpha}", 0)
        method_num_seqs = (f"{option}_alpha_{alpha}", sum(num_seqs) / len(num_seqs)) if num_seqs else (f"{option}_alpha_{alpha}", 0)
        method_errors = (f"{option}_alpha_{alpha}", sum(total_errors) / len(total_errors)) if total_errors else (f"{option}_alpha_{alpha}", 0)

        errors.append(method_error_rate)
        seq_count.append(method_num_seqs)
        error_count.append(method_errors)
print(seq_count)


plt.style.use('_mpl-gallery')
fig, ax = plt.subplots()

x = [1, 3, 5, 7, 9, 11, 13]
x.reverse()

y1 = [item[1] for item in error_count if "no" not in item[0]]
y1.reverse()

y2 = [item[1] for item in error_count if "no" in item[0]]
y2.reverse()

ax.plot(x, y1, label="Entropy Correction")
ax.plot(x, y2, linestyle="dashed", label="No Entropy")
ax.set_xlabel("alpha value")
ax.set_ylabel("Error Rate")
ax.set_ylim(0, 0.05)
ax.set_xticklabels(x)
ax.legend()

plt.show()
