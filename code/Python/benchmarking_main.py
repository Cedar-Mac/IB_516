import os
import pipeline
import true_errors

path_to_data = "../../data/test_data"

fasta_suffix = ".fasta"

errors = {}

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

        if option == "ent":
            args = [str(alpha), "3", "1", "1", "-y"]

        if option == "no_ent":
            args = [str(alpha), "3", "1", "1"]

        #pipeline.denoise(data_dir=path_to_data, output_dir=f"benchmarking/{option}_alpha_{alpha}_denoised", DnoisE_args=args)

        # Calculate errors for each fasta file
        for file in os.listdir(f"{path_to_data}/benchmarking/{option}_alpha_{alpha}_denoised/"):
            if (".csv" not in file) and ("log" not in file):
                file_error_rate = true_errors.check_seqs(data_dir = f"{path_to_data}/benchmarking/{option}_alpha_{alpha}_denoised",
                                        fasta = file,
                                        codon_1= c1,
                                        codon_2 = c2,
                                        codon_3 = c3)
                errors[f"{option}_alpha_{alpha}_{file[0:-len(fasta_suffix)]}"] = file_error_rate

print(errors)