import os
import pipeline
import true_errors

path_to_data = "../../data/test_data"

if "freq_filtered" not in os.listdir(path_to_data):
    pipeline.merge_pairs(data_dir=path_to_data)

    pipeline.trim_primers(data_dir=path_to_data)

    pipeline.length_filter(data_dir=path_to_data, amplicon_length=142)

    pipeline.chimera_filter(data_dir=path_to_data)

    pipeline.quality_filter(data_dir=path_to_data)

    pipeline.frequency_filter(data_dir=path_to_data, min_seq_count=3, min_site_occurance=3)

# Copied these from 
c1 = [60, 61, 62]
c2 = [69, 70, 71]
c3 = [96, 97, 98]

# Run through combinations of denoising parameters
all_stats = []
meth = "ratio"
for option in ["ent", "no_ent"]:
    for alpha in [1, 3, 5, 7, 9, 11, 13]:
        args = [str(alpha), "3", "1", "1", "4"] if option == "no_ent" else [str(alpha), "3", "1", "1", "-y", "4"]
        denoise_dir = f"{path_to_data}/benchmarking/"

        # Denoise with each parameter
        #pipeline.denoise(data_dir=path_to_data, output_dir=f"benchmarking/{option}_alpha_{alpha}_denoised", DnoisE_args=args)

        avg_err_rate, avg_seq_counts, avg_err_counts = true_errors.get_method_avgs(denoised_dir = denoise_dir, 
                                                                                   method = meth, 
                                                                                   o = option,
                                                                                   a = alpha,
                                                                                   c1 = c1, 
                                                                                   c2 = c2, 
                                                                                   c3 = c3)
        
        all_stats.append((f"{option}_alpha_{alpha}", [avg_err_rate, avg_seq_counts, avg_err_counts]))


true_errors.make_plot(method = meth,
                      stat = "sequence count",
                      stats = all_stats)
