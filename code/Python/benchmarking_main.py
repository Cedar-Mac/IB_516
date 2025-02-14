import os
import pipeline
import true_errors

path_to_data = "../../data/test_data/freq_filtered"
output_path = "../../data/test_data/benchmarking/"

if not os.listdir(path_to_data):
    pipeline.merge_pairs(data_dir=path_to_data)

    pipeline.trim_primers(data_dir=path_to_data)

    pipeline.length_filter(data_dir=path_to_data, amplicon_length=142)

    pipeline.chimera_filter(data_dir=path_to_data)

    pipeline.quality_filter(data_dir=path_to_data)

    pipeline.frequency_filter(data_dir=path_to_data, min_seq_count=3, min_site_occurance=3)

# Copied these from true_errors.get_conserved_base_positions()
c1 = [60, 61, 62]
c2 = [69, 70, 71]
c3 = [96, 97, 98]

# Run through combinations of denoising parameters
all_stats = []
for option in ["dnoise", "unoise"]:
    for alpha in [1, 3, 5, 7, 9, 11, 13]:
        dnoise_args = [str(alpha), "3", "1", "-y"] 
        unoise_args = ["1", str(alpha)]

        # Denoise with each parameter
#        if option == "dnoise":
#            pipeline.denoise(data_dir=path_to_data, 
#                             output_dir=f"{output_path}/{option}_alpha_{alpha}", 
#                             option=option,
#                             DnoisE_args=dnoise_args)

#        if option == "unoise":
#            pipeline.denoise(data_dir=path_to_data, 
#                             output_dir=f"{output_path}/{option}_alpha_{alpha}", 
#                             option=option,
#                             Unoise_args=unoise_args)

        avg_err_rate, avg_seq_counts, avg_err_counts = true_errors.get_method_avgs(denoised_dir = output_path,
                                                                                   o = option,
                                                                                   a = alpha,
                                                                                   c1 = c1, 
                                                                                   c2 = c2, 
                                                                                   c3 = c3)
        
        all_stats.append((f"{option}_alpha_{alpha}", [avg_err_rate, avg_seq_counts, avg_err_counts]))


true_errors.make_plot(stat = "error rate",
                      stats = all_stats)
