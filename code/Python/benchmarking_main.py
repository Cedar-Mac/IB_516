import os
import pipeline

path_to_data = "../../data/test_data"

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
            denoise_args = [str(alpha), "2", "1", "-y"]

        if option == "no_ent":
            denoise_args = [str(alpha), "2", "1"]

        pipeline.denoise(data_dir=path_to_data, output_dir=f"{option}_alpha_{alpha}_denoised", DnoisE_args=denoise_args)

