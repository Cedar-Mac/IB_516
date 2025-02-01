import merge_pairs, trim_primers, length_filter, chimera_filter, frequency_filter, quality_filter, entropy_denoise

path_to_data = "../../data/test_data"

merge_pairs.merge_pairs(data_dir=path_to_data)

trim_primers.trim_primers(data_dir=path_to_data)

length_filter.length_filter(data_dir=path_to_data, amplicon_length=142)

chimera_filter.chimera_filter(data_dir=path_to_data)

quality_filter.quality_filter(data_dir=path_to_data)

frequency_filter.frequency_filter(data_dir=path_to_data, min_seq_count=3, min_site_occurance=3)

#translation_filter()

entropy_denoise.denoise(data_dir=path_to_data)

#cluster_otus()

#get_taxon_assignments()

