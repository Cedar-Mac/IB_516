import pipeline

path_to_data = "../../data/test_data"

pipeline.merge_pairs(data_dir=path_to_data)

pipeline.trim_primers(data_dir=path_to_data)

pipeline.quality_filter(data_dir=path_to_data)

pipeline.length_filter(data_dir=path_to_data, amplicon_length=142)

pipeline.chimera_filter(data_dir=path_to_data)

pipeline.frequency_filter(data_dir=path_to_data, min_seq_count=3, min_site_occurance=3)

#translation_filter()

pipeline.denoise(data_dir=path_to_data)

#cluster_otus()

#get_taxon_assignments()

