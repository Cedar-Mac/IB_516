import os

trimmed_suffix = "_trimmed.fastq"

def length_filter(data_dir:str):
    # make output directory for length-filtered files
    if "len_filtered" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/len_filtered/")
    

    for in_file in os.listdir(f"{data_dir}/quality_filter/"):
        out_file = open(f"{data_dir}/len_filtered/{in_file[0:-len(trimmed_suffix)]}.fastq", "a")
        for line in in_file:
            if len(line) == 142:
                out_file.write(line)


length_filter("../../data/test_data")

