import os, subprocess, shutil

def quality_filter(data_dir:str, vsearch_args:list=[1, 0]):
    """
    Filter fastq filters that have already been trimmed.

    Input:
        - data_dir: String with main data directory

    vsearch_args:
        - [0] --fastq_maxee: Max expected cummulative error.
                            Set to 1 (probability of single error = )
        - [1] --fastq_maxns: Number of allowed N's in the sequence.
                            Set to 0

    Output:
        - quality filtered reads in subfolder of data directory
    """

    trimmed_suffix = "_trimmed.fastq"

    # Remove old directory and files 
    if "quality_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/quality_filtered/")

    # Make directory for quality filtered reads
    os.makedirs(f"{data_dir}/quality_filtered/")

    for data_file in os.listdir(f"{data_dir}/trimmed/"):
        # vsearch call for fastq filtering. See function description for arguments
        if ".fastq" in data_file:
            vsearch_ee_filter_call = ["vsearch",
                                      "--fastx_filter", f"{data_dir}/trimmed/{data_file}",
                                      "--fastqout", f"{data_dir}/quality_filtered/{data_file[0:-len(trimmed_suffix)]}.fastq",
                                      "--fastq_maxee", f"{vsearch_args[0]}",
                                      "--fastq_maxns", f"{vsearch_args[1]}"]
        
        subprocess.check_call(vsearch_ee_filter_call)


if __name__ == "__main__":
    quality_filter(data_dir="../../data/test_data")