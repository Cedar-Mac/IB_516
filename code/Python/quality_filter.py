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

    if len(vsearch_args) != 2:
        print(f"quality_filter vsearch_args must be list of length two. {len(vsearch_args)} arguments were provided. Exiting.")
        exit()

    trimmed_suffix = "_trimmed.fastq"

    # Remove old directory and files 
    if "quality_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/quality_filtered/")

    # Make directory for quality filtered reads
    os.makedirs(f"{data_dir}/quality_filtered/")

    # Make subdirectory for logs
    os.makedirs(f"{data_dir}/quality_filtered/logs/")

    for file in os.listdir(f"{data_dir}/trimmed/"):
        # vsearch call for fastq filtering. See function description for arguments
        if ".fastq" in file:
            vsearch_ee_filter_call = ["vsearch",
                                      "--fastx_filter", f"{data_dir}/trimmed/{file}",
                                      "--fastqout", f"{data_dir}/quality_filtered/{file[0:-len(trimmed_suffix)]}.fastq",
                                      "--fastq_maxee", f"{vsearch_args[0]}",
                                      "--fastq_maxns", f"{vsearch_args[1]}"]
        
            with open(f"{data_dir}/quality_filtered/logs/{file[0:-len(trimmed_suffix)]}.log", "a") as log:
                try:
                    subprocess.run(vsearch_ee_filter_call, stdout=log, stderr=log, check=True)
                    print(f"\nQuality filtered {file} successfully.\n")
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {file}: {e}\n")


if __name__ == "__main__":
    quality_filter(data_dir="../../data/test_data")