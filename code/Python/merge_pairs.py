import os, re, subprocess, concurrent

def merge_pairs(data_dir:str):
    data_files = os.listdir(data_dir)

    fastq_suffix = "_RX.fastq"

    # get file names of forward and reverse reads, strip suffix
    fwd_files = [fname for fname in data_files if re.search(pattern = "_R1.fastq", string = fname)]
    rvs_files = [fname for fname in data_files if re.search(pattern = "_R2.fastq", string = fname)]

    # Check that all reads are paired
    only_in_fwd_reads = [s for s in fwd_files if any(substr in s for substr in rvs_files)]
    only_in_rvs_reads = [s for s in fwd_files if any(substr in s for substr in rvs_files)]
    unmatched_files = only_in_fwd_reads + only_in_rvs_reads

    if unmatched_files:
        print("Not all reads are paired!")
        exit()

    # make output directory for merged files
    if "merged" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/merged/")
    #if "logs" not in os.listdir(data_dir):
        #os.makedirs(f"{data_dir}/logs/") 

    # Use VSEARCH fastq_mergepairs function
    for fastq_name in fwd_files:
        vsearch_merge_call = ["vsearch", "--fastq_mergepairs", f"{data_dir}/{fastq_name}", 
                                "--reverse", f"{data_dir}/{fastq_name}", 
                                #"--eetabbedout", f"{data_dir}/logs/{fastq_name}.log", 
                                "--fastqout", f"{data_dir}/merged/{fastq_name[0:-len(fastq_suffix)]}_merged.fastq"]
        subprocess.check_call(vsearch_merge_call)


merge_pairs("../../data/test_data")