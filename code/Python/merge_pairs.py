import os, re, subprocess

def merge_pairs(data_dir:str):
    """
    Merges paired fastq reads using vsearch algorithm.

    Arguments:

        data_dir: string providing the path to the data directory with raw input fastq files.

    vsearch arguments:

    --fastq_mergepairs: specify the maximum number of non-matching nucleotides 
                        allowed in the overlap region. That option has a strong 
                        influence on the merging success rate. Set to 99 a la JAMP.
    --fastq_minovlen: Minimum overlap between paired sequences to merge. JAMP default is 16,
                        base default is 10. Currently set to JAMP default.
    --fastq_diffpct: Maximum allowed percent difference between the two paired sequences. 
                        vsearch default is 100%, JAMP sets the default to 25%
                        I'm using vsearch default.
    --fastq_allowmergestagger: Allows for reads to be staggered, this is following JAMP convention.

    Output:

    Merged fastq files in merged directory within data directory.
    """

    data_files = os.listdir(data_dir)

    fastq_suffix = "_RX.fastq"

    # get file names of forward and reverse reads, strip suffix
    fwd_files = [fname for fname in data_files if re.search(pattern = "_R1.fastq", string = fname)]
    rvs_files = [fname for fname in data_files if re.search(pattern = "_R2.fastq", string = fname)]

    # Check that all reads are paired 
    fwd_bare_names = set([name[0:-len(fastq_suffix)] for name in fwd_files])
    rvs_bare_names = set([name[0:-len(fastq_suffix)] for name in rvs_files])

    if fwd_bare_names != rvs_bare_names:
        print("Not all reads are paired!")
        exit()

    # make output directory for merged files
    if "merged" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/merged/")
    #if "logs" not in os.listdir(data_dir):
        #os.makedirs(f"{data_dir}/logs/") 

    # Use VSEARCH fastq_mergepairs function
    for fastq_name in fwd_files:
        vsearch_merge_call = ["vsearch", 
                                "--fastq_mergepairs", f"{data_dir}/{fastq_name}", 
                                "--reverse", f"{data_dir}/{fastq_name}", 
                                "--fastqout", f"{data_dir}/merged/{fastq_name[0:-len(fastq_suffix)]}_merged.fastq",
                                #"-log", f"{data_dir}/logs/{fastq_name}[0:-len(fastq_suffix)].log", 
                                "--fastq_maxdiffs", "99", 
                                "--fastq_minovlen", "16",
                                "--fastq_maxdiffpct", "25",
                                "--fastq_allowmergestagger"]


        subprocess.check_call(vsearch_merge_call)


merge_pairs("../../data/test_data")
