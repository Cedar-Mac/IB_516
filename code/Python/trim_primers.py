import os, subprocess

def trim_primers(data_dir: str, primer_option:int=1):
    merged_suffix = "_merged.fastq"
    fwd_primer = ""
    rvs_primer = ""
    
    # default primer options uses the EPT COI primer set
    if primer_option == 1:
        fwd_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
        rvs_primer = "CAAACAAATARDGGTATTCGDTY"
    
    # make output directory for trimmed files
    if "trimmed" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/trimmed/")

    for file in os.listdir(f"{data_dir}/merged/"):
        cutadapt_call = ["cutadapt", 
                         "-a", f"^{fwd_primer}...{rvs_primer}$", 
                         "-o", f"{data_dir}/trimmed/{file[0:-len(merged_suffix)]}_trimmed.fastq", 
                         f"{data_dir}/merged/{file}"]
        
        subprocess.check_call(cutadapt_call)

trim_primers(data_dir="../../data/test_data/")