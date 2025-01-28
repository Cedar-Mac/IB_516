import os, subprocess

def trim_primers(data_dir: str, primer_option:int=1):
    complement = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "Y":"R", "R":"Y", "W":"W", "D":"H", "H":"D"} 
    merged_suffix = "_merged.fastq"

    
    # default primer option uses the EPTDr2n COI primer set
    if primer_option == 1:
        fwd_primer_1 = "GACACTCTTTCCCTACACGACGCTCTTCCGATCTGGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_2 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNGGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_3 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNGGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_4 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNGGDACWGGWTGAACWGTWTAYCCHCC"

        rvs_primer_1 = "CAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCAAACAAATARDGGTATTCGDTY"
        rvs_primer_2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNCAAACAAATARDGGTATTCGDTY"
        rvs_primer_3 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNCAAACAAATARDGGTATTCGDTY"
        rvs_primer_4 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNCAAACAAATARDGGTATTCGDTY"

        # Get reverse complement of the reverse primer for the merged reads
        rv_comp_1 = "".join(complement.get(base, base) for base in reversed(rvs_primer_1))
        rv_comp_2 = "".join(complement.get(base, base) for base in reversed(rvs_primer_2))
        rv_comp_3 = "".join(complement.get(base, base) for base in reversed(rvs_primer_3))
        rv_comp_4 = "".join(complement.get(base, base) for base in reversed(rvs_primer_4))
    
    # make output directory for trimmed files
    if "trimmed" not in os.listdir(data_dir):
        os.makedirs(f"{data_dir}/trimmed/")

    # cutadapt call to trim files in the "merged" subdirectory
    for file in os.listdir(f"{data_dir}/merged/"):
        if ".fastq" in file:
            cutadapt_call = ["cutadapt", 
                         "-g", f"^{fwd_primer_1}", 
                         "-g", f"^{fwd_primer_2}", 
                         "-g", f"^{fwd_primer_3}", 
                         "-g", f"^{fwd_primer_4}", 
                         "-a", f"{rv_comp_1}", 
                         "-a", f"{rv_comp_2}", 
                         "-a", f"{rv_comp_3}", 
                         "-a", f"{rv_comp_4}", 
                         "--discard-untrimmed",
                         "-n", "2",
                         "-o", f"{data_dir}/trimmed/{file[0:-len(merged_suffix)]}_trimmed.fastq", 
                         "--error-rate", "0.1",
                         f"{data_dir}/merged/{file}"]
        
            subprocess.check_call(cutadapt_call)

trim_primers(data_dir="../../data/test_data")