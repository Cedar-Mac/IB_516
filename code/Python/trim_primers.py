import subprocess

def trim_primers(merged_dir: str, primer_option:int):
    fwd_primer = ""
    rvs_primer = ""
    
    if primer_option == 1:
        fwd_primer = "GGDACWGGWTGAACWGTWTAYCCHCC"
        rvs_primer = "CAAACAAATARDGGTATTCGDTY"
    
    for file in merged_dir:
        cutadapt_call = [""]