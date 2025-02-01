import os, subprocess, shutil


def denoise(data_dir:str, DnoisE_args:list=["5", "2", "1", "-y"]):
    """
    Denoise using Antich's DnoisE algorithm.

    Inputs:
        - data_dir: string with path to data directory.

    DnoisE_args:
        - [0] --alpha: alpha value
        - [1] -x:
        - [2] --min_abund:
        - [3] -y: Use entropy? -y for yes, otherwise no flag.

    Outputs:
        - Denoised fasta files in subdirectory plus csv INFO files.
    """

    # Remove old directory and files
    if "denoised" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/denoised/")

    # make output directory for merged files
    os.makedirs(f"{data_dir}/denoised/")

    for file in os.listdir(f"{data_dir}/freq_filtered/"):
        if len(DnoisE_args) == 4: # With -y flag
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/denoised/{file}",
                            "--alpha", f"{DnoisE_args[0]}",
                            "-x", f"{DnoisE_args[1]}",
                            "--min_abund", f"{DnoisE_args[2]}",
                            "-y"]
            
        if len(DnoisE_args) == 3: # No -y flag
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/denoised/{file}",
                            "--alpha", f"{DnoisE_args[0]}",
                            "-x", f"{DnoisE_args[1]}",
                            "--min_abund", f"{DnoisE_args[2]}"]

        subprocess.check_call(DnoisE_call)


if __name__ == "__main__":
    denoise(data_dir="../../data/test_data")