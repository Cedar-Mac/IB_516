import os, subprocess, shutil


def denoise(data_dir:str, DnoisE_args:list=["5", "2", "1", "-y"]):
    """
    
    """

    # Remove old directory and files
    if "denoised" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/denoised/")

    # make output directory for merged files
    os.makedirs(f"{data_dir}/denoised/")

    for file in os.listdir(f"{data_dir}/freq_filtered/"):
        if len(DnoisE_args) == 4:
            DnoisE_call = ["DnoisE.py", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/denoised/{file}",
                            "--alpha", f"{DnoisE_args[0]}",
                            "-x", f"{DnoisE_args[1]}",
                            "-y"]
        if len(DnoisE_args) == 3:
            DnoisE_call = ["DenoisE.py", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/denoised/{file}",
                            "--alpha", f"{DnoisE_args[0]}",
                            "-x", f"{DnoisE_args[1]}",
                            "--min_abund", f"{DnoisE_args[2]}"]

    subprocess.check_call(DnoisE_call)


denoise(data_dir="../../data/test_data")