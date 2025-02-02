import os, subprocess, shutil


def denoise(data_dir:str, output_dir="denoised", DnoisE_args:list=["1", "5", "3", "1", "1", "-y"]):
    """
    Denoise using Antich's DnoisE algorithm.

    Inputs:
        - data_dir: string with path to data directory.

    DnoisE_args:
        - [0] --joining_criteria:
        - [1] --alpha: alpha value
        - [2] -x:
        - [3] --min_abund:
        - [4] --cores: 
        - [5] -y: Use entropy? -y for yes, otherwise no flag.

    Outputs:
        - Denoised fasta files in subdirectory plus csv INFO files.
    """

    if (len(DnoisE_args) < 5) or (len(DnoisE_args) > 6):
        print(f"denoise DnoisE_args must be list of length 5 or 6, length {len(DnoisE_args)} provided. Exiting.")
        exit()

    fasta_suffix = ".fasta"

    # Remove old directory and files
    if "denoised" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/{output_dir}/")

    # make output directory for denoised files
    os.makedirs(f"{data_dir}/{output_dir}/")

    # make log subdirectory
    os.makedirs(f"{data_dir}/{output_dir}/logs/")

    for file in os.listdir(f"{data_dir}/freq_filtered/"):
        if len(DnoisE_args) == 6: # With -y flag
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/{output_dir}/{file}",
                            "-j", str(DnoisE_args[0]),
                            "--alpha", str(DnoisE_args[1]),
                            "-x", str(DnoisE_args[2]),
                            "--min_abund", str(DnoisE_args[3]),
                            "-c", str(DnoisE_args[4]),
                            "-y"]
            
        if len(DnoisE_args) == 5: # No -y flag
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/{output_dir}/{file}",
                            "-j", str(DnoisE_args[0]),
                            "--alpha", str(DnoisE_args[1]),
                            "-x", str(DnoisE_args[2]),
                            "--min_abund", str(DnoisE_args[3]),
                            "-c", str(DnoisE_args[4])]

        if ".fasta" in file:
            with open(f"{data_dir}/{output_dir}/logs/{file[0:-len(fasta_suffix)]}.log", "a") as log:
                try:
                    subprocess.run(DnoisE_call, stdout=log, stderr=log, check=True)
                    print(f"\nDenoised {file} successfully.\n")
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {file}: {e}\n")


if __name__ == "__main__":

    args = ["1", "9", "3", "1", "1"]
    denoise(data_dir="../../data/test_data", DnoisE_args=args)