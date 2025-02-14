import os, subprocess, shutil


def denoise(data_dir:str, 
            output_dir:str="denoised", 
            option:str="dnoise", 
            Unoise_args:list=["1", "2"], 
            DnoisE_args:list=["2", "1", "3", "-y"]):
    """
    Denoise using Antich's DnoisE algorithm.

    Inputs:
        - data_dir: Path to data directory.
        - output_dir: Path to output directory.
        - option: "unoise" or "dnoise"

    DnoisE_args:
        - [0] --alpha: alpha value for Unoise distance calculation.
        - [1] -x: First codon position in sequence (3 for Jared's amplicon)
        - [2] --min_abund: Minimum abundance of sequences to include.
                Abundance filtering has already been applied, so set to 1.
        - [3] -y: Use entropy? -y for yes, otherwise no flag.

    Unoise_args:
        - [0] --minsize: Minimum abundance of sequences to include.
                Abundance filtering has already been applied, so set to 1.
        - [1] --unoise_alpha: alpha value for Unoise distance calculation.

    Outputs:
        - Denoised fasta files in subdirectory plus csv INFO files.
    """

    fasta_suffix = ".fasta"

    # Remove old directory and files
    if "denoised" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/{output_dir}/")

    # make output directory for denoised files
    os.makedirs(f"{data_dir}/{output_dir}/")

    # make log subdirectory
    os.makedirs(f"{data_dir}/{output_dir}/logs/")

    for file in os.listdir(f"{data_dir}/freq_filtered/"):

        if option == "dnoise":
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/{output_dir}/{file}",
                            "--alpha", str(DnoisE_args[0]),
                            "--min_abund", str(DnoisE_args[1]),
                            "-x", str(DnoisE_args[2]),
                            "-y"]
            
            if ".fasta" in file:
                with open(f"{data_dir}/{output_dir}/logs/{file[0:-len(fasta_suffix)]}.log", "a") as log:
                    try:
                        subprocess.run(DnoisE_call, stdout=log, stderr=log, check=True)
                        print(f"\nDenoised {file} successfully with DnoisE.\n")
                    except subprocess.CalledProcessError as e:
                        print(f"\nError processing {file}: {e}\n")
            
        if option == "unoise":
            Unoise_call = ["vsearch",
                           "--cluster_unoise", f"{data_dir}/freq_filtered/{file}",
                           "--centroids", f"{data_dir}/{output_dir}/{file}",
                           "--minsize", str(Unoise_args[0]),
                           "--unoise_alpha", str(Unoise_args[1]),
                           "--sizeout",
                           "--sizein"
                           ]

            if ".fasta" in file:
                with open(f"{data_dir}/{output_dir}/logs/{file[0:-len(fasta_suffix)]}.log", "a") as log:
                    try:
                        subprocess.run(Unoise_call, stdout=log, stderr=log, check=True)
                        print(f"\nDenoised {file} successfully with Unoise.\n")
                    except subprocess.CalledProcessError as e:
                        print(f"\nError processing {file}: {e}\n")


if __name__ == "__main__":
    args = [1, 2]
    denoise(data_dir="../../data/test_data", option="unoise", Unoise_args=args)