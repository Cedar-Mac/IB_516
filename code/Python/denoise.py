import os, subprocess, shutil


def denoise(data_dir:str, 
            output_dir:str, 
            option:str="dnoise", 
            Unoise_args:list=["1", "2"], 
            DnoisE_args:list=["2", "1", "3", "-y"]):
    """
    Denoise using Antich's DnoisE algorithm or Edgar's Unoise3 algorithm.

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
    if os.path.exists(output_dir):
        shutil.rmtree(f"{output_dir}/")

    # make output directory for denoised files
    os.makedirs(f"{output_dir}/")

    # make log subdirectory
    os.makedirs(f"{output_dir}/logs/")

    for file in os.listdir(f"{data_dir}/"):

        if option == "dnoise":
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/{file}",
                            "--fasta_output", f"{output_dir}/{file}",
                            "--alpha", str(DnoisE_args[0]),
                            "--min_abund", str(DnoisE_args[1]),
                            "-x", str(DnoisE_args[2]),
                            "-y"]
            
            if ".fasta" in file:
                with open(f"{output_dir}/logs/{file[0:-len(fasta_suffix)]}.log", "a") as log:
                    try:
                        subprocess.run(DnoisE_call, stdout=log, stderr=log, check=True)
                        print(f"\nDenoised {file} successfully with DnoisE.\n")
                    except subprocess.CalledProcessError as e:
                        print(f"\nError processing {file}: {e}\n")
            
        if option == "unoise":
            Unoise_call = ["usearch",
                           "--unoise3", f"{data_dir}/{file}",
                           "--ampout", f"{output_dir}/{file}",
                           "--minsize", str(Unoise_args[0]),
                           "--unoise_alpha", str(Unoise_args[1])
                           ]

            if ".fasta" in file:
                with open(f"{output_dir}/logs/{file[0:-len(fasta_suffix)]}.log", "a") as log:
                    try:
                        subprocess.run(Unoise_call, stdout=log, stderr=log, check=True)
                        print(f"\nDenoised {file} successfully with Unoise.\n")
                    except subprocess.CalledProcessError as e:
                        print(f"\nError processing {file}: {e}\n")
                
    # DnoisE fasta output gets stupid names, fix
    stupid_suffix = ".fasta_Adcorr_denoised_ratio_d.fasta"
    if option == "dnoise":
        for file in os.listdir(output_dir):
            if "denoised" in file:
                os.rename(os.path.join(output_dir, file), 
                          os.path.join(output_dir, f"{file[0:-len(stupid_suffix)]}_denoised.fasta"))
    
    # Unoise Sequence getting split over two lines. Make sequences whole.
    if option == "unoise":
        for file in os.listdir(output_dir):
            if os.path.isfile(os.path.join(output_dir, file)):
                with open(os.path.join(output_dir, file), "r") as file_io:
                    with open(f"{output_dir}/{file[0:-len(fasta_suffix)]}_denoised.fasta", "w") as outfile:
                        for line in file_io:
                            if line.startswith(">"):
                                outfile.write(f"\n{line}")
                            else:
                                outfile.write(line.strip())
                os.remove(os.path.join(output_dir, file))


if __name__ == "__main__":
    args = [1, 2]
    denoise(data_dir="../../data/test_data/freq_filtered", 
            output_dir="../../data/test_data/denoised",
            option="unoise", 
            Unoise_args=args)