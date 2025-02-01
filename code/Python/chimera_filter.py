import os, subprocess, shutil

def chimera_filter(data_dir:str, vsearch_args:list=["1.4", "8", "3", "1.2", "0.2"]):
    """
    Input:
        - data_dir: string of data directory.

    vsearch arguments:
        - [0] -dn: Pseudo-count prior for "no" votes. 
                Increasing this value tends to decrease the number of false positives (and also sensitivity)
                Set to default of 1.4
        - [1] -xn: Weight of "no" vote (β).  
                Increasing this value tends to decrease the number of false positives (and also sensitivity). 
                Must be > 1, default is 8.
        - [2] --mindiffs: Minimum number of diffs in a segment. 
                Increasing this value tends to reduce the number of false positives while
                reducing sensitivity to very low-divergence chimeras. Must be > 0.
        - [3] --mindiv: Minimum divergence, i.e. 100% - identity between the query and closest reference database sequence. 
                Expressed as a percentage, so the default is 0.8%, which allows chimeras that are up to 
                99.2% similar to a reference sequence.
        - [4] --minh: Minimum score (h) to be classified as chimera. 
                Increasing this value tends to increase the number of false positives (and also sensitivity).
    
    Output:
        - Chimera filtered subfolder in data directory.
    """

    if len(vsearch_args) != 5:
        print(f"chimera_filter vsearch_args must be list of length 5, {len(vsearch_args)} were provided. Exiting.")
        exit()

    fasta_suffix = ".fasta"

    # Remove old directory and files 
    if "chimera_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/chimera_filtered/")

    # Make directory for quality filtered reads
    os.makedirs(f"{data_dir}/chimera_filtered/")

    os.makedirs(f"{data_dir}/chimera_filtered/logs")

    for file in os.listdir(f"{data_dir}/length_filtered/"):
        if "fasta" in file:
            vsearch_chimera_call = ["vsearch",
                                 "--uchime_denovo", f"{data_dir}/length_filtered/{file}",
                                 "--nonchimeras", f"{data_dir}/chimera_filtered/{file}",
                                 "--dn", f"{vsearch_args[0]}",
                                 "--xn", f"{vsearch_args[1]}",
                                 "--mindiffs", f"{vsearch_args[2]}",
                                 "--mindiv", f"{vsearch_args[3]}",
                                 "--minh", f"{vsearch_args[4]}"]

            with open(f"{data_dir}/chimera_filtered/logs/{file[0:-len(fasta_suffix)]}.log", "a") as log:
                try:
                    subprocess.run(vsearch_chimera_call, stdout=log, stderr=log, check=True)
                    print(f"\nChimera filtered {file} successfully.\n")
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {file}: {e}\n")


if __name__ == "__main__":
    chimera_filter(data_dir="../../data/test_data")