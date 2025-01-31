import os, subprocess, shutil

def chimera_filter(data_dir:str, vsearch_args:list=["1.4", "8"]):
    """
    
    """

    # Remove old directory and files 
    if "chimera_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/chimera_filtered/")

    # Make directory for quality filtered reads
    os.makedirs(f"{data_dir}/chimera_filtered/")

    for file in os.listdir(f"{data_dir}/length_filtered/"):
        if "fasta" in file:
            vsearch_chimera_call = ["vsearch",
                                 "--uchime_denovo", f"{data_dir}/length_filtered/{file}",
                                 "--nonchimeras", f"{data_dir}/chimera_filtered/{file}",
                                 "--dn", f"{vsearch_args[0]}",
                                 "--xn", f"{vsearch_args[1]}"]

            subprocess.check_call(vsearch_chimera_call)

chimera_filter(data_dir="../../data/test_data")