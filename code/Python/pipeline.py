import os, re, subprocess, shutil
from itertools import pairwise
from collections import Counter


###### MERGE ######
def merge_pairs(data_dir:str, vsearch_args:list=["99", "16", "25", "--fastq_allowmergestagger"]):
    """
    Merges paired fastq reads using the vsearch algorithm.

    Input:

        - data_dir: string providing the path to the data directory with raw input fastq files.

    vsearch arguments:

        - [0] --fastq_mergepairs: specify the maximum number of non-matching nucleotides 
                            allowed in the overlap region. That option has a strong 
                            influence on the merging success rate. Set to 99 a la JAMP.
        - [1] --fastq_minovlen: Minimum overlap between paired sequences to merge. JAMP default is 16,
                            base default is 10. Currently set to JAMP default.
        - [2] --fastq_diffpct: Maximum allowed percent difference between the two paired sequences. 
                            vsearch recommended default is 100%, JAMP sets the default to 25%
                            I'm using the JAMP default.
        - [3] --fastq_allowmergestagger: Allows for reads to be staggered, this is following JAMP convention.

    Output:

        - Merged fastq files in merged directory within data directory.
    """

    if len(vsearch_args) != 4:
        print(f"merge_pairs vsearch_args must be list of length 4, {len(vsearch_args)} were provided. Exiting.")
        exit()

    data_files = os.listdir(data_dir)

    fastq_suffix = "_RX.fastq"

    # get file names of forward and reverse reads, strip suffix
    fwd_files = [fname for fname in data_files if re.search(pattern = "_R1.fastq", string = fname)]
    rvs_files = [fname for fname in data_files if re.search(pattern = "_R2.fastq", string = fname)]

    # Check that all reads are paired 
    fwd_bare_names = set([name[0:-len(fastq_suffix)] for name in fwd_files])
    rvs_bare_names = set([name[0:-len(fastq_suffix)] for name in rvs_files])

    if fwd_bare_names != rvs_bare_names:
        print("Not all reads are paired!")
        exit()

    # Remove old directory and files
    if "merged" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/merged/")

    # make output directory for merged files
    os.makedirs(f"{data_dir}/merged/")

    #Make output directory for logs
    if "logs" not in os.listdir(f"{data_dir}/merged/"):
        os.makedirs(f"{data_dir}/merged/merge_logs/") 

    # Use VSEARCH fastq_mergepairs function
    for name in fwd_bare_names:

        vsearch_merge_call = ["vsearch", 
                                "--fastq_mergepairs", f"{data_dir}/{name}_R1.fastq", 
                                "--reverse", f"{data_dir}/{name}_R2.fastq", 
                                "--fastqout", f"{data_dir}/merged/{name}_merged.fastq",
                                "--fastq_maxdiffs", f"{vsearch_args[0]}", 
                                "--fastq_minovlen", f"{vsearch_args[1]}",
                                "--fastq_maxdiffpct", f"{vsearch_args[2]}",
                                f"{vsearch_args[3]}"]

        with open(f"{data_dir}/merged/merge_logs/{name}.log", "a") as log:
            try:
                subprocess.run(vsearch_merge_call, stdout=log, stderr=log, check=True)
                print(f"\nMerged {name} successfully.\n")
            except subprocess.CalledProcessError as e:
                print(f"\nError processing {name}: {e}\n")


###### TRIM ######
def trim_primers(data_dir: str, primer_option:int=1):
    """
    Inputs:
        - data_dir: String specifying data directory
        - primer_option: Integer specifying the primer set to use
            - 1 = invertebrate COI
            - Anything else: not supported yet
        
    cutadapt arguments:
        - [0] -g: The 5' primer sequence. I anchor the fwd primer with ^
        - [1] -a: The 3' primer sequence. Not anchored
        - [2] --discard_untrimmed: Get rid of sequences with no matching primer
        - [3] -n: Number of times to repeat. Cutadapt only removes one primer at a time.
            Repeat twice to remove fwd and rvs primers.
        - [4] -o: output file
        - [5] --error-rate: What fraction of bases can not match the primer sequence.
                        Set to 0.1 in line with JAMP.
    
    Outputs:
        - Trimmed sequences output in sub-folder of data directory in fastq format.
    """

    # Complementary base pairings for making the reverse complement of a primer sequence.
    # Necessary for the reverse primer after merging paired-end reads
    complement = {"A":"T", 
                  "T":"A",
                  "C":"G",
                  "G":"C",
                  "N":"N",
                  "Y":"R",
                  "R":"Y",
                  "W":"W",
                  "D":"H",
                  "H":"D"} 
    
    merged_suffix = "_merged.fastq"

    # default primer option uses the EPTDr2n COI primer set.
    # Check with Jared about which are necessary (seems like only main primers required for filtering)
    if primer_option == 1:
        fwd_main = "GGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_1 = "GACACTCTTTCCCTACACGACGCTCTTCCGATCTGGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_2 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNGGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_3 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNGGDACWGGWTGAACWGTWTAYCCHCC"
        fwd_primer_4 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNGGDACWGGWTGAACWGTWTAYCCHCC"

        rvs_main = "CAAACAAATARDGGTATTCGDTY"
        rvs_primer_1 = "CAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCAAACAAATARDGGTATTCGDTY"
        rvs_primer_2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNCAAACAAATARDGGTATTCGDTY"
        rvs_primer_3 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNCAAACAAATARDGGTATTCGDTY"
        rvs_primer_4 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNCAAACAAATARDGGTATTCGDTY"

        # Get reverse complement of the reverse primer for the merged reads
        rv_comp_main = "".join(complement.get(base, base) for base in reversed(rvs_main))
        rv_comp_1 = "".join(complement.get(base, base) for base in reversed(rvs_primer_1))
        rv_comp_2 = "".join(complement.get(base, base) for base in reversed(rvs_primer_2))
        rv_comp_3 = "".join(complement.get(base, base) for base in reversed(rvs_primer_3))
        rv_comp_4 = "".join(complement.get(base, base) for base in reversed(rvs_primer_4))
    
    # Remove old directory and files
    if "trimmed" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/trimmed/")

    # make output directory for trimmed files
    os.makedirs(f"{data_dir}/trimmed/")

    # make log subdirectory
    os.makedirs(f"{data_dir}/trimmed/logs")

    # cutadapt call to trim files in the "merged" subdirectory
    for file in os.listdir(f"{data_dir}/merged/"):
        if ".fastq" in file:
            cutadapt_call = ["cutadapt", 
                         "-g", f"^{fwd_main}", 
                         "-g", f"^{fwd_primer_1}", 
                         "-g", f"^{fwd_primer_2}", 
                         "-g", f"^{fwd_primer_3}", 
                         "-g", f"^{fwd_primer_4}", 
                         "-a", f"{rv_comp_main}",
                         "-a", f"{rv_comp_1}", 
                         "-a", f"{rv_comp_2}", 
                         "-a", f"{rv_comp_3}", 
                         "-a", f"{rv_comp_4}", 
                         "--discard-untrimmed",
                         "-n", "2",
                         "-o", f"{data_dir}/trimmed/{file[0:-len(merged_suffix)]}_trimmed.fastq", 
                         "--error-rate", "0.1",
                         f"{data_dir}/merged/{file}"]
        
            with open(f"{data_dir}/trimmed/logs/{file[0:-len(merged_suffix)]}.log", "a") as log:
                try:
                    subprocess.run(cutadapt_call, stdout=log, stderr=log, check=True)
                    print(f"\nTrimmed {file} successfully.\n")
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {file}: {e}\n")


###### QUAL. FILTER ######
def quality_filter(data_dir:str, vsearch_args:list=[1, 0]):
    """
    Filter fastq filters that have already been trimmed.

    Input:
        - data_dir: String with main data directory

    vsearch_args:
        - [0] --fastq_maxee: Max expected cummulative error.
                            Set to 1 (probability of single error = )
        - [1] --fastq_maxns: Number of allowed N's in the sequence.
                            Set to 0

    Output:
        - quality filtered reads in subfolder of data directory
    """

    if len(vsearch_args) != 2:
        print(f"quality_filter vsearch_args must be list of length two. {len(vsearch_args)} arguments were provided. Exiting.")
        exit()

    trimmed_suffix = "_trimmed.fastq"

    # Remove old directory and files 
    if "quality_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/quality_filtered/")

    # Make directory for quality filtered reads
    os.makedirs(f"{data_dir}/quality_filtered/")

    # Make subdirectory for logs
    os.makedirs(f"{data_dir}/quality_filtered/logs/")

    for file in os.listdir(f"{data_dir}/trimmed/"):
        # vsearch call for fastq filtering. See function description for arguments
        if ".fastq" in file:
            vsearch_ee_filter_call = ["vsearch",
                                      "--fastx_filter", f"{data_dir}/trimmed/{file}",
                                      "--fastqout", f"{data_dir}/quality_filtered/{file[0:-len(trimmed_suffix)]}.fastq",
                                      "--fastq_maxee", f"{vsearch_args[0]}",
                                      "--fastq_maxns", f"{vsearch_args[1]}"]
        
            with open(f"{data_dir}/quality_filtered/logs/{file[0:-len(trimmed_suffix)]}.log", "a") as log:
                try:
                    subprocess.run(vsearch_ee_filter_call, stdout=log, stderr=log, check=True)
                    print(f"\nQuality filtered {file} successfully.\n")
                except subprocess.CalledProcessError as e:
                    print(f"\nError processing {file}: {e}\n")


###### LENGTH FILTER ######
def length_filter(data_dir:str, amplicon_length:int):
    """
    Filter sequences to match amplicon length.

    Inputs:
        - data_dir: path to data directory as a string
        - amplicon_length: fixed length of amplicon sequence.

    Outputs:
        - Trimmed sequences as fasta files in subfolder of data directory
    """

    fastq_suffix = ".fastq"

    # Remove old files and make output directory for length-filtered files
    if "length_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/length_filtered/")

    # make output directory
    os.makedirs(f"{data_dir}/length_filtered/")

    # make log subdirectory
    os.makedirs(f"{data_dir}/length_filtered/log")


    for data_file in os.listdir(f"{data_dir}/quality_filtered/"):
        if "fastq" in data_file:
            rm_counts = 0 # Keep track of removed lines (includes metadata lines at the moment)
            kepper_counts = 0 # Keep track of kept sequences

            log_file = open(f"{data_dir}/length_filtered/log/length_filter.log", "a")
            in_file = open(f"{data_dir}/quality_filtered/{data_file}", "r")
            out_file = open(f"{data_dir}/length_filtered/{data_file[0:-len(fastq_suffix)]}.fasta", "a")

            lines = in_file.readlines()

            # Iterate over pairs of lines.
            for i, pair in enumerate(pairwise(lines)):
                if (i % 4 == 0) and (len(pair[1].rstrip()) == amplicon_length): # if current line is a header and next line is a keeper sequence
                    out_file.write(f">{pair[0]}") # write header line
                if i % 4 == 1:
                    if len(pair[0].rstrip()) == amplicon_length: # if current line is a keeper sequence
                        out_file.write(f"{pair[0]}") # write keeper sequence
                        kepper_counts += 1
                    if len(pair[0]) != 142: # sequence was not correct length, keep track of skipped sequences.
                        rm_counts += 1
            out_file.close()
            in_file.close()
            log_file.write(f"{data_file} had {rm_counts} sequences removed and {kepper_counts} sequences were kept.\n\n")
    log_file.close()


###### CHIMERA FILTER ######
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


###### FREQ. FILTER ######
def frequency_filter(data_dir:str, min_seq_count:int, min_site_occurance:int):
    """
    Filter sequences based on their frequency of occurances within and between sites.

    Input:
        - data_dir: string with path to data directory
        - min_seq_count: Integer, minimum sequence occurance at a site
        - min_site_occurance: minimum number of sites a low abundance sequence
                            must occur at to be retained.

    Output: 
        - fasta files with size (seq count) in header in freq_filtered subdirectory.
        - Formated for input into DnoisE filtering algorithm.

    Details:
        A sequence occurance at a site is only dropped if two conditions are met:

        - The sequence occurs less than a minimum number of times (min_seq_count)
        - And the sequence occurs at less than a minimum number of sites (min_site_occurance).

        A temporary file for each site is created to group duplicate sequences and order them by abundance.
        Currently the only information stored in the sequence header is the sequence name (site_name_sequence_id)
        and count of the sequence ("size")
    """

    fasta_suffix = ".fasta"

    # Remove old directory and files
    if "freq_filtered" in os.listdir(data_dir):
        shutil.rmtree(f"{data_dir}/freq_filtered/")

    # make output directory for frequency filtered files
    os.makedirs(f"{data_dir}/freq_filtered/")

    # Make subdirectory for log
    os.makedirs(f"{data_dir}/freq_filtered/log/")

    # make temp directory for sorted counts
    os.makedirs(f"{data_dir}/freq_filtered/temp/")
    os.makedirs(f"{data_dir}/freq_filtered/fixed_infiles/")

    too_few_list = []

    for file in os.listdir(f"{data_dir}/chimera_filtered/"):

        if "fasta" in file:
            # Sequence getting split over two lines. Make sequences whole.
            with open(f"{data_dir}/chimera_filtered/{file}", "r") as infile, open(f"{data_dir}/freq_filtered/fixed_infiles/{file}", "w") as outfile:
                for line in infile:
                    if line.startswith(">"):
                        outfile.write(f"\n{line}")
                    else:
                        outfile.write(line.strip())

            temp_file = open(f"{data_dir}/freq_filtered/temp/{file}", "a")
            in_file = open(f"{data_dir}/freq_filtered/fixed_infiles/{file}", "r")

            # Read lines in file into list
            lines = in_file.readlines()

            # Count occurance of sequences
            seq_counts = Counter(lines[2::2]) # my fasta file fix introduces an empty first line

            for key in seq_counts:
                if seq_counts[key] < min_seq_count:
                    too_few_list.append(key)

            seq_counts_sorted = seq_counts.most_common()

            for seq in seq_counts_sorted:
                temp_file.write(f"{seq[1]}\n")
                temp_file.write(f"{seq[0]}")

            temp_file.close()
            in_file.close()
    
    # Ones list now contains all sequences seen a minimum number of times.
    # From these limited sequences, find ones that only occur at a minimum number of sites.
    site_occurances = Counter(too_few_list)

    # Track sequence index
    i = 1

    for tmp_file in os.listdir(f"{data_dir}/freq_filtered/temp/"):
        
        if "fasta" in tmp_file:

            with open(f"{data_dir}/freq_filtered/temp/{tmp_file}") as f:
                out_file = open(f"{data_dir}/freq_filtered/{tmp_file}", "a")
                log_file = open(f"{data_dir}/freq_filtered/log/freq_filter.log", "a")

                for count, seq in zip(f, f):
                    if (site_occurances[seq] >= min_site_occurance) or (int(count) >= min_seq_count):
                        out_file.write(f">id={tmp_file[0:-len(fasta_suffix)]}_{i}; size={count.rstrip()};\n")
                        out_file.write(seq)
                        i += 1
                    else:
                        log_file.write(f"{seq} occured less than {min_seq_count} times at {tmp_file[0:-len(fasta_suffix)]} and present at {site_occurances[seq] - 1} other sites.\n")

                out_file.close()
    log_file.close()


    # Remove temp directories and files
    if "temp" in os.listdir(f"{data_dir}/freq_filtered"):
        shutil.rmtree(f"{data_dir}/freq_filtered/temp/")

    if "fixed" in os.listdir(f"{data_dir}/freq_filtered"):
        shutil.rmtree(f"{data_dir}/freq_filtered/fixed_infiles/") 


###### DENOISE ######
def denoise(data_dir:str, output_dir="denoised", DnoisE_args:list=["1", "5", "2", "1", "10", "-y"]):
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
        if len(DnoisE_args) == 7: # With -y flag
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/{output_dir}/{file}",
                            "-j", str(DnoisE_args[0]),
                            "--alpha", str(DnoisE_args[1]),
                            "-x", str(DnoisE_args[2]),
                            "--min_abund", str(DnoisE_args[3]),
                            "-c", str(DnoisE_args[4]),
                            "-y"]
            
        if len(DnoisE_args) == 6: # No -y flag
            DnoisE_call = ["dnoise", 
                            "--fasta_input", f"{data_dir}/freq_filtered/{file}",
                            "--fasta_output", f"{data_dir}/{output_dir}/{file}",
                            "--j", str(DnoisE_args[0]),
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