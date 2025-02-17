import os, re
import matplotlib.pyplot as plt
import alignments

def check_seqs(data_dir:str, fasta:str, codon_1:list, codon_2:list, codon_3:list)->list:
    """
    Get number of errors and total sequences in a fasta file.

    Inputs:
        - data_dir: path to directory
        - fasta: Name of fasta (no path)
        - codon_1, 2, & 3: Lists of three base positions of each codon.

    Outputs:
        - List containing: 
            [0]: num_errors
            [1]: num_seqs
            [2]: error_rate
    """

    codon_dict = {"H": ["CAC", "CAT"], "G": ["GGA", "GGT", "GGC", "GGG"], "N": ["AAC", "AAT"]}
    stop_seqs = ["TAA", "TGA", "TAG"]

    with open(f"{data_dir}/{fasta}") as file:
        num_errors = 0
        num_seqs = 0
        for rawline in file.readlines():
            num_seqs += 1
            if (rawline.strip()) and (">" not in rawline):
                # Two N's added to match reading frame
                seq = "NN" + rawline.rstrip()
                first_codon = "".join([seq[i] for i in codon_1])
                second_codon = "".join([seq[i] for i in codon_2])
                third_codon = "".join([seq[i] for i in codon_3])

                if ((first_codon != codon_dict["H"][0]) and 
                    (first_codon != codon_dict["H"][1])):
                    num_errors += 1

                if ((second_codon != codon_dict["G"][0]) and 
                    (second_codon != codon_dict["G"][1]) and
                    (second_codon != codon_dict["G"][2]) and
                    (second_codon != codon_dict["G"][3])):
                    num_errors += 1

                if ((third_codon != codon_dict["N"][0]) and 
                    (third_codon != codon_dict["N"][1])):
                    num_errors += 1
                
                for codon in [seq[i:i+3] for i in range(0, len(seq), 3)]:
                    if (codon == stop_seqs[0] or 
                        codon == stop_seqs[1] or 
                        codon == stop_seqs[2]):
                        num_errors += 1

        error_rate = num_errors / num_seqs

        return [num_errors, num_seqs, error_rate]
    

def get_method_avgs(denoised_dir, o, a, c1, c2, c3):
    method_stats = []
    for file in os.listdir(f"{denoised_dir}/{o}_alpha_{a}"):

        if (".csv" not in file) and ("log" not in file):
            file_stats = check_seqs(data_dir=f"{denoised_dir}/{o}_alpha_{a}", 
                                                     fasta=file, 
                                                     codon_1=c1, 
                                                     codon_2=c2, 
                                                     codon_3=c3
                )
            method_stats.append(file_stats)

    err_counts = [site[0] for site in method_stats]
    seq_counts = [site[1] for site in method_stats]
    avg_err_count = sum(err_counts) / len(err_counts)
    avg_seq_counts = sum(seq_counts) / len(seq_counts)
    avg_err_rate = avg_err_count / avg_seq_counts

    return avg_err_rate, avg_seq_counts, avg_err_count


def make_plot(stat, stats):
    fig, ax = plt.subplots()

    x = [int(*re.findall(r'\d+', item[0])) for item in stats if "dnoise" in item[0]]

    if stat == "error rate":
        y1 = [item[1][0] for item in stats if "dnoise" in item[0]]
        y2 = [item[1][0] for item in stats if "unoise" in item[0]]
    if stat == "sequence count":
        y1 = [item[1][1] for item in stats if "dnoise" in item[0]]
        y2 = [item[1][1] for item in stats if "unoise" in item[0]]
    if stat == "error count":
        y1 = [item[1][2] for item in stats if "dnoise" in item[0]]
        y2 = [item[1][2] for item in stats if "unoise" in item[0]]

    ax.plot(x, y1, label="Entropy Correction")
    ax.plot(x, y2, linestyle="dashed", label="No Entropy")
    ax.set_title(f"{stat}")
    ax.set_xlabel("alpha value")
    ax.set_ylabel(f"{stat}")
    ax.legend()
    plt.savefig(f"../../tmp_plots/{stat}.png")
    print(f"{stat} plot saved in tmp_plots directory.")

if __name__ == "__main__":

    alignment_dir = "../../data/alignments"
    denoised_dir = "../../data/test_data/benchmarking"
    ref_file = "coi_fixed.fasta"
    ref_aligned = "coi_aligned.fasta"
    fully_aligned = "bug_aligned_coi.fasta"
    file_suffix = ".fasta_Adcorr_denoised_d.fasta"

    
    c_1, c_2, c_3 = alignments.get_conserved_base_positions(alignment_dir, fully_aligned_file = fully_aligned)

    all_stats = []
    for option in ["dnoise", "unoise"]:
        for alpha in [1,3,5,7,9]:
            avg_err_rate, avg_seq_counts, avg_err_counts = get_method_avgs(denoised_dir=denoised_dir,
                                                                    o=option,
                                                                    a=alpha,
                                                                    c1=c_1,
                                                                    c2=c_2,
                                                                    c3=c_3)
            all_stats.append((f"{option}_alpha_{alpha}", [avg_err_rate, avg_seq_counts, avg_err_counts]))

    make_plot("error rate", all_stats)
    

