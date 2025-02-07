import subprocess, os
from Bio import AlignIO
from Bio.Seq import Seq

# Got amino acid sequences from supplementary materials of
# Pentinsaari paper "Molecular evolution of a widely-adopted 
# taxonomic marker (COI) across the animal tree of life"

def align_sequences(data_dir, outfile):

    with open(f"{data_dir}/{outfile}", "w+") as handle:
        mafft_call = ["mafft", 
                      f"{data_dir}/coi_fixed.fasta"]
        subprocess.check_call(mafft_call, stdout = handle)
        print("COI reference sequences aligned.")


def get_conserved_aas(data_dir):

    # Load the alignment file
    alignment = AlignIO.read(f"{data_dir}/coi_aligned.fasta", "fasta")

    # Initialize a list to store conserved aa's
    conserved_aas = dict()

    # Iterate over each column in the alignment
    for column in range(0, len(alignment[0,:])):
        # Check if all bases in the column are the same
        if len(set(alignment[:,column])) == 1:
            # If true, add the base to the conserved_bases list
            conserved_aas[column] = set(alignment[:,column]).pop()
    return conserved_aas



def align_amplicon(data_dir, ex_file, aligned_file):
    # Example sequence of an amplicon from study
    nucleotide_sequence = "NNCTTATCATCTGGAATTGCCCATGCCGGGGCTTCCGTTGATTTAGCAATTTTTTCACTTCACCTAGCAGGAATTTCATCTATTCTAGGGGCTGTAAATTTTATTACTACAATTATTAATATACGATCTAATGGAATTACTTTT"

    # Turn into BioPython Seq object
    coding_dna = Seq(nucleotide_sequence)

    # Translate to amino acid sequence
    amino_acid_sequence = coding_dna.translate()

    # Create fasta file of single AA example sequence
    with open(f"{data_dir}/{ex_file}", "w") as out_file:
        out_file.write(">ex_amplicon_amino_acid_seq\n")
        out_file.write(str(amino_acid_sequence))
    
    # Align example fasta from study with reference alignment.
    with open(aligned_file, "w+") as handle:

        mafft_call = ["mafft", 
                      "--add", f"{data_dir}/{ex_file}", 
                        "--out", f"{data_dir}/{aligned_file}",
                        f"{data_dir}/coi_aligned.fasta", ]

        subprocess.check_call(mafft_call, stdout = handle)


def get_conserved_base_positions(data_dir, alignment_file):

    # AA sequence from the bugs lines up with the 114th AA position in the COI 5' fragment (Folmer fragment)
    # Amplicon is 142 bp, so goes up until AA 162
    # The 3 conserved AA's in our amplicon are:
    # 134: H
    # 137: G
    # 146: N
    aa_1_pos = 134
    aa_2_pos = 137
    aa_3_pos = 146
    amplicon_length = 142 + 2 # plus two for two N's added to match reading frame

    aa_seq = ""
    target_seq = False
    with open(f"{data_dir}/{alignment_file}") as align_handle:
        for line in align_handle.readlines():
            bare_line = line.rstrip()
            if target_seq == True:
                aa_seq += bare_line
            if ">ex_amplicon" in bare_line:
                target_seq = True

    for i, let in enumerate(aa_seq):
        seq_start = -1
        if let.isalpha():
            seq_start = i
            break
    
    codon_1_start = (aa_1_pos - seq_start) * 3
    codon_2_start = (aa_2_pos - seq_start) * 3
    codon_3_start = (aa_3_pos - seq_start) * 3

    codon_1 = [codon_1_start, codon_1_start + 1, codon_1_start + 2]
    codon_2 = [codon_2_start, codon_2_start + 1, codon_2_start + 2]
    codon_3 = [codon_3_start, codon_3_start + 1, codon_3_start + 2]

    





codon_dict = {"H": ["CAC", "CAT"], "G": ["GGA", "GGT", "GGC", "GGG"], "N": ["AAC", "AAT"]}
stop_seqs = ["TAA", "TGA", "TAG"]



if __name__ == "__main__":
   # align_sequences(data_dir = "../../data/alignments", outfile="coi_aligned.fasta")
    conserved_aas =  get_conserved_aas(data_dir="../../data/alignments")
    print(conserved_aas)
   # align_amplicon(data_dir = "../../data/alignments", 
   #                ex_file = "coi_bug.fasta", 
   #                aligned_file = "bug_aligned_coi.fasta")
    get_conserved_base_positions("../../data/alignments", alignment_file = "bug_aligned_coi.fasta")