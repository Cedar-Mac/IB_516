import subprocess
from Bio import AlignIO
from Bio.Seq import Seq

# Got amino acid sequences from supplementary materials of
# Pentinsaari paper "Molecular evolution of a widely-adopted 
# taxonomic marker (COI) across the animal tree of life"

# Note, find conserved amino acid positions before incorporating
# example sequence from the study. Example sequence could have
# errors at conserved regions. Only use of the example sequence
# is to align the amplicon with the reference alignment.

def align_ref_sequences(data_dir:str, in_file, outfile:str):
    """
    Align amino acid reference sequences from Pentinsaari paper.

    Inputs:
        - data_dir: string, with path to data directory
        - outfile: string, name of output file (no path)

    mafft arguments:
        - Running mafft with no additional arguments
            only the fasta of amino acid sequences.
    
    Output:
        - Aligned fasta file.
    """

    with open(f"{data_dir}/{outfile}", "w+") as handle:
        mafft_call = ["mafft", 
                      f"{data_dir}/{in_file}"]
        subprocess.check_call(mafft_call, stdout = handle)
        print("COI reference sequences aligned.")


def get_conserved_aas(data_dir:str, ref_file:str) -> dict:
    """
    From aligned reference sequences, find conserved amino acids.

    Inputs:
        - data_dir: string, path to data.
        - ref_file: string, name of reference alignment file

    Outputs:
        - Dict: Key = position of conserved amino acid
                Value = Amino acid letter abrr.
    """

    # Load the alignment file
    alignment = AlignIO.read(f"{data_dir}/{ref_file}", "fasta")

    # Initialize a dict to store conserved aa's
    conserved_aas = dict()

    # Iterate over each column in the alignment
    for column in range(0, len(alignment[0,:])):
        # Check if all bases in the column are the same
        if len(set(alignment[:,column])) == 1:
            # If true, add the base to the conserved_bases list
            conserved_aas[column] = set(alignment[:,column]).pop()
    return conserved_aas



def align_amplicon(data_dir, ex_file, aligned_file, fully_aligned_file):
    """
    Align amplicon with reference alignment

    Inputs:
        - data_dir: path to data
        - aligned_file: file with 
        - ex_file: name of file to write example amino acid sequence
        - fully_aligned_file: name of file to write the aligned reference and example sequences.

    Outputs:
        - File with fully aligned reference and example sequences.
    """
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
    with open(fully_aligned_file, "w+") as handle:

        mafft_call = ["mafft", 
                      "--add", f"{data_dir}/{ex_file}", 
                        "--out", f"{data_dir}/{fully_aligned_file}",
                        f"{data_dir}/{aligned_file}"]

        subprocess.check_call(mafft_call, stdout = handle)


def get_conserved_base_positions(data_dir:str, fully_aligned_file:str):
    """
    Get the base positions in target amplicon of conserved amino acids.

    Inputs:
        - data_dir: Path to data
        - fully_aligned_file: Name of aligned fasta with example and refs.

    Outputs:
        - 3 lists containing the base positions of codons 1, 2, and 3.

    Details:
        The three conserved amino acids in the COI amplicon are:
        - 134: H
        - 137: G
        - 146: N
    """
    # AA sequence from the bugs lines up with the 114th AA position in the COI 5' fragment (Folmer fragment)
    # Amplicon is 142 bp, so goes up until AA 162

    aa_1_pos = 134
    aa_2_pos = 137
    aa_3_pos = 146

    # sequence spread across multiple lines.
    # concatenate each line to this empty string.
    aa_seq = ""

    # Only start concatenating lines once we reach the ex sequence
    target_seq = False

    with open(f"{data_dir}/{fully_aligned_file}") as align_handle:
        for line in align_handle.readlines():
            bare_line = line.rstrip()
            if target_seq == True: # once true concatenate lines to make full seq
                aa_seq += bare_line
            if ">ex_amplicon" in bare_line: # Check for header after
                target_seq = True

    # with aa_seq find position of first amino acid (skip all "-")
    for i, let in enumerate(aa_seq):
        seq_start = -1
        if let.isalpha():
            seq_start = i
            break
    
    # Subtract the start position of amplicon
    # from the amino acid position in ref alignment.
    # Multiply by three and we get the codon start position
    # in our amplicon sequence.
    codon_1_start = (aa_1_pos - seq_start) * 3
    codon_2_start = (aa_2_pos - seq_start) * 3
    codon_3_start = (aa_3_pos - seq_start) * 3

    # the base positions of each codon
    codon_1 = [codon_1_start, codon_1_start + 1, codon_1_start + 2]
    codon_2 = [codon_2_start, codon_2_start + 1, codon_2_start + 2]
    codon_3 = [codon_3_start, codon_3_start + 1, codon_3_start + 2]

    return codon_1, codon_2, codon_3