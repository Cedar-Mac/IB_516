import requests, subprocess, os
from Bio import AlignIO
from Bio.Seq import Seq

def get_animalia_COI(data_dir:str, params:list=["COI-5P", "Animalia"]): 
    # Define API URL for COI-5P sequences in all animals
    url = "http://www.boldsystems.org/index.php/API_Public/sequence"
    params = {
        "marker": str(params[0]),  # Gene marker for COI
        "taxon": str(params[1]),  # Get all animal sequences
        "format": "fasta"     # Request FASTA format
    }

    with requests.get(url, params=params, timeout=120, stream=True) as response:
        response.raise_for_status()
        with open(f"{data_dir}/COI_sequences.fasta", "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):  # Stream in chunks
                f.write(chunk)
        print("COI sequences downloaded in chunks!")

#get_animalia_COI(data_dir="../../data")





# Load the alignment file
alignment = AlignIO.read("../../data/coi_aligned.fasta", "fasta")

# Initialize a list to store conserved aa's
conserved_aas = dict()

# Iterate over each column in the alignment
for column in range(0, len(alignment[0,:])):
    # Check if all bases in the column are the same
    if len(set(alignment[:,column])) == 1:
        # If true, add the base to the conserved_bases list
        conserved_aas[column] = set(alignment[:,column]).pop()

# Print the conserved bases
print(conserved_aas)



# Example sequence of an amplicon from study
nucleotide_sequence = "NNCTTATCATCTGGAATTGCCCATGCCGGGGCTTCCGTTGATTTAGCAATTTTTTCACTTCACCTAGCAGGAATTTCATCTATTCTAGGGGCTGTAAATTTTATTACTACAATTATTAATATACGATCTAATGGAATTACTTTT"

coding_dna = Seq(nucleotide_sequence)

# Translate to AA sequence
amino_acid_sequence = coding_dna.translate()

# Create fasta file with this single AA sequence for alignment to ref align file (coi_aligned.fasta)
print(amino_acid_sequence)

# AA sequence from the bugs lines up with the 118th AA position in the COI 5' fragment (Folmer fragment)
# Amplicon is 142 bp, so goes up until AA 166
# The 3 conserved AA's in our amplicon are:
# 134: H
# 137: G
# 146: N

codon_dict = {"H": ["CAC", "CAT"], "G": ["GGA", "GGT", "GGC", "GGG"], "N": ["AAC", "AAT"]}
stop_seqs = ["TAA", "TGA", "TAG"]

