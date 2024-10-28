import argparse
from Bio import PDB
"""
This script extracts the protein sequence from a PDB file and saves it in FASTA format.
"""
def extract_sequence_from_pdb(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    sequences = []
    for model in structure:
        for i, chain in enumerate(model):
            chain_sequence = [] # extract sequence from one chain
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    chain_sequence.append(PDB.Polypeptide.three_to_one(residue.get_resname()))
            sequences.append(''.join(chain_sequence)) #append sequence of one chain to the list of sequences in PDB file
    return sequences

def sequences_to_fasta(sequences: list[str], pdb_file: str, output: str) -> None:
    #get the file name without the extension
    file_name = pdb_file.split('/', 1)[-1].rsplit('.', 1)[0]
    fasta_file = output + '/' + file_name + '.fasta'
    with open(fasta_file, 'w') as f:
        for i, sequence in enumerate(sequences):
            f.write(f">{file_name}_chain{i+1}\n") #Header of the sequence including the file name and chain number
            f.write(sequence + '\n') #Sequence of the chain


def main():
    parser = argparse.ArgumentParser(description='Extract protein sequence from a PDB file and save it in FASTA format.')
    parser.add_argument('--pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('--output', type=str, help='Path to the output FASTA file')
    args = parser.parse_args()

    sequences = extract_sequence_from_pdb(args.pdb_file)
    sequences_to_fasta(sequences, args.pdb_file, args.output)
    print(f"Protein sequence saved to {args.output}/{args.pdb_file.split('/', 1)[-1].split('.', 1)[0] + '.fasta'}")

if __name__ == '__main__':
    main()