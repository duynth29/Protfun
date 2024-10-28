import pandas as pd
from Bio import SeqIO

def fasta_files_to_csv(fasta_files: list[str], output_name: str = "", write_csv: bool = False) -> pd.DataFrame:
    """
    Convert a list of fasta files to a csv file as the input for ColabFold
    
    Args:
        fasta_files: list of fasta files
        output_name: name of the csv file
        write_csv: if True, write the csv file
    Returns:
        pd.DataFrame: a dataframe with two columns: id and sequence
    csv file format:
    id, sequence
    prot_a, chain_1:chain_2:chain_n
    prot_b, chain_1:chain_2:chain_n
    ...
    """
    protein_ids = []
    protein_sequences = []
    
    for fasta_file in fasta_files:
        # Extract the protein id from the fasta file
        protein_id = fasta_file.split('/', 1)[-1].split('.', 1)[0]
        protein_ids.append(protein_id)
        
        # Parse the fasta file to get the sequences
        sequences = SeqIO.parse(fasta_file, "fasta")
        sequence_list = [str(sequence.seq) for sequence in sequences]
        protein_sequences.append(':'.join(sequence_list))
    
    df = pd.DataFrame({'id': protein_ids, 'sequence': protein_sequences})
    
    if write_csv:
        df.to_csv(output_name, index=False)
    
    return df
