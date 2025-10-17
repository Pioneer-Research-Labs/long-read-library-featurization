import pandas as pd

def read_fasta_to_dataframe(fasta_path):
    """
    Reads a FASTA file and returns a pandas DataFrame with columns:
    'id' (sequence identifier) and 'sequence'.

    Assumes each FASTA header only contains the sequence ID.
    
    Parameters:
        fasta_path (str): Path to the FASTA file.

    Returns:
        pd.DataFrame: DataFrame with FASTA records.
    """
    records = []
    with open(fasta_path, 'r') as f:
        seq_id = None
        sequence_lines = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    records.append({
                        'id': seq_id,
                        'sequence': ''.join(sequence_lines)
                    })
                seq_id = line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(line)

        # Add the last record
        if seq_id:
            records.append({
                'id': seq_id,
                'sequence': ''.join(sequence_lines)
            })

    return pd.DataFrame(records)