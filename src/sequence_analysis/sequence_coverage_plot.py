# Import necessary libraries
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

# Function to read A3M file and convert to numpy array
def read_a3m(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return np.array([list(seq) for seq in sequences])

# Read A3M file
msa = read_a3m("your_file.a3m")

# Convert MSA to numerical format (e.g., A=0, C=1, ..., gap=-1)
# This is a simplified example; you may need a more complex mapping
msa_numeric = np.array([[ord(char) for char in seq] for seq in msa])

# Calculate final coverage matrix
final = (msa_numeric != ord('-')).astype(float)

# Plotting
plt.figure(figsize=(14, 4), dpi=100)
plt.subplot(1, 2, 1)
plt.title(f"Sequence coverage ({name})")
plt.imshow(final,
           interpolation='nearest', aspect='auto',
           cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
plt.plot((msa_numeric != ord('-')).sum(0), color='black')
plt.xlim(-0.5, msa_numeric.shape[1] - 0.5)
plt.ylim(-0.5, msa_numeric.shape[0] - 0.5)
plt.colorbar(label="Sequence identity to query")
plt.xlabel("Positions")
plt.ylabel("Sequences")
plt.show()