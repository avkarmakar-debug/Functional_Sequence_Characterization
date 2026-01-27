import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# --- PATH CONFIGURATION ---
# This code automatically finds the path relative to where this script is saved.
# It looks for: ../data/input_sequence.fasta
base_dir = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(base_dir, "..", "data", "input_sequence.fasta")
output_file = os.path.join(base_dir, "..", "results", "qc_summary.txt")

print(f"Reading sequence from: {input_file}")

# Check if file exists
if not os.path.exists(input_file):
    print(f"ERROR: File not found!")
    print("Please make sure 'input_sequence.fasta' is inside the 'data' folder.")
    sys.exit()

# --- STEP 2: QUALITY ANALYSIS ---
try:
    # 1. Read the sequence
    record = SeqIO.read(input_file, "fasta")
    
    # 2. Calculate Length
    seq_len = len(record.seq)
    
    # 3. Calculate GC Content
    gc_percent = gc_fraction(record.seq) * 100
    
    # 4. Generate the Report Text
    report_lines = [
        "=" * 40,
        "SEQUENCE QUALITY CONTROL REPORT",
        "=" * 40,
        f"Sequence ID:   {record.id}",
        f"Description:   {record.description}",
        "-" * 40,
        f"Length:        {seq_len} bp",
        f"GC Content:    {gc_percent:.2f}%",
        "=" * 40
    ]
    
    report_text = "\n".join(report_lines)
    
    # Print to the Terminal at the bottom
    print("\n" + report_text + "\n")
    
    # Save to the 'results' folder
    with open(output_file, "w") as f:
        f.write(report_text)
    print(f"[SUCCESS] Report saved to: results/qc_summary.txt")

except Exception as e:
    print(f"An error occurred: {e}")