import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# --- PATH CONFIGURATION ---
base_dir = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(base_dir, "..", "data", "input_sequence.fasta")
output_file = os.path.join(base_dir, "..", "results", "filtering_log.txt")

if not os.path.exists(input_file):
    print("ERROR: Input file not found.")
    sys.exit()

# --- STEP 3: FILTERING & VALIDATION ---
print("Running Step 3: Sequence Filtering & Validation...")

record = SeqIO.read(input_file, "fasta")
seq_len = len(record.seq)
gc_percent = gc_fraction(record.seq) * 100

# DEFINING THE RULES (Criteria)
# Rule 1: Sequence must be longer than 200 bp (to avoid fragments)
# Rule 2: GC content must be between 30% and 70% (normal for coding genes)
min_length = 200
min_gc = 30
max_gc = 70

decision = "FAIL"
reason = ""

if seq_len < min_length:
    reason = f"Too short (Length {seq_len} < {min_length})"
elif gc_percent < min_gc or gc_percent > max_gc:
    reason = f"Abnormal GC Content ({gc_percent:.2f}%)"
else:
    decision = "PASS"
    reason = "Sequence quality is within normal biological ranges."

# --- REPORT ---
report_text = f"""
========================================
STEP 3: FILTERING DECISION REPORT
========================================
Sequence ID: {record.id}
Length:      {seq_len} bp
GC Content:  {gc_percent:.2f}%
----------------------------------------
DECISION:    {decision}
REASON:      {reason}
========================================
"""

print(report_text)

# Save the decision to the results folder
with open(output_file, "w") as f:
    f.write(report_text)
print(f"[SUCCESS] Decision saved to: results/filtering_log.txt")