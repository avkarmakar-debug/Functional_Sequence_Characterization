import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord

# CONFIGURATION
INPUT_FILE = "rat.fasta"
OUTPUT_PROTEIN = "protein.fasta"
BLAST_XML = "blast_result.xml"
REPORT_FILE = "functional_annotation.txt"


# STEP 1: Read DNA sequence

print(f"Reading {INPUT_FILE}...")

try:
    dna_record = SeqIO.read(INPUT_FILE, "fasta")
except ValueError:
    print(" Error: Input file contains multiple sequences. Use SeqIO.parse().")
    sys.exit()
except FileNotFoundError:
    print(f" Error: File '{INPUT_FILE}' not found.")
    sys.exit()

dna_seq = dna_record.seq
print(f"Sequence ID: {dna_record.id}")
print(f" DNA Length: {len(dna_seq)} bp")


# STEP 2: Quality check (GC%)

gc = gc_fraction(dna_seq) * 100
print(f" GC Content: {gc:.2f}%")


# STEP 3: DNA validation

if len(dna_seq) < 100:  # Lowered limit for testing
    print(" Sequence too short. Analysis stopped.")
    sys.exit()


# STEP 4: Translate DNA -> Protein (ORF Finder)

print("\nTranslating DNA to protein (searching for longest ORF)...")

# Logic: Check all 3 reading frames, split by stop codons (*), find longest string starting with M
frames = [0, 1, 2]
candidates = []

for frame in frames:
    # Translate full length of this frame
    # (pad sequence with N if needed to make length multiple of 3)
    seq_len = len(dna_seq)
    remainder = (seq_len - frame) % 3
    if remainder != 0:
        # Just slice off the end for translation purposes
        trimmed_seq = dna_seq[frame : seq_len - remainder]
    else:
        trimmed_seq = dna_seq[frame:]
        
    full_trans = trimmed_seq.translate(table=1)
    
    # Split translation at stop codons
    parts = full_trans.split("*")
    
    # Check each part for a Start Codon (M)
    for part in parts:
        if "M" in part:
            start_index = part.find("M")
            candidates.append(part[start_index:])

if candidates:
    protein_seq = max(candidates, key=len)
    print(f" ORF Found! Protein Length: {len(protein_seq)} aa")
else:
    print("⚠️ No valid ORF found. Translating raw sequence.")
    protein_seq = dna_seq.translate(to_stop=True)
    print(f"Protein Length: {len(protein_seq)} aa")

# Save protein FASTA
protein_record = SeqRecord(
    protein_seq,
    id=f"{dna_record.id}_protein",
    description="Translated_ORF"
)
SeqIO.write(protein_record, OUTPUT_PROTEIN, "fasta")
print(f" Saved protein sequence to: {OUTPUT_PROTEIN}")


# STEP 5: Homology search (ONLINE BLASTp)

print(f"\n Submitting {len(protein_seq)} aa to NCBI BLASTp...")
print("   (Please wait... this commonly takes 1-5 minutes)")

try:
    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="nr",
        sequence=protein_seq
    )
    
    with open(BLAST_XML, "w") as out_handle:
        out_handle.write(result_handle.read())
    
    print(f" BLAST search completed. Results saved to {BLAST_XML}")

except Exception as e:
    print(f" BLAST Connection Error: {e}")
    print("Check your internet connection or try again later.")
    sys.exit()


# STEP 6: Functional Annotation Report

print("\nParsing BLAST results and generating report...")

try:
    with open(BLAST_XML) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    with open(REPORT_FILE, "w") as out:
        out.write("FUNCTIONAL ANNOTATION REPORT\n")
        out.write("="*50 + "\n\n")
        out.write(f"Query ID: {blast_record.query}\n")
        out.write(f"Query Length: {blast_record.query_length} aa\n")
        out.write(f"Database: {blast_record.database}\n\n")

        if not blast_record.alignments:
            out.write("No significant hits found.\n")
            print("⚠️ No hits found in BLAST database.")
        else:
            out.write(f"Found {len(blast_record.alignments)} hits. Showing Top 5:\n\n")
            
            count = 0
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    count += 1
                    out.write(f"Hit #{count}: {alignment.title[:80]}...\n") # Truncate long titles
                    out.write(f"Accession: {alignment.accession}\n")
                    out.write(f"E-value: {hsp.expect}\n")
                    out.write(f"Identity: {hsp.identities}/{hsp.align_length} ({hsp.identities/hsp.align_length:.1%})\n")
                    out.write("-" * 40 + "\n")
                    
                    if count >= 5: break 
                if count >= 5: break
                
    print(f" Analysis Complete! Open '{REPORT_FILE}' to see results.")

except Exception as e:
    print(f" Error parsing XML: {e}")


    