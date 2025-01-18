from Bio import SeqIO
import time
import sys

reference_file = "project1a_reference_genome.fasta"
reads_file = "project1a_with_error_paired_reads.fasta"

print("Reference Genome:")
for record in SeqIO.parse(reference_file, "fasta"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}\n")

# Reading the paired reads
print("Paired Reads:")
for record in SeqIO.parse(reads_file, "fasta"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}\n")