import time
from collections import defaultdict
# from collections import Counter

def read_ref(file_path):
    """
    read reference genome
    """
    with open(file_path, 'r') as file:
        # skip first line
        next(file)
        # remove new line characters
        return file.read().replace('\n', '')

def read_reads(file_path):
    """
    read single reads file
    """
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if not line.startswith('>') and len(line.strip()) == 50]

    
def ham_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences."""
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def build_kmer_index(ref, k): 
    """
    building kmer index for reference genome
    """

    kmer_index = defaultdict(list) 
    for i in range(len(ref) - k + 1):
        kmer = ref[i:i + k]
        kmer_index[kmer].append(i)
    return kmer_index

def split_parts(seq, num_parts):
    """
    split given sequence into parts
    """
    part_length = len(seq) // num_parts
    return [seq[i * part_length: (i + 1) * part_length] for i in range(num_parts)]

def find_best_match_kmer(ref, read, k, kmer_index, max_mismatches):
    """
    find best matching kmer
    """
    possible_pos = set()
    # split into parts of L/3
    parts = split_parts(read, 3)

    # Collect candidate positions from k-mer matches
    for part in parts:
        # print(f"Read Part: {part}")
        for i in range(len(part) - k + 1):
            kmer = part[i:i + k]
            # print(f"Extracted K-mer: {kmer}")
            if kmer in kmer_index:
                possible_pos.update(kmer_index[kmer])

    # print(possible_pos)
    best_match = None
    min_mismatches = float('inf')

    # Evaluate candidate positions
    for pos in possible_pos:
        if pos + len(read) > len(ref):
            continue

        mismatches = ham_distance(ref[pos:pos + len(read)], read)

        if mismatches < min_mismatches and mismatches <= max_mismatches:
            min_mismatches = mismatches
            best_match = pos

    return best_match

def consensus_with_dict(ref, position_bases):
    """
    voting algorithm to determine substitution
    """
    res = set()
    for pos, bases in position_bases.items():
        # Initialize a dictionary to count occurrences of each base
        base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

        # Count occurrences of each base
        for base in bases:
            if base in base_counts:
                base_counts[base] += 1

        # Determine the consensus base 
        consensus_base = max(base_counts, key=base_counts.get)

        # Get the reference base
        ref_base = ref[pos]

        # If the consensus differs from the reference, it's a substitution
        if consensus_base != ref_base:
            res.add(f">S{pos} {ref_base} {consensus_base}")

    return res
    
'''
Result format: 
S[INDEX] [REFERENCE] [DONOR]
I[INDEX] [REFERENCE]
D[INDEX] [REFERENCE]
'''

def main():
    ref_path = 'project1a_reference_genome.fasta'
    reads_path = 'project1a_with_error_paired_reads.fasta'

    ref = read_ref(ref_path)
    reads = read_reads(reads_path)

    start_time = time.time()
    
    k = 16  # Fixed k-mer length
    max_mismatches = 3  # Maximum allowed mismatches
    
    kmer_index = build_kmer_index(ref, k)
    
    position_bases = defaultdict(list)

    for read in reads:
        best_match = find_best_match_kmer(ref, read, k, kmer_index, max_mismatches)
        if best_match is not None:
            for j in range(len(read)):
                position_bases[best_match + j].append(read[j])

    res = consensus_with_dict(ref, position_bases)

    for result in res:
        print(result)

    end_time = time.time()
    print(f"Execution time: {end_time - start_time} seconds")

    return res

    
if __name__ == "__main__":
    main()