from Bio import SeqIO

def find_gaps(genome_file):
    gap_coordinates = []
    with open(genome_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = str(record.seq)
            for i in range(len(sequence)):
                if sequence[i] == 'N':
                    start = i
                    end = i
                    while end + 1 < len(sequence) and sequence[end + 1] == 'N':
                        end += 1
                    gap_coordinates.append((start + 1, end + 1, record.id))

    return gap_coordinates

if __name__ == "__main__":
    genome_file = "your_genome.fasta"  # Replace with the path to your genome file
    gap_coordinates = find_gaps(genome_file)

    if gap_coordinates:
        print("Gap Coordinates:")
        for start, end, record_id in gap_coordinates:
            print(f"Sequence ID: {record_id}, Start: {start}, End: {end}")
    else:
        print("No gaps found in the genome.")

