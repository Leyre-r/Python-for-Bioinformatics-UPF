def FASTA_iterator (fasta_filename):
    """
    Efficiently parses a multiline FASTA file using a generator.

    Args:
        fasta_filename (str): Path to the FASTA file.

    Yields:
        tuple: A tuple containing the (identifier, sequence) for each entry.
    """
    sequence = []
    identifier = None
    with open(fasta_filename, "r") as fastafile:
        for line in fastafile:
            line = line.strip("\n")
            if line[0] == ">":
                if sequence:
                    full_seq = "".join(sequence)
                    yield (identifier, full_seq)
                identifier = line[1:]
                sequence = []
                
            else:
                sequence.append(line)
        if sequence:
            full_seq = "".join(sequence)
            yield (identifier, full_seq)


if __name__ == "__main__":
    # Example of use
    for identifier, sequence in FASTA_iterator("example.fa"):
        print(f"ID: {identifier} | Length: {len(sequence)}")
        


