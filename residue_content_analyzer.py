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


def get_proteins_ratio_by_residue_threshold(filename, residue, relative_threshold, 
                                            absolute_threshold):
    """  
    Calculates the proportion of proteins from the fasta file that meet both relative and absolute 
    frequency thresholds for a specific residue.

    Args:
        filename (str): Path to the protein FASTA file.
        residue (str): The amino acid residue to analyze.
        relative_threshold (float): Minimum relative frequency required (0.0 to 1.0).
        absolute_threshold (int): Minimum absolute occurrences required.

    Returns:
        float: The ratio (0.0 to 1.0) of proteins that passed the filters.
    """
    proteins_total = 0
    proteins_pass = 0
    for nombre, secuencia in FASTA_iterator(filename):
        proteins_total += 1
        nres = secuencia.count(residue)
        nres_rel = nres/len(secuencia)
        if nres >= absolute_threshold and nres_rel >= relative_threshold:
            proteins_pass += 1
    ratio = proteins_pass/proteins_total
    return ratio


def aminoacid_frequency(filename, N, M,output_file_name):
  """
   Generates a report file with protein fragments and amino acid composition.

   The output file includes the protein identifier, the first N residues, the last M 
    residues, and the absolute frequency of each amino acid found.

   Args:
    filename (str): Path to the protein FASTA file.
    N (int): Number of residues to extract from the N-terminus.
    M (int): Number of residues to extract from the C-terminus.
    output_file_name (str) : Name or path of the file where results will be saved.
  """
  with open(output_file_name, "w") as out:
    for nombre, secuencia in FASTA_iterator(filename):
        dic_freq = {}
        nfirst = secuencia[:N]
        mlast = secuencia[-M:]
        out.write(f'{nombre}\t{nfirst}\t{mlast}\t')
        for aa in secuencia:
            if aa not in dic_freq:
                dic_freq[aa]= 1
            else: 
                dic_freq[aa] += 1
        lista1 = []
        for aa in dic_freq:
            naa = dic_freq[aa]
            lista1.append(f'{aa}:{naa}')
        line1 = ",".join(lista1)
        out.write(line1 + "\n")
        