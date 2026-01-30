def FASTA_iterator(fasta_filename):
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


def compare_fasta_file_identifiers(fasta_filenames_list):
    """
    Compares identifiers across multiple FASTA files to identify commonalities and differences.

    Args:
        fasta_filenames_list (list): A list of paths to FASTA files.
    Returns:
        dict: A dictionary containing:
            - "intersection" (set): only Identifiers present in ALL input files.
            - "union" (set): all Identifiers (unique) found across all files.
            - "frequency" (dict): Identifier count (key: ID, value: number of files where it appears).
            - "specific" (dict): Exclusive identifiers (key: filename, value: set of IDs unique 
               to that file, i.e. identifiers that are exclusive in that fasta file).         
     """
    dic_results = {}
    dic_identifiers = {}
    dic_specific = {}
    dic_files = {}
    identifier_total_set = set()
    identifier_common_set = None

    for file in fasta_filenames_list:
        set_actual = set()
        for nombre, secuencia in FASTA_iterator(file):
            nombre = nombre.lower()
            set_actual.add(nombre)
           
            if nombre not in dic_identifiers:
                dic_identifiers[nombre] = 1
            else:
                dic_identifiers[nombre] += 1

        dic_files[file] = set_actual.copy()
        #union
        identifier_total_set.update(set_actual)
        #intersection
        if identifier_common_set == None:
            identifier_common_set = set_actual.copy()
        else: 
            identifier_common_set.intersection_update(set_actual)
    
    for file in dic_files:
        dic_specific[file] = set()
        identifiers = dic_files[file]
        for name in identifiers:
            if dic_identifiers[name] == 1:
                dic_specific[file].add(name)

    dic_results["intersection"] = identifier_common_set
    dic_results["union"] = identifier_total_set
    dic_results["frequency"] = dic_identifiers
    dic_results["specific"] = dic_specific

    return dic_results



