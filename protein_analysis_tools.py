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


def get_max_sequence_length_from_FASTA_file(fasta_filename):
    """
    Calculates the length of the sequence with the maximum length from a multiline fasta file
    
    Args: 
        fasta_filename (str): Path to the protein FASTA file.
    """
    max_length = max(len(seq) for nombre, seq in FASTA_iterator(fasta_filename))
    return max_length


def get_min_sequence_length_from_FASTA_file (fasta_filename):
    """
     Calculates the length of the sequence with the minimum length from a multiline fasta file
    
    Args: 
        fasta_filename (str): Path to the protein FASTA file.
    """
    min_length = min(len(seq) for nombre, seq in FASTA_iterator(fasta_filename))
    return min_length

get_min_sequence_length_from_FASTA_file("example_fasta_file.fa")

def get_longest_sequences_from_FASTA_file(fasta_filename):
    """
    Returns a sorted list of tuples (identifier, sequence) corresponding to the sequence(s) 
    with maximum length. The list is sorted (case insensitive) by the identifier.
    
    Args: 
        fasta_filename (str): Path to the protein FASTA file.
    """
    lista_tuples = []
    lista_max_tuples = []
    for name, seq in FASTA_iterator(fasta_filename):
        tuple_actual = (name,seq)
        lista_tuples.append(tuple_actual)
    max_len = max(len(tup[1]) for tup in lista_tuples) 
    for tup in lista_tuples:
        if len(tup[1]) == max_len:
            lista_max_tuples.append(tup)
    lista_sorted = sorted(lista_max_tuples, key=lambda x: x[0].lower())
    return lista_sorted


def get_shortest_sequences_from_FASTA_file(fasta_filename):
    """
    Returns a sorted list of tuples (identifier, sequence) corresponding to the sequence(s) 
    with minimum length. The list is sorted (case insensitive) by the identifier.
    
    Args: 
        fasta_filename (str): Path to the protein FASTA file.
    """
    list_tuples = []
    list_min_tuples =[]
    for name, seq in FASTA_iterator(fasta_filename):
        tuple_actual = (name,seq)
        list_tuples.append(tuple_actual)
    min_len = min((len(tup[1])) for tup in list_tuples)
    for tup in list_tuples:
        if len(tup[1]) == min_len:
            list_min_tuples.append(tup)
    list_sorted = sorted(list_min_tuples, key=lambda x:x[0].lower())
    return list_sorted


def get_molecular_weights(fasta_filename):
    """
    Calculates the molecular weight (Da) for each protein in a FASTA file.
    
    It accounts for the water molecule loss during peptide bond formation.
    
    Args:
        fasta_filename (str): Path to the protein FASTA file.
    Returns:
        dict: A dictionary (key: identifier, value: molecular weight as float).
    """
    dic_seq = {}
    dic_weight = {"G":75.1, "A":89.1, "S":105.1,"P":115.1,"V":117.1,"T":119.1,"C":121.2,"I":131.2,"L":131.2,"N":132.1,"D":133.1,"Q":146.2,"K":146.2,"E":147.1,"M":149.2,"H":155.2,"F":165.2,"R":174.2,"Y":181.2,"W":204.2}
    for nombre, seq in FASTA_iterator(fasta_filename):
        dic_seq[nombre] = 0
        peso_actual = sum(dic_weight[aa] for aa in seq)
        peso_actual = peso_actual - (len(seq)-1)*18.02
        dic_seq[nombre] = round(peso_actual, 2)
    return dic_seq


def get_sequence_with_min_molecular_weight(fasta_filename):
    """
    Returns a tuple (identifier, sequence) of the protein with the lowest molecular weight in a fasta file.
    If there are two or more proteins having the minimum molecular weight, it just returns the first one.
   
    Args: 
        fasta_filename (str): Path to the protein FASTA file.
    """
    dic_seqs = {}
    for nombre, seq in FASTA_iterator(fasta_filename):
        dic_seqs[nombre] = seq
    dic_pesos = get_molecular_weights(fasta_filename)
    min_peso = min(dic_pesos.items(), key=lambda x:x[1])
    id_min_peso = min_peso[0]
    min_tuple = (id_min_peso, dic_seqs[id_min_peso])
    return min_tuple


def get_mean_molecular_weight(fasta_filename):
    """
    Given a protein fasta file, calculates the mean of the molecular weights of all the proteins in the file
    
    Args: 
        fasta_filename (str): Path to the protein FASTA file.
    """
    dic_prote_peso = get_molecular_weights(fasta_filename)
    peso_tot = sum(dic_prote_peso[protein] for protein in dic_prote_peso)
    nprotein = len(dic_prote_peso)
    return(round(peso_tot/nprotein,2))

