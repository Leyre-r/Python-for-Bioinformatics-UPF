def calculate_aminoacid_frequencies(fasta_filename,subsequences_filename,number_of_repetitions, output_filename):
    """
    Calculates the absolute and relative frequencies of specific sub-sequences within a set of proteins from a FASTA file. 
    Only proteins containing the sub-sequence at least N times are considered. 

    Args:
        fasta_filename (str): Path to the protein FASTA file.
        subsequences_filename (str): Path to the subsequences file.
        number_of_repetitions (int): Number of repetitions of the subsequences
        output_filename (str): Name or path of the file where results will be saved.        
    
    Results:
        saved to an output file, sorted by descending frequency.
    """
    list_subseq = []
    dic_prote = {}
    dic_nsubseq = {}
    with open(subsequences_filename, "r") as subseq_file:
        for line in subseq_file:
            line = line.strip("\n")
            list_subseq.append(line)

    with open(fasta_filename, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip("\n")
            if line.startswith(">"):
                header = line[1:]
                dic_prote[header]= ""
            else:
                dic_prote[header]+= line
        
    for subseq in list_subseq:
        dic_nsubseq[subseq] = 0
        for key in dic_prote:
            proteina = dic_prote[key]
            nveces_subseq = proteina.count(subseq)
            if nveces_subseq >= number_of_repetitions:
                dic_nsubseq[subseq] += 1
    dic_nsubseq_sorted = dict(sorted(dic_nsubseq.items(), key=lambda item:item[1], reverse=True))
    
    with open(output_filename, "w") as output_file:
        nprote_keys = len(dic_prote)
        nsubseq_tot = len(list_subseq)
        line1 = f"#Number of proteins: {nprote_keys:>19} \n"
        line2 = f"#Number of subsequences: {nsubseq_tot:>15} \n"
        line3 = f"#Subsequence proportions: \n"
        output_file.write(line1 + line2 + line3)
        for key in dic_nsubseq_sorted:
            nsubseq = dic_nsubseq_sorted[key]
            proportion = nsubseq/nprote_keys
            line4 = f"{key:<10}{nsubseq:>10}{proportion:>20.4f}\n"
            output_file.write(line4)


if __name__ == "__main__":
    # Example of use
    calculate_aminoacid_frequencies(
        fasta_filename="example_fasta_file.fa",
        subsequences_filename="sequence_fragments.txt",
        number_of_repetitions=10,
        output_filename="results.txt"
    )
