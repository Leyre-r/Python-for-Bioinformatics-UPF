import sys
import math

def calculate_pdb_chain_mean_minimum_distances(pdb_file):
        """
        Calculates the mean of the minimum distances between all residue pairs 
        within each chain of a PDB file.

        Args:
            pdb_file (file object): file with Legacy PDB Format
        
        Returns:
            A dictionary where keys are chain identifiers (str) and 
              values are the calculated mean minimum distances (float).

        "Note: This script uses a nested loop approach for distance calculation. 
        For large complexes, performance could be improved using spatial indexing 
        (KD-Trees) or NumPy broadcasting."
        """
        dic_final = {}
        dic_chain_res = {}
        for line in pdb_file:
            line = line.strip("\n")
            if line[0:4] != "ATOM" and line[0:4] != "TER": continue
            res_id = line[22:26].strip()
            chain = line[21]
            #coordinates of each atom     
            coords = [float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())]
            if chain not in dic_chain_res:
                dic_chain_res[chain]={}
            if res_id not in dic_chain_res[chain]:
                dic_chain_res[chain][res_id] = []
            dic_chain_res[chain][res_id].append(coords)

        #Comparisons between atoms
        for cadena in dic_chain_res:
            res = dic_chain_res[cadena]
            list_res = res.keys()  
            lista_dist = []       
            for res1 in list_res:
                list_coords1 = dic_chain_res[cadena][res1]
                for res2 in list_res:
                    list_coords2 = dic_chain_res[cadena][res2]
                    min_dis = float('inf')
                    if res1 == res2: continue
                    for coord1 in list_coords1:
                        for coord2 in list_coords2:
                            restax = coord1[0] - coord2[0]
                            restay = coord1[1] - coord2[1]
                            restaz = coord1[2] - coord2[2]
                            #module of resulting vector
                            modulo = math.sqrt(restax**2+restay**2+restaz**2)
                            if modulo < min_dis: 
                                min_dis = modulo
                    lista_dist.append(min_dis)
            media = sum(lista_dist)/len(lista_dist)                                     
            dic_final[cadena] = round(media, 4)
        return dic_final

if __name__ == "__main__":
    if len(sys.argv) > 1:
        pdb_file = open(sys.argv[1], "r")
    else:
        pdb_file = sys.stdin

    resultado = calculate_pdb_chain_mean_minimum_distances(pdb_file)
    print(resultado)

    if pdb_file is not sys.stdin:
        pdb_file.close()

