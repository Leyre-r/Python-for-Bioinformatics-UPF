# Bioinformatics Python Toolbox
This repository contains a collection of Python scripts developed during my Introduction to Python course (MSc Bioinformatics for Health Sciences, UPF/UB). Currently organizing the content by topics.

## 🛠️ Tools Included

| Script | Description | 
| :--- | :--- | 
| `fasta_iterator.py` | Efficient multiline FASTA parser using generators. | 
| `protein_analysis_tools.py` | Calculates molecular weights and sequence metrics. | 
| `pdb_residue_distances.py` | Computes mean minimum distances between residues in PDB chains. |
| `fasta_comparator.py` | Set-theory based comparison between multiple FASTA files (Union, Intersection). | 
| `residue_content_analyzer.py` | Filters proteins by residue thresholds and generates composition reports. |
| `calculate_subseq_frequencies.py`| Calculates the absolute and relative frequencies of specific sub-sequences within a set of proteins. |

## How to use
The scripts can be run via command line and some of them can be imported as modules:
```bash
python3 pdb_residue_distances.py input_structure.pdb
