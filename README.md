# PASM- phylogenetic analysis-MSA statistics-mutation identification

PASM is a Python-based program that builds phylogenetic trees, computes multiple sequence alignment (MSA) statistics, and identifies key mutations. It performs pairwise global alignment using the Needleman–Wunsch algorithm, constructs a distance matrix based on the Kimura 2-parameter (K80) model, and generates phylogenetic trees using the UPGMA algorithm.

The program also calculates conservation scores and detects key mutations from multiple sequence alignments.

PASM is designed as a pipeline for biologists and bioinformaticians, particularly for applications in  evolution and sequence analysis.

## Python packages
```
install with pip
Biopython
numpy
matplotlib
```

# How to run the program:
Input: The user can input a FASTA file of sequences and use the file_to_folder() function to convert the file to a folder of FASTA files of one sequence or directly use a folder as input. 

In Line 13
path = "specify path to save the folder"
In Line 42
file_path = "enter the fasta file path"
In Line 43
file_to_folder(file_path,name of the folder)In line 
In line 130:
os.chdir("enter the path of where you want to save the pairwise alignment output files")
In line 201: 
dataset_path ="enter the path of the folder that contain the sequences"
In line 647:
for i in os.listdir("enter the path of the folder that contain the sequences"):
In line 656:
stats=Statistics("enter the path of the multiple sequence alignment file to calculate the alignment statistics")

# Project structure: 
The program is composed of 3 main classes 
Global alignment class (calculate the distance matrix)
2. UPGMA class (construct phylogenetic tree)
3. Statistics class (calculate all statistics, identify mutation, visualize the tree 

# Output: 
phylogenetic tree files
alignment statistics file 
key mutation file 
FASTA file for the pairwise alignment 


