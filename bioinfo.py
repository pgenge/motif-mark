#!/usr/bin/env python
# Author: Palak Genge pgenge@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.8"        # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

iupac_notation_dict = {
    "A":"[Aa]",
    "C":"[Cc]",
    "G":"[Gg]",
    "T":"[Tt]",
    "U":"[UuTt]",
    "W":"[AaTtUu]",
    "S":"[CcGg]",
    "M":"[AaCc]",
    "K":"[GgTtUu]",
    "R":"[AaGg]",
    "Y":"[CcTtUu]",
    "B":"[CcGgTt]",
    "D":"[AaGgTtUu]",
    "H":"[AaCcTtUu]",
    "V":"[AaCcGg]",
    "N":"[AaCcGgTtUu]",
    "Z":"[]"
}

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """This function calculates the average quality score of the whole phred string"""
    sum = 0
    for x in phred_score:
        sum += convert_phred(x)
    return sum/len(phred_score)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA: str):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    DNA = DNA.upper()         #Make sure sequence is all uppercase
    Gs = DNA.count("G")       #count the number of Gs
    Cs = DNA.count("C")       #count the number of Cs
    return (Gs+Cs)/len(DNA)

# def oneline_fasta(fasta):
#     seq = ""
#     newfasta = open("onelinefasta.fa", "w")
#     with open (fasta, "r") as fa:
#         for line in fa:
#             line = line.strip("\n")
#             if line.startswith(">"):
#                 if seq != "":
#                     newfasta.write(f"{seq}\n")
#                     newfasta.write(f"{line}\n")
#                 else:
#                     seq += line
#                     if  seq:
#                         print(seq)
#                         newfasta.write(f"{seq}\n")
#                         newfasta.write(f"{line}\n")

def oneline_fasta(fasta):
    '''Makes FASTA sequences on one line. Writes out to the file
    fa_one_line.fa. Returns the number of records so they can be
    manually compared to the number of header lines in the output file,
    to confirm the output file is accurate.'''
    # make dict with headers as keys and sequences as values
    seq_dict = {}
    with open(fasta, 'r') as fa:
        line_count = 0
        for line in fa:
            line_count +=1
            line = line.strip('\n')
            # only get header lines
            if line[0] == '>':
                header_line = line
            # populate dict with seq lines (non-header lines)
            else:
                if header_line not in seq_dict:
                    seq_dict[header_line] = line
                else:
                    seq_dict[header_line] += line
    # write out to file
    onelinefasta = open('onelinefasta.fa', 'w')
    for keys,vals in seq_dict.items():
        onelinefasta.write(f"{str(keys)}\n{str(vals)}\n")
    onelinefasta.close()
    return onelinefasta


#unit tests
if __name__ == "__main__":
    assert gc_content("GCGCGC") == 1, "messed up calc when all GC"
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")
    print("correctly calculated GC content")
