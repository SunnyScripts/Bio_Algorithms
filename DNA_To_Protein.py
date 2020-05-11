# Created by Ryan Berg 5/8/20
# Synopsis: Convert DNA sequence to Protein Sequences

dnaSequence = "ATGCTAAGGTACCAGCATTAA"

Translation_Table = {
    "UUU": "Phenylalanine", "UCU": "Serine",    "UAU": "Tyrosine",      "UGU": "Cysteine",
    "UUC": "Phenylalanine", "UCC": "Serine",    "UAC": "Tyrosine",      "UGC": "Cysteine",
    "UUA": "Leucine",       "UCA": "Serine",    "UAA": "Stop",          "UGA": "Stop",
    "UUG": "Leucine",       "UCG": "Serine",    "UAG": "Stop",          "UGG": "Tryptophan",

    "CUU": "Leucine",       "CCU": "Proline",   "CAU": "Histidine",     "CGU": "Arginine",
    "CUC": "Leucine",       "CCC": "Proline",   "CAC": "Histidine",     "CGC": "Arginine",
    "CUA": "Leucine",       "CCA": "Proline",   "CAA": "Glutamine",     "CGA": "Arginine",
    "CUG": "Leucine",       "CCG": "Proline",   "CAG": "Glutamine",     "CGG": "Arginine",

    "AUU": "Isoleucine",    "ACU": "Threonine", "AAU": "Asparagine",    "AGU": "Serine",
    "AUC": "Isoleucine",    "ACC": "Threonine", "AAC": "Asparagine",    "AGC": "Serine",
    "AUA": "Isoleucine",    "ACA": "Threonine", "AAA": "Lysine",        "AGA": "Arginine",
    "AUG": "Methionine",    "ACG": "Threonine", "AAG": "Lysine",        "AGG": "Arginine",

    "GUU": "Valine",        "GCU": "Alanine",   "GAU": "Aspartic Acid", "GGU": "Glycine",
    "GUC": "Valine",        "GCC": "Alanine",   "GAC": "Aspartic Acid", "GGC": "Glycine",
    "GUA": "Valine",        "GCA": "Alanine",   "GAA": "Glutamic Acid", "GGA": "Glycine",
    "GUG": "Valine",        "GCG": "Alanine",   "GAG": "Glutamic Acid", "GGG": "Glycine"}

RNA_Compliment_Table = {
    "A": "U", "C": "G", "G": "C", "U": "A"}

def main():
    rnaSequence = ""
    rnaReverseCompliment = ""

    for nucleotide in dnaSequence:
        if nucleotide == "T":
            nucleotide = "U"
        rnaSequence += nucleotide
        rnaReverseCompliment = RNA_Compliment_Table[nucleotide] + rnaReverseCompliment

    proteinList = []
    for readingFrame in range(3):
        proteinList.extend(peptideChainToProteinList( rnaToPeptideChain( rnaSequence, readingFrame)))
        proteinList.extend(peptideChainToProteinList( rnaToPeptideChain( rnaReverseCompliment, readingFrame)))

    proteinList = sorted([protein for protein in proteinList if len(protein) > 0], key=len, reverse=True)
    print proteinList

def rnaToPeptideChain(rnaSequence, readingFrame):
    return [Translation_Table[rnaSequence[codonIndex:codonIndex+3]] for codonIndex in range(readingFrame, len(rnaSequence)-2, 3)]

def peptideChainToProteinList(peptideChain):
    putativeProteins = []
    completeProteins = []

    for residue in peptideChain:
        if residue == "Methionine":
            putativeProteins.append([])
        if len(putativeProteins) > 0:
            for eachList in putativeProteins:
                eachList.append(residue)
                if residue == "Stop":
                    completeProteins.append(putativeProteins.pop())
    return completeProteins

main()