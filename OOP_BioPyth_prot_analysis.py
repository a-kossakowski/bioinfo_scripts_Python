from Bio.SeqUtils.ProtParam import ProteinAnalysis

class Protein:
    """
    A class utilizng BioPython's ProteinAnalysis.
    """
    
    def __init__(self, sequence):
        """
        Initialize the instance with a protein sequence.
        
        :param sequence: A string for the protein sequence.
        """
        if not sequence.isalpha():
            raise ValueError("Invalid protein sequence. Sequence should contain only alphabetic characters.")
        
        self.sequence = sequence
        self.analysis = ProteinAnalysis(sequence)

    def molecular_weight(self):
        """
        Calculate molecular weight of the protein.
        
        :return: Molecular weight in Daltons (Da).
        """
        return self.analysis.molecular_weight()

    def get_domains(self, domain_seq):
        """
        Find the starting indices of domains within the protein sequence.
        
        :param domain_seq: The sequence of the domain to search for.
        :return: List of indices representing domain positions in the sequence.
        """
        return [i for i in range(len(self.sequence)) if self.sequence.startswith(domain_seq, i)]

    def isoelectric_point(self):
        """
        Calculate the isoelectric point (pI) of the protein.
        
        :return: Isoelectric point (pI) as pH value.
        """
        return self.analysis.isoelectric_point()

    def amino_acid_composition(self):
        """
        Calculate the amino acid composition as a percentage.
        
        :return: Dictionary with amino acid composition percentage.
        """
        return self.analysis.get_amino_acids_percent()

    def count_amino_acid(self, amino_acid):
        """
        Count the occurrences of a specific amino acid in the protein sequence.
        
        :param amino_acid: Single-letter code of the amino acid.
        :return: Count of the amino acid.
        """
        if len(amino_acid) != 1 or not amino_acid.isalpha():
            raise ValueError("Invalid amino acid. Amino acid should be a single alphabetic character.")
        
        return self.analysis.count_amino_acids().get(amino_acid.upper(), 0)

# Example Usage
try:
    protein1 = Protein("MLLLIIAGYY")
    print("Protein Sequence:", protein1.sequence)
    print("Molecular Weight:", protein1.molecular_weight(), "Da")
    print("Isoelectric Point:", protein1.isoelectric_point())
    print("Amino Acid Composition:", protein1.amino_acid_composition())
    print("Count of 'A' Amino Acid:", protein1.count_amino_acid("A"))
    print("Domain's Index:", protein1.get_domains("LII"))
except ValueError as e:
    print(e)

