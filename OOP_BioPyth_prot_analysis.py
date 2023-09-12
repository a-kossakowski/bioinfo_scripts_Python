from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Define a class to encapsulate protein sequence analysis
class Protein:
    def __init__(self, sequence):
        # Constructor initializes the instance with a protein sequence
        self.sequence = sequence
        # Create an analysis object for the protein sequence
        self.analysis = ProteinAnalysis(sequence)

    def molecular_weight(self):
        # Method to calculate and return the molecular weight of the protein
        # Units: Daltons (Da)
        return self.analysis.molecular_weight()

    def get_domains(self, domain_seq):
        # Method to find and return the starting indices of domains within the sequence
        # Parameter 'domain_seq' is the sequence of the domain to search for
        # Units: Indices representing positions in the sequence
        return [i for i in range(len(self.sequence)) if self.sequence.startswith(domain_seq, i)]

    def isoelectric_point(self):
        # Method to calculate and return the isoelectric point (pI) of the protein
        # Units: pH
        return self.analysis.isoelectric_point()

    def amino_acid_composition(self):
        # Method to calculate and return the amino acid composition as a percentage
        return self.analysis.get_amino_acids_percent()

    def count_amino_acid(self, amino_acid):
        # Method to count the occurrences of a specific amino acid in the protein sequence
        # Parameter 'amino_acid' is the single-letter code of the amino acid
        return self.analysis.count_amino_acids()[amino_acid]

    def charge_at_pH(self, pH):
        # Method to calculate and return the charge of the protein at a given pH
        # Parameter 'pH' is the pH value at which to calculate the charge
        # Units: Charge, a value indicating the net electrical charge of the protein at the specified pH

# Create an instance of the 'Protein' class with a protein sequence
protein1 = Protein("MLLLIIAGYY")

# Perform various analyses and display the results
print("Protein Sequence:", protein1.sequence)
print("Molecular Weight:", protein1.molecular_weight(), "Da")
print("Isoelectric Point:", protein1.isoelectric_point())
print("Amino Acid Composition:", protein1.amino_acid_composition())
print("Count of 'A' Amino Acid:", protein1.count_amino_acid("A"))
print("Charge at pH 7.0:", protein1.charge_at_pH(7.0))
print("Domain's Index:", protein1.get_domains("LII"))
