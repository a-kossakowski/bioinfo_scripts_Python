import unittest
from OOP_BioPyth_prot_analysis import Protein

class TestProteinAnalysis(unittest.TestCase):

    def setUp(self):
        # Valid protein instance for testing
        self.valid_protein = Protein("MLLLIIAGYY")

    def test_valid_protein_sequence(self):
        # Test with a valid sequence
        self.assertIsInstance(self.valid_protein.molecular_weight(), float)
        self.assertIsInstance(self.valid_protein.isoelectric_point(), float)
        self.assertIsInstance(self.valid_protein.amino_acid_composition(), dict)
        self.assertIsInstance(self.valid_protein.count_amino_acid("A"), int)
        self.assertIsInstance(self.valid_protein.get_domains("LII"), list)

    def test_invalid_protein_sequence(self):
        # Test initialization with an invalid sequence
        with self.assertRaises(ValueError):
            Protein("MLL#IIAGYY")

    def test_valid_amino_acid(self):
        # Test counting a valid amino acid
        self.assertEqual(self.valid_protein.count_amino_acid("A"), 2)

    def test_invalid_amino_acid(self):
        # Test counting an invalid amino acid
        with self.assertRaises(ValueError):
            self.valid_protein.count_amino_acid("O")

    def test_domain_search_valid(self):
        # Test domain search with a valid domain
        self.assertEqual(self.valid_protein.get_domains("LII"), [2])

    def test_domain_search_invalid(self):
        # Test domain search with an invalid domain
        self.assertEqual(self.valid_protein.get_domains("XYZ"), [])

if __name__ == '__main__':
    unittest.main()

