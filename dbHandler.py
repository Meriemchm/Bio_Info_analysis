import json
from collections import defaultdict

class Handler:
    def __init__(self, json_path):
        with open(json_path, "r") as file:
            self.data = json.load(file)
    
    def count_sequences_per_year(self):
        # Dictionnaire pour stocker le nombre de séquences par année
        sequences_per_year = defaultdict(int)
        
        for variant, years in self.data.items():
            for year, sequences in years.items():
                sequences_per_year[year] += len(sequences)
        
        return dict(sequences_per_year)
    
    def calculate_sequence_lengths(self):
        # Dictionnaire pour stocker les longueurs de séquences par variant et par année
        sequence_lengths = {}
        
        for variant, years in self.data.items():
            sequence_lengths[variant] = {}
            for year, sequences in years.items():
                sequence_lengths[variant][year] = {seq_id: len(details["séquence"]) for seq_id, details in sequences.items()}
        
        return sequence_lengths
    
    def calculate_base_percentages(self):
        # Dictionnaire pour stocker les pourcentages de A, T, C, G par variant et par année
        base_percentages = {}
        
        for variant, years in self.data.items():
            base_percentages[variant] = {}
            for year, sequences in years.items():
                base_percentages[variant][year] = {}
                for seq_id, details in sequences.items():
                    sequence = details["séquence"]
                    total_length = len(sequence)
                    counts = {
                        "A": sequence.count("A"),
                        "T": sequence.count("T"),
                        "C": sequence.count("C"),
                        "G": sequence.count("G")
                    }
                    # Calculer les pourcentages pour chaque base
                    base_percentages[variant][year][seq_id] = {
                        base: (count / total_length) * 100 for base, count in counts.items()
                    }
        
        return base_percentages
    
    def display_statistics(self):
        # Afficher toutes les statistiques
        print("Nombre de séquences par année:", self.count_sequences_per_year())
        print("\nLongueur des séquences par variant et année:", self.calculate_sequence_lengths())
        print("\nPourcentages de bases par variant, année et séquence:", self.calculate_base_percentages())


