import json
from collections import defaultdict
import random

class Handler:
    def __init__(self, json_path):
        with open(json_path, "r") as file:
            self.data = json.load(file)
    
    def count_sequences_per_year(self):
        ''' Dictionnaire pour stocker le nombre de séquences par année'''
        sequences_per_year = defaultdict(int)
        
        for variant, years in self.data.items():
            for year, sequences in years.items():
                sequences_per_year[year] += len(sequences)
        
        return dict(sequences_per_year)
    
    def calculate_sequence_lengths(self):
        ''' Dictionnaire pour stocker les longueurs de séquences par variant et par année'''
        sequence_lengths = {}
        
        for variant, years in self.data.items():
            sequence_lengths[variant] = {}
            for year, sequences in years.items():
                sequence_lengths[variant][year] = {seq_id: len(details["séquence"]) for seq_id, details in sequences.items()}
        
        return sequence_lengths
    
    def calculate_base_percentages(self):
        ''' Dictionnaire pour stocker les pourcentages de A, T, C, G par variant et par année'''
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
        ''' Afficher toutes les statistiques'''
        print("Nombre de séquences par année:", self.count_sequences_per_year())
        print("\nLongueur des séquences par variant et année:", self.calculate_sequence_lengths())
        print("\nPourcentages de bases par variant, année et séquence:", self.calculate_base_percentages())

    def choisir_sequences_random(self):
        toutes_les_sequences = []
        for variant, annees in self.data.items():
            for annee, sequences in annees.items():
                for seq_id, details in sequences.items():
                    toutes_les_sequences.append({
                        "id": seq_id,
                        "variant": variant,
                        "annee": annee,
                        "longueur": details["longueur"],
                        "séquence": details["séquence"]
                    })

        sequences_uniques = list({seq["séquence"]: seq for seq in toutes_les_sequences}.values())

        if len(sequences_uniques) < 3:
            print("Pas assez de séquences uniques pour en choisir trois.")
            return []

        sequences_choisies = random.sample(sequences_uniques, 3)
        tableau_sequences = []
        
        for seq in sequences_choisies:

            sequence_extraite = seq["séquence"][random.randint(20, 35): random.randint(20, 35) + random.randint(10, 15)]
            
   
            while len(sequence_extraite) < 5 or sequence_extraite == '':
                print("Séquence vide ou trop courte trouvée, génération d'une nouvelle séquence.")
                sequence_extraite = seq["séquence"][random.randint(20, 25): random.randint(20, 25) + random.randint(10, 15)]
            
            tableau_sequences.append(sequence_extraite)

        print("Séquences choisies :")
        for i, seq in enumerate(sequences_choisies, start=1):
            print(f"  Séquence {i}:")
            print(f"    ID : {seq['id']}")
            print(f"    Variant : {seq['variant']}")
            print(f"    Année : {seq['annee']}")
            print(f"    Longueur : {seq['longueur']}")
            print(f"    Séquence : {seq['séquence']}\n")

        return tableau_sequences
    
    def get_sequences_by_variant_and_year(self, variant, year):
        """
        Récupère toutes les séquences pour un variant donné et une année donnée.
        
        :param variant: Le nom du variant à rechercher.
        :param year: L'année à rechercher.
        :return: Un dictionnaire contenant les séquences ou un message si aucun résultat n'est trouvé.
        """
        if variant not in self.data:
            return f"Variant '{variant}' introuvable dans les données."
        
        if year not in self.data[variant]:
            return f"Aucune séquence trouvée pour le variant '{variant}' en {year}."
        
        # Récupérer toutes les séquences pour ce variant et cette année
        sequences = self.data[variant][year]
        return {seq_id: details["séquence"] for seq_id, details in sequences.items()}
    
    def export_to_fasta(self, variant, output_fasta, year=None):
        """
        Exporte les séquences d'un variant donné (et éventuellement une année spécifique)
        dans un fichier FASTA.
        """
        if variant not in self.data:
            print(f"Variant '{variant}' introuvable dans les données.")
            return

        sequences = []
        for year_key, seq_data in self.data[variant].items():
            if year and year_key != str(year):
                continue  # Ignore les années non sélectionnées
            for seq_id, details in seq_data.items():
                sequences.append((seq_id, details["séquence"]))

        if not sequences:
            print("Aucune séquence trouvée pour ce variant et cette année.")
            return

        # Écrire les séquences dans un fichier FASTA
        with open(output_fasta, "w") as fasta_file:
            for seq_id, sequence in sequences:
                print(f"Écriture de {seq_id} dans le fichier")
                fasta_file.write(f">{seq_id}\n{sequence}\n")

        print(f"Fichier FASTA '{output_fasta}' généré avec succès.")

