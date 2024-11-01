import os
import json
from Bio import SeqIO

database_dir = "C://Users//Trust_pc_dz//Documents//BIOINFORMATIQUE//Database"
output_json = "Sars-cov-2.json"

data = {}

# Liste de tous les sous-dossiers dans "database"
try:
    variants = [d for d in os.listdir(database_dir) if os.path.isdir(os.path.join(database_dir, d))]
except FileNotFoundError:
    print(f"Erreur : Le dossier '{database_dir}' n'existe pas.")
    variants = []

# Parcourir chaque variant dans la liste
for variant_name in variants:
    print(f"Traitement du variant : {variant_name}")
    variant_path = os.path.join(database_dir, variant_name)
    data[variant_name] = {}  # Initialiser le dictionnaire pour chaque variant
    
    # Parcourir les sous-dossiers (années) dans chaque variant
    try:
        years = [y for y in os.listdir(variant_path) if os.path.isdir(os.path.join(variant_path, y))]
    except FileNotFoundError:
        print(f"Erreur : Le dossier '{variant_path}' n'existe pas.")
        continue
    
    for year in years:
        print(f"  Traitement de l'année : {year}")
        year_path = os.path.join(variant_path, year)
        sequences = {}  # Stocker les séquences pour chaque année
        
        try:
            # Parcourir chaque fichier dans le dossier de l'année
            for file_name in os.listdir(year_path):
                file_path = os.path.join(year_path, file_name)
                
                # Vérifie si le fichier est un fichier FASTA
                if file_path.endswith(".fasta"):
                    print(f"    Traitement du fichier : {file_name}")
                    for record in SeqIO.parse(file_path, "fasta"):
                        # Vérifier la longueur de la séquence
                        if len(record.seq) <= 4000:
                            # Stocke la séquence si elle est de longueur <= 4000
                            sequences[record.id] = {
                                "longueur": len(record.seq),
                                "séquence": str(record.seq)
                            }
                        else:
                            print(f"    Séquence ignorée (longueur > 4000) : {record.id}")
        
        except FileNotFoundError:
            print(f"Erreur : Le dossier '{year_path}' ou le fichier '{file_name}' n'existe pas.")
        except PermissionError:
            print(f"Erreur : Accès refusé pour le fichier '{file_path}'.")
        except Exception as e:
            print(f"Une erreur inattendue est survenue avec '{file_name}': {e}")
        
        # Ajouter les séquences de l'année au variant
        data[variant_name][year] = sequences

# Enregistrement dans un fichier JSON
try:
    with open(output_json, "w") as json_file:
        json.dump(data, json_file, indent=4)
    print(f"Les séquences FASTA et leurs attributs ont été enregistrés dans {output_json}")
except IOError:
    print("Erreur : Impossible d'écrire dans le fichier JSON.")
