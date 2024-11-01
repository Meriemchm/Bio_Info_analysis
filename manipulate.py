import os
import json
from Bio import SeqIO, Entrez

database_dir = "C://Users//Trust_pc_dz//Documents//BIOINFORMATIQUE//Database"
output_json = "Sars-cov-2.json"

data = {}

# Liste de tous les sous-dossiers dans "database"
try:
    subdirs = [d for d in os.listdir(database_dir) if os.path.isdir(os.path.join(database_dir, d))]
except FileNotFoundError:
    print(f"Erreur : Le dossier '{database_dir}' n'existe pas.")
    subdirs = []

# Parcourir chaque sous-dossier dans la liste
for dir_name in subdirs:
    print(f"Traitement du dossier : {dir_name}")
    dir_path = os.path.join(database_dir, dir_name)
    
    # Dictionnaire pour stocker les séquences et leurs attributs dans le sous-dossier
    sequences = {}
    
    try:
        # Parcourir chaque fichier dans le dossier
        for file_name in os.listdir(dir_path):
            print(f"Traitement du fichier : {file_name}")
            file_path = os.path.join(dir_path, file_name)
            
            # Vérifie si le fichier est un fichier FASTA
            if file_path.endswith(".fasta"):
                for record in SeqIO.parse(file_path, "fasta"):
                    # Stocke chaque séquence avec ses attributs
                    sequences[record.id] = {
                        "longueur": len(record.seq),
                        "séquence": str(record.seq)
                    }
    
    except FileNotFoundError:
        print(f"Erreur : Le dossier '{dir_path}' ou le fichier '{file_name}' n'existe pas.")
    except PermissionError:
        print(f"Erreur : Accès refusé pour le fichier '{file_path}'.")
    except Exception as e:
        print(f"Une erreur inattendue est survenue avec '{file_name}': {e}")
    
    # Ajouter les séquences au dictionnaire principal sous le nom du dossier
    data[dir_name] = sequences

# Enregistrement dans un fichier JSON
try:
    with open(output_json, "w") as json_file:
        json.dump(data, json_file, indent=4)
    print(f"Les séquences FASTA et leurs attributs ont été enregistrés dans {output_json}")
except IOError:
    print("Erreur : Impossible d'écrire dans le fichier JSON.")