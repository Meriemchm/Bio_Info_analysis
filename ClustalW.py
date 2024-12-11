from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import os

class ClustalHandler:
    def __init__(self, clustalw_path):
        self.clustalw_path = clustalw_path

    def align_sequences(self, input_fasta, output_aln="output.aln"):
        """Utilise ClustalW pour aligner des séquences."""
        # Appel à ClustalW
        clustalw_cline = ClustalwCommandline(self.clustalw_path, infile=input_fasta)
        stdout, stderr = clustalw_cline()

        # Vérification des erreurs
        if stderr:
            print("Erreur lors de l'exécution de ClustalW :", stderr)
            return None

        # Lecture des alignements
        alignment = AlignIO.read(output_aln, "clustal")
        print("Alignement terminé. Résultat :")
        print(alignment)
        return alignment
