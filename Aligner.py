from Bio import SeqIO, Entrez
import os
from itertools import combinations
from collections import defaultdict

class AlignmentHandler:
    def __init__(self, match=1, mismatch=-1, gap=-2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap


    def score(self, char1, char2):
        """Calcule le score entre deux caractères."""
        if char1 == char2:
            return self.match
        else:
            return self.mismatch
    
    
    def explore_alignments(self, i, j, current_align1, current_align2,alignments,score_matrix,seq1,seq2):
    # Condition d'arrêt
        if i < 2 and j < 2:
            alignments.append((current_align1, current_align2))
            return

        diag = score_matrix[i - 1][j - 1] if i > 0 and j > 0 else float('-inf')
        up = score_matrix[i - 1][j] if i > 0 else float('-inf')
        left = score_matrix[i][j - 1] if j > 0 else float('-inf')

        max_value = max(diag, up, left)

        if max_value == diag:
            self.explore_alignments(i - 1, j - 1, seq1[i - 2] + current_align1, seq2[j - 2] + current_align2,alignments,score_matrix,seq1,seq2)
        if max_value == up:
            self.explore_alignments(i - 1, j, seq1[i - 2] + current_align1, '-' + current_align2,alignments,score_matrix,seq1,seq2)
        if max_value == left:
            self.explore_alignments(i, j - 1, '-' + current_align1, seq2[j - 2] + current_align2,alignments,score_matrix,seq1,seq2)

        return alignments    

    def global_alignment(self, seq1, seq2):

        n = len(seq1) + 1
        m = len(seq2) + 1

        score_matrix = [[0] * m for _ in range(n)]
        #print(score_matrix)
        # traceback_matrix = [[None] * m for _ in range(n)]

        # Rempli gap
        for i in range(1, n):
            score_matrix[i][0] = score_matrix[i - 1][0] + self.gap
            # traceback_matrix[i][0] = 'up'
        for j in range(1, m):
            score_matrix[0][j] = score_matrix[0][j - 1] + self.gap
            # traceback_matrix[0][j] = 'left'


        for i in range(1, n):
            for j in range(1, m):
                if seq1[i - 1] == seq2[j - 1]:
                    score = self.match
                else:
                    score = self.mismatch

                #best score
                diag = score_matrix[i - 1][j - 1] + score
                up = score_matrix[i - 1][j] + self.gap
                left = score_matrix[i][j - 1] + self.gap
                max_score = max(diag, up, left)
                score_matrix[i][j] = max_score

                # if max_score == diag:
                #     traceback_matrix[i][j] = 'diag'
                # elif max_score == up:
                #     traceback_matrix[i][j] = 'up'
                # else:
                #     traceback_matrix[i][j] = 'left'

        alignments = []  
        i, j = n - 1, m - 1
        initial_align1 = seq1[i-1]  
        initial_align2 = seq2[j-1]  

        alignments = self.explore_alignments(i, j, initial_align1, initial_align2,alignments,score_matrix,seq1,seq2)
                        
        #print("Matrice de score :")
        #for row in score_matrix:
        #    print(row)    
        
        return alignments, score_matrix[n - 1][m - 1]

    def local_alignment(self, seq1, seq2):
        n = len(seq1) + 1
        m = len(seq2) + 1

        score_matrix = [[0] * m for _ in range(n)]

        # max_pos = (0, 0)
        
        for i in range(1, n):
            for j in range(1, m):
                if seq1[i - 1] == seq2[j - 1]:
                    score = self.match
                else:
                    score = self.mismatch

                diag = score_matrix[i - 1][j - 1] + score
                up = score_matrix[i - 1][j] + self.gap
                left = score_matrix[i][j - 1] + self.gap
                
                score_matrix[i][j] = max(0,diag, up, left)
                
                # # Track maximum score position
                # if score_matrix[i][j] > max_score:
                #     max_score = score_matrix[i][j]

        End = score_matrix[n - 1][m - 1]
        diag = score_matrix[n - 2][m - 2]
        up = score_matrix[n - 2][m-1]
        left = score_matrix[n-1][m - 2]
        score_final = max(diag,up,left,End) 

        if  score_final == End:
            i, j = n - 1, m - 1
            initial_align1 = seq1[i-1]  
            initial_align2 = seq2[j-1]
        elif   score_final == diag:
            i, j = n - 2, m - 1
            initial_align1 = seq1[i-2]  
            initial_align2 = seq2[j-2]        
        elif   score_final == up:
            i, j = n - 2, m - 1
            initial_align1 = seq1[i-2]  
            initial_align2 = seq2[j-1]
        elif score_final == left: 
            i, j = n - 1, m - 2 
            initial_align1 = seq1[i-1]  
            initial_align2 = seq2[j-2]   


        alignments = []  
        

        alignments = self.explore_alignments(i, j, initial_align1, initial_align2,alignments,score_matrix,seq1,seq2)
        
        #print("Matrice de score :")
        #for row in score_matrix:
        #    print(row)
        
        return alignments, score_final


    def calculate_score(self, align1, align2):
        score = 0
        for a, b in zip(align1, align2):
            if a == '-' or b == '-':
                score += self.gap
            elif a == b:
                score += self.match
            else:
                score += self.mismatch
        return score
    
    

    def calculate_distance(self, seq1, seq2):
        """
        Calcule la distance entre deux séquences en comptant le nombre de différences
        caractère par caractère jusqu'à la longueur de la séquence la plus courte.
        Les caractères restants de la séquence la plus longue sont ignorés.
        """
        # Longueur minimale des deux séquences
        min_len = min(len(seq1), len(seq2))
        
        # Comparer caractère par caractère jusqu'à la longueur minimale
        differences = sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a != b)
        
        return differences

    def compute_distance_matrix(self, sequences):
        """Construit une matrice des distances entre toutes les séquences."""
        n = len(sequences)
        distance_matrix = [[0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                distance = self.calculate_distance(sequences[i], sequences[j])
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance

        print(distance_matrix)
        return distance_matrix

    def progressive_alignment(self, sequences):
        """Effectue un alignement multiple progressif basé sur la matrice des distances."""
        n = len(sequences)
        distance_matrix = self.compute_distance_matrix(sequences)

        # Liste pour suivre les séquences alignées
        aligned_sequences = []
        unaligned_indices = list(range(n))

        # Étape 1 : Trouver les deux séquences les plus proches
        min_distance = float('inf')
        first, second = 0, 0
        for i in range(n):
            for j in range(i + 1, n):
                if distance_matrix[i][j] < min_distance:
                    min_distance = distance_matrix[i][j]
                    first, second = i, j

        # Aligner les deux premières séquences
        alignments, score_align = self.global_alignment(sequences[first], sequences[second])
        print("les premiers alignements")
        print(alignments)
        aligned1, aligned2 = alignments[0]
        print("aligned1 : ")
        print(aligned1)
        print("aligned2 : ")
        print(aligned2)
        #aligned_sequences.extend([(first, aligned1), (second, aligned2)])
        aligned_sequences.extend([aligned1, aligned2])
        print("aligned séquences for now :")
        print(aligned_sequences)
        unaligned_indices.remove(first)
        unaligned_indices.remove(second)
        print("les indices qui restent : ")
        print(unaligned_indices)

        # Étape 2 : Ajouter les autres séquences progressivement
        while unaligned_indices:
            best_candidate = None
            best_alignment = None
            min_distance = float('inf')

            for i in unaligned_indices:
                for aligned_seq in aligned_sequences:
                    distance = self.calculate_distance(sequences[i], aligned_seq)
                    if distance < min_distance:
                        min_distance = distance
                        best_candidate = i

            if best_candidate is not None:
                # Align only once after finding the best candidate
                alignments, score_align = self.global_alignment(sequences[best_candidate], aligned_sequences[0])
                aligned1, aligned2 = alignments[0]
                aligned_sequences.append(aligned1)
                unaligned_indices.remove(best_candidate)

        return aligned_sequences

# Exemple d'utilisation :
#alignment_handler = AlignmentHandler()
#sequences = ["ACTGG", "ACTTGG", "ACTGC", "CTTG"]  # Exemple de séquences
#aligned = alignment_handler.progressive_alignment(sequences)

#for seq in aligned:
#    print(seq)
