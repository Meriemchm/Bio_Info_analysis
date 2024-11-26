from Bio import SeqIO, Entrez
import os
from itertools import combinations
from collections import defaultdict
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import numpy as np
from visualize import Visualize  

class AlignmentHandler:
    def __init__(self,visualize=None, match=2, mismatch=-1, gap=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.visualize = visualize() if callable(visualize) else visualize


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
            if i>1 and j>1 :
        
                    self.explore_alignments(i - 1, j - 1, seq1[i -2] + current_align1, seq2[j -2] + current_align2,alignments,score_matrix,seq1,seq2)
            elif j<1 : 
                self.explore_alignments(i - 1, j - 1, seq1[i -2] + current_align1, '-' + current_align2,alignments,score_matrix,seq1,seq2)
            elif i<1:
                self.explore_alignments(i - 1, j - 1, '-' + current_align1, seq2[j -2] + current_align2,alignments,score_matrix,seq1,seq2)
        if max_value == up:
            self.explore_alignments(i - 1, j, seq1[i - 2] + current_align1, '-' + current_align2,alignments,score_matrix,seq1,seq2)
        if max_value == left:
            self.explore_alignments(i, j - 1, '-' + current_align1, seq2[j - 2] + current_align2,alignments,score_matrix,seq1,seq2)

        return alignments    

    def global_alignment(self, seq1, seq2):

        n = len(seq1) + 1
        m = len(seq2) + 1

        score_matrix = [[0] * m for _ in range(n)]

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
        
        return alignments, score_matrix[n - 1][m - 1],score_matrix

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
            #print(row)
        
        return alignments, score_final,score_matrix


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



    def pair_par_pair(self, sequences,align='matrix', method='global'):
        n = len(sequences)
        alignment_scores = [[0] * n for _ in range(n)]
        alignement_best = []

        for i in range(n):
            for j in range(i + 1, n):
                if method == 'global':
                    alignement, score,_ = self.global_alignment(sequences[i], sequences[j])
                    alignments_with_scores = []
                    # print('matrix score : ',score) 
                    for align1 ,align2 in alignement:
                        score_alignement = self.calculate_score(align1, align2)
                        alignments_with_scores.append([align1, align2, score_alignement])
                    alignments_with_scores.sort(key=lambda x: x[2], reverse=True)

                    #retourne alignement a la pace de la matrice
                    if align != 'matrix':
                        new_align = alignments_with_scores[0] 
                        current_score = new_align[2]
                        # print(alignement_best)
                        if len(alignement_best) !=0:
                            # print(current_score,alignement_best[0][2])
                            if not alignement_best or current_score > int(alignement_best[0][2]):
                                alignement_best.clear()
                                alignement_best.append(new_align)
                        else :
                            alignement_best.append(new_align)

                    else :
                        print(alignments_with_scores[0])        
                elif method == 'local':
                    alignement, score,_ = self.local_alignment(sequences[i], sequences[j])
                    alignments_with_scores = []
                    for align1 ,align2 in alignement:
                        score_alignement = self.calculate_score(align1, align2)
                        alignments_with_scores.append([align1, align2, score_alignement])
                    alignments_with_scores.sort(key=lambda x: x[2], reverse=True)
                   
                alignment_scores[i][j] = score
                alignment_scores[j][i] = score

        if align == 'matrix':       
            print("\nMatrice d'alignement :")
            for row in alignment_scores:
                print(row)  
            return alignment_scores
        else :
            return alignement_best   


    def construct_guide_tree(self, sequences):
        # distance_matrix = self.pair_par_pair(sequences,'matrix','local')
        distance_matrix = self.pair_par_pair(sequences)
        num_seqs = len(distance_matrix)
        data = {i: [i] for i in range(num_seqs)}
        print(data)
        guide_tree = []

        while len(data) > 1:
            # max
            max_dist = float('-inf')
            max_pair = None

            for i in data:
                print(i)
                
                for j in data:
                   
                    if i < j:
                        print(data[i],data[j] )
                        for x in data[i] :
                            for y in data[j]:
                                print(x,y)
                                print(distance_matrix[x][y])
                        print("--------")   
                        distmean = np.mean([distance_matrix[x][y] for x in data[i] for y in data[j]])     
                        dist = max([distance_matrix[x][y] for x in data[i] for y in data[j]])
                        print(dist,distmean)
                        print("--------") 

                        if dist > max_dist: 
                            max_dist = dist
                            max_pair = (i, j)

            # Fusionner 
            i, j = max_pair
            new_data = data[i] + data[j]
            print(new_data)
            
            guide_tree.append((data[i], data[j]))
            print(guide_tree)

            # Supprime add data
            del data[i]
            del data[j]

            # add nv dt
            new_data_id = max(data.keys(), default=-1) + 1
            data[new_data_id] = new_data
        
        return guide_tree
    

    def progressive_alignment(self, sequences):
        """
        Aligne les séquences en suivant l'arbre guide.
        Chaque alignement est effectué uniquement avec la première séquence alignée du résultat précédent.
        """
        guide_tree = self.construct_guide_tree(sequences)
        print(guide_tree)
        self.visualize.visualize_tree(guide_tree)
        aligned_sequences = {} 
        final_alignments= [] 

        for step, (cluster1, cluster2) in enumerate(guide_tree):
            print(f"\nÉtape {step + 1} : Alignement des data {cluster1} et {cluster2}")

            seq1 = aligned_sequences[cluster1[0]] if cluster1[0] in aligned_sequences else sequences[cluster1[0]]
            seq2 = aligned_sequences[cluster2[0]] if cluster2[0] in aligned_sequences else sequences[cluster2[0]]

            alignments, _, _ = self.global_alignment(seq1, seq2)
            alignments_with_scores = []

            for align1, align2 in alignments:
                score_alignement = self.calculate_score(align1, align2)
                alignments_with_scores.append((align1, align2, score_alignement))

            alignments_with_scores.sort(key=lambda x: x[2], reverse=True)
            align1, align2 = alignments_with_scores[0][0], alignments_with_scores[0][1]
            if step == 0 :
                    final_alignments.append(align1)
                    final_alignments.append(align2)
            else :
                final_alignments.append(align1)


            print(f"Alignement retenu :\n{align1}\n{align2}")

            for i in cluster1:
                aligned_sequences[i] = align1  # add
            for i in cluster2:
                aligned_sequences[i] = align2  
        return final_alignments
    

    def calculate_distance(self, seq1, seq2):
        """
        Calcule la distance entre deux séquences.
        """
        min_len = min(len(seq1), len(seq2))
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

        print("Matrice des distances :")
        for row in distance_matrix:
            print(row)
           
        return distance_matrix

    def print_distance_matrix(self,matrix, labels):
        print("\nMatrice de distances:")
        print("   ", "  ".join(labels))
        for i, row in enumerate(matrix):
            print(f"{labels[i]:<3} {row}")

    # nvl matrice apr regroupement
    def upgma(self,matrix, data_a, data_b, new_label, labels):
        new_matrix = []
        new_labels = [label for i, label in enumerate(labels) if i not in (data_a, data_b)]
        new_labels.append(new_label)
        
        # nvl matrice de distances
        for i in range(len(matrix)):
            if i in (data_a, data_b):
                continue
            new_row = []
            for j in range(len(matrix)):
                if j in (data_a, data_b):
                    continue
                new_row.append(matrix[i][j])
            new_matrix.append(new_row)
        
        # new data cl , rw
        new_data_distances = []
        for i in range(len(matrix)):
            if i not in (data_a, data_b):
                dist = (matrix[i][data_a] + matrix[i][data_b]) / 2
                new_data_distances.append(dist)
        new_data_distances.append(0)  # dist luimme 0
        
        for i, row in enumerate(new_matrix):
            row.append(new_data_distances[i])
        new_matrix.append(new_data_distances)
        
        
        new_matrix = np.array(new_matrix)
        return new_matrix, new_labels

    def upgma_processing(self,sequences):
        
        dist_matrix = self.compute_distance_matrix(sequences)
        labels = [f"I{i+1}" for i in range(len(sequences))]
        step = 1
        while len(dist_matrix) > 1:

            print('len dist ',len(dist_matrix))
            # deux dataproches
            min_dist = float("inf")
            cluster_a, cluster_b = -1, -1
            for i in range(len(dist_matrix)):
                for j in range(i + 1, len(dist_matrix)):
                    if dist_matrix[i][j] < min_dist:
                        min_dist = dist_matrix[i][j]
                        cluster_a, cluster_b = i, j
            
            # Afficheage
            print(f"\nÉtape {step} : Regrouper {labels[cluster_a]} et {labels[cluster_b]} (distance = {min_dist})")
            new_label = f"({labels[cluster_a]},{labels[cluster_b]})"
            dist_matrix, labels = self.upgma(dist_matrix, cluster_a, cluster_b, new_label, labels)
            self.print_distance_matrix(dist_matrix, labels)
            step += 1

        print("\nClustering terminé!")
        print(labels[0])
        self.visualize.tree_vis(labels[0])


    # R_i 
    def calculate_Ri(self,matrix):
        n = len(matrix)
        Ri = np.sum(matrix, axis=1) / (n - 2)
        return Ri

    # S_ij = dij - ri -rj
    def calculate_Sij(self,matrix, Ri):
        n = len(matrix)
        Sij = np.zeros_like(matrix)
        for i in range(n):
            for j in range(n):
                if i != j:
                    Sij[i][j] = matrix[i][j] - Ri[i] - Ri[j]
        return Sij

    # Maj 
    def NJ(self,matrix, cluster_a, cluster_b):
        n = len(matrix)
        new_matrix = []
        for i in range(n):
            if i != cluster_a and i != cluster_b:
                new_row = []
                for j in range(n):
                    if j != cluster_a and j != cluster_b:
                        new_row.append(matrix[i][j])
                new_matrix.append(new_row)
        # nvl distances
        new_cluster_distances = []
        for i in range(n):
            if i not in (cluster_a, cluster_b):
                dist = (matrix[i][cluster_a] + matrix[i][cluster_b]-matrix[cluster_a][cluster_b]) / 2
                new_cluster_distances.append(dist)
        new_cluster_distances.append(0)  # Distance à lui-même

        for i, row in enumerate(new_matrix):
            row.append(new_cluster_distances[i])
        new_matrix.append(new_cluster_distances)
        return np.array(new_matrix)

    def NJ_processing(self,sequences):
        dist_matrix = np.array(self.compute_distance_matrix(sequences))
        labels = [f"I{i+1}" for i in range(len(sequences))]
        step = 1
        while len(dist_matrix) > 1:

            print(f"\n=== Étape {step} ===")
            print("\nMatrice de distances:")
            print("\t" + "\t".join(labels))
            for label, row in zip(labels, dist_matrix):
                print(f"{label}\t{row}")
            if len(dist_matrix) > 2:    
                # Calcul de R_i
                Ri = self.calculate_Ri(dist_matrix)
                print(f"\nR_i = {Ri}")
                
                # Calcul de S_ij
                Sij = self.calculate_Sij(dist_matrix, Ri)
                print("\nS_ij =")
                print(Sij)
                
                # data proche(min Sij)
                min_value = np.inf
                cluster_a, cluster_b = -1, -1
                for i in range(len(Sij)):
                    for j in range(i + 1, len(Sij)):
                        if Sij[i][j] < min_value:
                            min_value = Sij[i][j]
                            cluster_a, cluster_b = i, j
            if len(dist_matrix) == 2:
                cluster_a, cluster_b = 0, 1
                min_value = None
            print(f"Regroupement des data {labels[cluster_a]} et {labels[cluster_b]} (S_ij = {min_value})")
            
            if len(dist_matrix) > 2:
                # Calcul des nouvelles distances
                d_a = (dist_matrix[cluster_a, cluster_b] + Ri[cluster_a] - Ri[cluster_b]) / 2
                d_b = (dist_matrix[cluster_a, cluster_b] + Ri[cluster_b] - Ri[cluster_a]) / 2
                print(f"Distance {labels[cluster_a]} -> arbre = {d_a}")
                print(f"Distance {labels[cluster_b]} -> arbre = {d_b}")
            
            new_label = f"({labels[cluster_a]},{labels[cluster_b]})"
            labels = [label for i, label in enumerate(labels) if i not in (cluster_a, cluster_b)] + [new_label]
            dist_matrix = self.NJ(dist_matrix, cluster_a, cluster_b)
            
            step += 1

        print("\nClustering terminé!")
        print(labels[0])
        self.visualize.tree_vis(labels[0])

# visualize = Visualize()
# alignment_handler = AlignmentHandler(visualize)
# sequences = ["ACTGG", "ACTTGG", "CTTG", "ACTGC", "ATGCT"]  
# print("sequences : ",sequences)
# # alignment_handler.upgma_processing(sequences)
# aligned = alignment_handler.pair_par_pair(sequences,'k')
# # print('aligned sequences : ')
# # for seq in aligned:
# #     print(seq)
# # Visualize.visualize_tree(sequences)