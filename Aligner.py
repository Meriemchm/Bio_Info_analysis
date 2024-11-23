from Bio import SeqIO, Entrez
import os
from itertools import combinations
from collections import defaultdict
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from ete3 import Tree
from Bio import Phylo
from io import StringIO

class AlignmentHandler:
    def __init__(self, match=2, mismatch=-1, gap=-1):
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



    def pair_par_pair(self, sequences, method='global'):
        n = len(sequences)
        alignment_scores = [[0] * n for _ in range(n)]

        for i in range(n):
            for j in range(i + 1, n):
                if method == 'global':
                    alignement, score,_ = self.global_alignment(sequences[i], sequences[j])
                    alignments_with_scores = []
                    for align1 ,align2 in alignement:
                        score_alignement = self.calculate_score(align1, align2)
                        alignments_with_scores.append((align1, align2, score_alignement))
                    alignments_with_scores.sort(key=lambda x: x[2], reverse=True)
                    print(alignments_with_scores[0])
                elif method == 'local':
                    _, score,_ = self.local_alignment(sequences[i], sequences[j])
                alignment_scores[i][j] = score
                alignment_scores[j][i] = score

        print("\nMatrice d'alignement :")
        for row in alignment_scores:
            print(row)     
        return alignment_scores   


    def construct_guide_tree(self, sequences):
        distance_matrix = self.pair_par_pair(sequences)
        num_seqs = len(distance_matrix)
        clusters = {i: [i] for i in range(num_seqs)}
        print(clusters)
        guide_tree = []

        while len(clusters) > 1:
            # max
            max_dist = float('-inf')
            max_pair = None

            for i in clusters:
                print(i)
                
                for j in clusters:
                   
                    if i < j:
                        # Calculer la distance moyenne entre les clusters i et j
                        # dist = np.mean([distance_matrix[x][y] for x in clusters[i] for y in clusters[j]])
                        print(clusters[i],clusters[j] )
                        for x in clusters[i] :
                            for y in clusters[j]:
                                print(x,y)
                                print(distance_matrix[x][y])
                        print("--------")   
                        distmean = np.mean([distance_matrix[x][y] for x in clusters[i] for y in clusters[j]])     
                        dist = max([distance_matrix[x][y] for x in clusters[i] for y in clusters[j]])
                        print(dist,distmean)
                        print("--------") 

                        if dist > max_dist: 
                            max_dist = dist
                            max_pair = (i, j)

            # Fusionner 
            i, j = max_pair
            new_cluster = clusters[i] + clusters[j]
            
            guide_tree.append((clusters[i], clusters[j]))

            # Supprimer les anciens clusters
            del clusters[i]
            del clusters[j]

            # Ajouter le nouveau cluster
            new_cluster_id = max(clusters.keys(), default=-1) + 1
            clusters[new_cluster_id] = new_cluster
        print(guide_tree)
        return guide_tree
    

    def progressive_alignment(self, sequences):
        """
        Aligne les séquences en suivant l'arbre guide.
        Chaque alignement est effectué uniquement avec la première séquence alignée du résultat précédent.
        """
        guide_tree = self.construct_guide_tree(sequences)
        aligned_sequences = {} 
        final_alignments= [] 

        for step, (cluster1, cluster2) in enumerate(guide_tree):
            print(f"\nÉtape {step + 1} : Alignement des clusters {cluster1} et {cluster2}")

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
                aligned_sequences[i] = align1  # Mettre à jour avec align1
            for i in cluster2:
                aligned_sequences[i] = align2  # Mettre à jour avec align2
        return final_alignments
    

    def visualize_tree(self,sequences):
        """
        Visualise un arbre à partir de clusters.
        :param clusters: Liste de tuples représentant les clusters, ex: [([0], [1]), ([2], [0, 1]), ([3], [2, 0, 1])]
        """

        clusters = self.construct_guide_tree(sequences)
        G = nx.DiGraph()

        cluster_to_node = {}

        for i, (left, right) in enumerate(clusters):
            current_node = f"Align {i + 1}"
            G.add_node(current_node)
            
            if tuple(left) not in cluster_to_node:
                cluster_to_node[tuple(left)] = f"Seq {', '.join(map(str, left))}"
            if tuple(right) not in cluster_to_node:
                cluster_to_node[tuple(right)] = f"Seq {', '.join(map(str, right))}"

            G.add_edge(cluster_to_node[tuple(left)], current_node)
            G.add_edge(cluster_to_node[tuple(right)], current_node)
            
            cluster_to_node[tuple(left + right)] = current_node

        try:
            pos = nx.nx_agraph.graphviz_layout(G, prog='dot')  
        except ImportError:
            pos = nx.spring_layout(G, seed=25)

        plt.figure(figsize=(12, 8))
        nx.draw(
            G,
            pos,
            with_labels=True,
            node_size=2000,
            node_color="lightgreen",
            font_size=10,
            font_weight="bold",
        )
        plt.title("Arbre Guide - Alignements")
        plt.show()



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
    def upgma(self,matrix, cluster_a, cluster_b, new_label, labels):
        new_matrix = []
        new_labels = [label for i, label in enumerate(labels) if i not in (cluster_a, cluster_b)]
        new_labels.append(new_label)
        
        # nvl matrice de distances
        for i in range(len(matrix)):
            if i in (cluster_a, cluster_b):
                continue
            new_row = []
            for j in range(len(matrix)):
                if j in (cluster_a, cluster_b):
                    continue
                new_row.append(matrix[i][j])
            new_matrix.append(new_row)
        
        # Ajouter une nouvelle ligne/colonne pour le nouveau cluster
        new_cluster_distances = []
        for i in range(len(matrix)):
            if i not in (cluster_a, cluster_b):
                dist = (matrix[i][cluster_a] + matrix[i][cluster_b]) / 2
                new_cluster_distances.append(dist)
        new_cluster_distances.append(0)  # La distance à lui-même est 0
        
        # Insérer la nouvelle ligne/colonne
        for i, row in enumerate(new_matrix):
            row.append(new_cluster_distances[i])
        new_matrix.append(new_cluster_distances)
        
        
        new_matrix = np.array(new_matrix)
        return new_matrix, new_labels

    def upgma_processing(self,sequences):
        
        dist_matrix = self.compute_distance_matrix(sequences)
        labels = [f"I{i+1}" for i in range(len(sequences))]
        step = 2
        while len(dist_matrix) > 1:
            # deux clustersproches
            min_dist = float("inf")
            cluster_a, cluster_b = -1, -1
            for i in range(len(dist_matrix)):
                for j in range(i + 1, len(dist_matrix)):
                    if dist_matrix[i][j] < min_dist:
                        min_dist = dist_matrix[i][j]
                        cluster_a, cluster_b = i, j
            
            # Afficher l'étape
            print(f"\nÉtape {step} : Regrouper {labels[cluster_a]} et {labels[cluster_b]} (distance = {min_dist})")
            new_label = f"({labels[cluster_a]},{labels[cluster_b]})"
            dist_matrix, labels = self.upgma(dist_matrix, cluster_a, cluster_b, new_label, labels)
            self.print_distance_matrix(dist_matrix, labels)
            step += 1

        print("\nClustering terminé!")
        print(labels[0])
        self.tree_vis(labels[0])


    # R_i 
    def calculate_Ri(self,matrix):
        n = len(matrix)
        Ri = np.sum(matrix, axis=1) / (n - 2)
        return Ri

    # S_ij
    def calculate_Sij(self,matrix, Ri):
        n = len(matrix)
        Sij = np.zeros_like(matrix)
        for i in range(n):
            for j in range(n):
                if i != j:
                    Sij[i][j] = matrix[i][j] - Ri[i] - Ri[j]
        return Sij

    # Maj matrice après regroupement
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
                dist = (matrix[i][cluster_a] + matrix[i][cluster_b]) / 2
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
            
            # Calcul de R_i
            Ri = self.calculate_Ri(dist_matrix)
            print(f"\nR_i = {Ri}")
            
            # Calcul de S_ij
            Sij = self.calculate_Sij(dist_matrix, Ri)
            print("\nS_ij =")
            print(Sij)
            
            # clusters les plus proches (min Sij)
            min_value = np.inf
            cluster_a, cluster_b = -1, -1
            for i in range(len(Sij)):
                for j in range(i + 1, len(Sij)):
                    if Sij[i][j] < min_value:
                        min_value = Sij[i][j]
                        cluster_a, cluster_b = i, j

            print(f"Regroupement des clusters {labels[cluster_a]} et {labels[cluster_b]} (S_ij = {min_value})")
            
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
        self.tree_vis(labels[0])

    def tree_vis(self,newick):
        handle = StringIO(newick)
        tree = Phylo.read(handle, "newick")

        # Dessiner l'arbre dans une fenêtre graphique
        Phylo.draw(tree)

        # Ou dessiner dans une console texte (arbre ASCII)
        Phylo.draw_ascii(tree)

alignment_handler = AlignmentHandler()
sequences = ["ACTGG", "ACTTGG", "CTTG", "ACTGC", "ATGCT"]  
print("sequences : ",sequences)
alignment_handler.upgma_processing(sequences)

# print('aligned sequences : ')
# for seq in aligned:
#     print(seq)
# alignment_handler.visualize_tree(sequences)