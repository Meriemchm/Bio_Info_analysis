import networkx as nx
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo


class Visualize:

    def visualize_tree(self, data):
        """
        Visualise un arbre à partir de data.
        :param sequences: Liste de séquences utilisées pour construire l'arbre.
        """
        G = nx.DiGraph()

        cluster_to_node = {}

        for i, (left, right) in enumerate(data):
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

    def tree_vis(self, newick):
        """
        Visualise un arbre phylogénétique à partir d'une chaîne Newick.
        :param newick: Chaîne au format Newick représentant l'arbre.
        """
        handle = StringIO(newick)
        tree = Phylo.read(handle, "newick")

        # Dessin graphique
        Phylo.draw(tree)

        # Dessin ASCII
        print("\nReprésentation ASCII de l'arbre :\n")
        Phylo.draw_ascii(tree)


