import json
import tkinter as tk
from tkinter import ttk, messagebox
from Aligner import AlignmentHandler  

class AlignmentGUI:
    def __init__(self, root, json_path, alignment_handler):
        self.root = root
        self.json_path = json_path
        self.alignment_handler = alignment_handler
        self.data = self.load_data()

        self.create_widgets()

    def load_data(self):
        try:
            with open(self.json_path, "r") as file:
                return json.load(file)
        except FileNotFoundError:
            messagebox.showerror("Erreur", f"Le fichier {self.json_path} est introuvable.")
            return {}

    def create_widgets(self):
        self.root.title("Interface d'Alignement (ɔ◔‿◔)ɔ ♥")

        # Choix du variant pour la première séquence
        ttk.Label(self.root, text="Séquence 1 - Choisissez un variant :").grid(row=0, column=0, padx=10, pady=5, sticky="w")
        self.variant1_combobox = ttk.Combobox(self.root, values=list(self.data.keys()))
        self.variant1_combobox.grid(row=0, column=1, padx=10, pady=5)
        self.variant1_combobox.bind("<<ComboboxSelected>>", self.update_years1)

        # Choix de l'année pour la première séquence
        ttk.Label(self.root, text="Séquence 1 - Choisissez une année :").grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.year1_combobox = ttk.Combobox(self.root)
        self.year1_combobox.grid(row=1, column=1, padx=10, pady=5)
        self.year1_combobox.bind("<<ComboboxSelected>>", self.update_sequences1)

        # Choix de la séquence spécifique pour la première séquence
        ttk.Label(self.root, text="Séquence 1 - Sélectionnez une séquence :").grid(row=2, column=0, padx=10, pady=5, sticky="w")
        self.sequence1_combobox = ttk.Combobox(self.root)
        self.sequence1_combobox.grid(row=2, column=1, padx=10, pady=5)

        # Entrées pour les positions de début et de fin de la première séquence
        ttk.Label(self.root, text="Séquence 1 - Position de départ :").grid(row=3, column=0, padx=10, pady=5, sticky="w")
        self.start1_entry = ttk.Entry(self.root)
        self.start1_entry.grid(row=3, column=1, padx=10, pady=5)

        ttk.Label(self.root, text="Séquence 1 - Position de fin :").grid(row=4, column=0, padx=10, pady=5, sticky="w")
        self.end1_entry = ttk.Entry(self.root)
        self.end1_entry.grid(row=4, column=1, padx=10, pady=5)

        # Choix du variant pour la deuxième séquence
        ttk.Label(self.root, text="Séquence 2 - Choisissez un variant :").grid(row=5, column=0, padx=10, pady=5, sticky="w")
        self.variant2_combobox = ttk.Combobox(self.root, values=list(self.data.keys()))
        self.variant2_combobox.grid(row=5, column=1, padx=10, pady=5)
        self.variant2_combobox.bind("<<ComboboxSelected>>", self.update_years2)

        # Choix de l'année pour la deuxième séquence
        ttk.Label(self.root, text="Séquence 2 - Choisissez une année :").grid(row=6, column=0, padx=10, pady=5, sticky="w")
        self.year2_combobox = ttk.Combobox(self.root)
        self.year2_combobox.grid(row=6, column=1, padx=10, pady=5)
        self.year2_combobox.bind("<<ComboboxSelected>>", self.update_sequences2)

        # Choix de la séquence spécifique pour la deuxième séquence
        ttk.Label(self.root, text="Séquence 2 - Sélectionnez une séquence :").grid(row=7, column=0, padx=10, pady=5, sticky="w")
        self.sequence2_combobox = ttk.Combobox(self.root)
        self.sequence2_combobox.grid(row=7, column=1, padx=10, pady=5)

        # Entrées pour les positions de début et de fin de la deuxième séquence
        ttk.Label(self.root, text="Séquence 2 - Position de départ :").grid(row=8, column=0, padx=10, pady=5, sticky="w")
        self.start2_entry = ttk.Entry(self.root)
        self.start2_entry.grid(row=8, column=1, padx=10, pady=5)

        ttk.Label(self.root, text="Séquence 2 - Position de fin :").grid(row=9, column=0, padx=10, pady=5, sticky="w")
        self.end2_entry = ttk.Entry(self.root)
        self.end2_entry.grid(row=9, column=1, padx=10, pady=5)

        # Choix du type d'alignement
        ttk.Label(self.root, text="Type d'alignement :").grid(row=10, column=0, padx=10, pady=5, sticky="w")
        self.alignment_type = tk.StringVar(value="global")
        ttk.Radiobutton(self.root, text="Global", variable=self.alignment_type, value="global").grid(row=10, column=1, sticky="w")
        ttk.Radiobutton(self.root, text="Local", variable=self.alignment_type, value="local").grid(row=10, column=1, sticky="e")

        # Bouton pour lancer l'alignement
        self.run_button = ttk.Button(self.root, text="Lancer l'alignement", command=self.run_alignment)
        self.run_button.grid(row=11, column=0, columnspan=2, pady=10)

        self.clear_button = tk.Button(self.root, text="Clear", command=self.clear_fields)
        self.clear_button.grid(row=11, column=1, columnspan=2, pady=10)  # Utilisation de grid ici aussi
        # Zone pour afficher les résultats
        self.result_text = tk.Text(self.root, height=15, width=60)
        self.result_text.grid(row=12, column=0, columnspan=2, pady=10, padx=10)

    def update_years1(self, event):
        variant = self.variant1_combobox.get()
        if variant in self.data:
            years = list(self.data[variant].keys())
            self.year1_combobox['values'] = years

    def update_years2(self, event):
        variant = self.variant2_combobox.get()
        if variant in self.data:
            years = list(self.data[variant].keys())
            self.year2_combobox['values'] = years

    def update_sequences1(self, event):
        variant = self.variant1_combobox.get()
        year = self.year1_combobox.get()
        if variant and year in self.data[variant]:
            sequences = list(self.data[variant][year].values())
            self.sequence1_combobox['values'] = sequences

    def update_sequences2(self, event):
        variant = self.variant2_combobox.get()
        year = self.year2_combobox.get()
        if variant and year in self.data[variant]:
            sequences = list(self.data[variant][year].values())
            self.sequence2_combobox['values'] = sequences

    def run_alignment(self):
        seq1 = self.sequence1_combobox.get()
        seq2 = self.sequence2_combobox.get()
        align_type = self.alignment_type.get()

        if not seq1 or not seq2:
            messagebox.showerror("Erreur", "Veuillez sélectionner des séquences pour l'alignement.")
            return
        
        try:
            start1 = int(self.start1_entry.get() or 1) - 1
            end1 = int(self.end1_entry.get() or len(seq1))
            start2 = int(self.start2_entry.get() or 1) - 1
            end2 = int(self.end2_entry.get() or len(seq2))
        except ValueError:
            messagebox.showerror("Erreur", "Les positions de départ et de fin doivent être des nombres entiers.")
            return

        # Slice the sequences based on user input
        seq1 = seq1[start1:end1]
        seq2 = seq2[start2:end2]


        if align_type == "global":
            alignments, score = self.alignment_handler.global_alignment(seq1, seq2)
        else:
            alignments, score = self.alignment_handler.local_alignment(seq1, seq2)

        self.display_results(alignments, score)

    def display_results(self, alignments, score):
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, f"Score d'alignement : {score}\n\n")
        for align1, align2 in alignments:
            self.result_text.insert(tk.END, f"{align1}\n{align2}\n\n")

    def clear_fields(self):
        """Effacer les champs de texte et autres éléments"""
        self.variant1_combobox.delete(0, tk.END)  
        self.variant2_combobox.delete(0, tk.END)  
        self.year1_combobox.delete(0, tk.END)  
        self.year2_combobox.delete(0, tk.END)  
        self.sequence1_combobox.delete(0, tk.END)  
        self.sequence2_combobox.delete(0, tk.END) 
        self.end1_entry.delete(0, tk.END)  
        self.start1_entry.delete(0, tk.END) 
        self.end2_entry.delete(0, tk.END)  
        self.start2_entry.delete(0, tk.END) 
        self.result_text.delete('1.0', tk.END)
        

# Utilisation
if __name__ == "__main__":
    json_path = "Sars-cov-2.json"
    alignment_handler = AlignmentHandler(match=2, mismatch=-1, gap=-2)

    root = tk.Tk()
    app = AlignmentGUI(root, json_path, alignment_handler)
    root.mainloop()
