import sys
import os
import matplotlib.pyplot as plt
import pandas as pd

# Add the directory containing your modules to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from DataStats import DataStats
from Aligner import AlignmentHandler
from dbHandler import Handler
from GeneticAlgorithm import GeneticAlgorithm

def main():

    def visualize_performance(csv_path):
        # Charger les données à partir du fichier CSV
        df = pd.read_csv(csv_path)

        # Exemple 1 : Graphe d'exécution par méthode
        plt.figure(figsize=(10, 6))
        for method in df['Method'].unique():
            subset = df[df['Method'] == method]
            plt.plot(subset['Execution Time (s)'], marker='o', label=method)

        plt.title('Temps d\'exécution par méthode et par année')
        #plt.xlabel('Année')
        plt.ylabel('Temps d\'exécution (secondes)')
        plt.legend(title="Méthode")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

        # Exemple 2 : Mémoire utilisée par méthode
        plt.figure(figsize=(10, 6))
        for method in df['Method'].unique():
            subset = df[df['Method'] == method]
            plt.bar( subset['Memory Used (MB)'], marker='o',label=method)

        plt.title('Mémoire utilisée par méthode et par année')
        plt.xlabel('Année')
        plt.ylabel('Mémoire utilisée (MB)')
        plt.legend(title="Méthode")
        plt.tight_layout()
        plt.show()

        # Exemple 3 : Conservation moyenne par variante
        plt.figure(figsize=(10, 6))
        variant_conservation = df.groupby('Variant')['Conservation Score'].mean()
        variant_conservation.plot(kind='bar', color='skyblue')

        plt.title('Score de conservation moyen par variante')
        plt.xlabel('Variante')
        plt.ylabel('Score de conservation')
        plt.grid(axis='y')
        plt.tight_layout()
        plt.show()
    # Initialize necessary components
    json_path = os.path.join(os.path.dirname(__file__), "Sars-cov-2.json") 
    clustalw_path = "C://Program Files (x86)//ClustalW2//clustalw2.exe" 

    handler = Handler(json_path)
    alignment_handler = AlignmentHandler()

    # Create DataStats instance
    data_stats = DataStats(
        handler=handler,
        alignment_handler=alignment_handler,
        genetic_algorithm=GeneticAlgorithm,
        clustalw_path=clustalw_path
    )

    #variant = "Variant Beta"  
    #year = "2024"      

    #variants = ["Variant Alpha", "Variant Beta", "Variant Gamma", "Variant Delta"]
    variants = ["Variant Gamma"]
    years = ["2022", "2023", "2024"]
    methods = ["Pair", "ClustalW", "Genetic Algorithm"]

    # Run the performance test for a single method
    #performance_df = data_stats.run_single_method_test(variant, year, "ClustalW")
    #performance_df = data_stats.run_single_method_test(variant, year, "Genetic Algorithm")
    #performance_df = data_stats.run_single_method_test(variant, year, "Pair")

    performance_df = data_stats.run_comprehensive_tests(variants, years, methods, 3, 0, 15)
    print(performance_df)

    if performance_df is not None:
        # Save the DataFrame to a CSV file
        output_path = "performance_results_gamma.csv"
        performance_df.to_csv(output_path, index=False)
        print(f"Performance results saved to {output_path}")
    else:
        print("No performance data to save.")

    csv_file_path = "performance_results_15.csv"
    visualize_performance(csv_file_path)
    
    

if __name__ == "__main__":
    main()
