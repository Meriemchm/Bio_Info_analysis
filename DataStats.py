import time
import psutil
import pandas as pd
import numpy as np
from ClustalW import ClustalHandler
from Bio import AlignIO
import os

class DataStats:
    def __init__(self, handler, alignment_handler, genetic_algorithm, clustalw_path):
        self.handler = handler
        self.alignment_handler = alignment_handler
        self.genetic_algorithm_class = genetic_algorithm
        self.performances = []
        self.clustalw_path = clustalw_path

    def measure_performance(self, sequences, method_name, variant):
        process = psutil.Process()
        
        # Prepare temporary fasta file for ClustalW input
        temp_input_fasta = "temp_input.fasta"
        self.handler.export_to_fasta(variant, temp_input_fasta)  # Pass variant and output file
        
        # Start performance tracking
        start_time = time.time()
        start_cpu = process.cpu_percent()
        start_memory = process.memory_info().rss

        try:
            # Execute the alignment method
            if method_name == "Genetic Algorithm":
                model = self.genetic_algorithm_class(self.alignment_handler, sequences)
                scores, aligned_sequences, best_score = model.processing(15)
            elif method_name == "ClustalW":
                # Use ClustalW alignment method
                aligned_sequences = self.alignment_handler.progressive_alignment(sequences)
            
            elif method_name == "NJ":
                aligned_sequences = self.alignment_handler.NJ_processing(sequences)

            elif method_name == "UPGMA":
                aligned_sequences = self.alignment_handler.upgma_processing(sequences)

            elif method_name =="Pair":
                aligned_sequences = self.alignment_handler.pair_par_pair(sequences)
            else:
                # Use other methods
                print("cette méthode n'est pas implémentée")

            # End performance tracking
            end_time = time.time()
            end_cpu = process.cpu_percent()
            end_memory = process.memory_info().rss

            # Calculate metrics
            execution_time = end_time - start_time
            cpu_usage = end_cpu
            memory_used = (end_memory - start_memory) / (1024 * 1024)  # Convert to MB

            # Calculate alignment quality metrics
            alignment_length = max(len(seq) for seq in aligned_sequences)
            gap_percentage = self._calculate_gap_percentage(aligned_sequences)
            conservation_score = self._calculate_conservation_score(aligned_sequences)

            # Cleanup temporary files
            if os.path.exists(temp_input_fasta):
                os.remove(temp_input_fasta)
            if os.path.exists("output.aln"):
                os.remove("output.aln")

            # Create performance record
            performance_record = {
                "Method": method_name,
                "Execution Time (s)": execution_time,
                "CPU Usage (%)": cpu_usage,
                "Memory Used (MB)": memory_used,
                "Alignment Length": alignment_length,
                "Gap Percentage (%)": gap_percentage,
                "Conservation Score": conservation_score
            }

            if all(value is not None for value in performance_record.values()):
                self.performances.append(performance_record)
            else:
                print(f"Incomplete performance record for {method_name}")

            return performance_record

        except Exception as e:
            print(f"Error in {method_name} performance measurement: {e}")
            return None

    def _calculate_gap_percentage(self, aligned_sequences):
        total_chars = sum(len(seq) for seq in aligned_sequences)
        total_gaps = sum(seq.count('-') for seq in aligned_sequences)
        return (total_gaps / total_chars) * 100 if total_chars > 0 else 0

    def _calculate_conservation_score(self, aligned_sequences):
        if not aligned_sequences:
            return 0

        # Transpose the alignment
        alignment_columns = list(zip(*aligned_sequences))
        
        # Calculate conservation for each column
        conservation_columns = []
        for column in alignment_columns:
            # Count unique non-gap characters
            unique_chars = set(char for char in column if char != '-')
            if len(unique_chars) == 1:
                conservation_columns.append(1)
            elif len(unique_chars) > 1:
                conservation_columns.append(0)

        # Return average conservation
        return np.mean(conservation_columns) if conservation_columns else 0

    def run_single_method_test(self, variant, year, method_name, num_sequences=3, start=0, end=5):
        # Get sequences
        sequences = self.handler.get_sequences_by_variant_and_year(variant, year)
        if isinstance(sequences, dict):
            sequences = [seq[start:end] for seq in list(sequences.values())[:num_sequences]]  # Limit the number of sequences
            sequences = [str(seq) if isinstance(seq, (str, bytes)) else "" for seq in sequences]

            if len(sequences) < num_sequences:
                print(f"Not enough sequences for variant '{variant}' in {year}.")
                return None

            # Run performance measurement for the selected method
            result = self.measure_performance(sequences, method_name, variant)

            if result:
                # Create a DataFrame for a single method's results
                df = pd.DataFrame([result])
                return df

        else:
            print(sequences)  # Error message
            return None
        

    def run_comprehensive_tests(self, variants, years, methods, num_sequences=3, start=0, end=5):
        all_results = []

        for variant in variants:
            for year in years:
                for method_name in methods:
                    print(f"Testing {method_name} for {variant} in {year}")
                    result = self.run_single_method_test(variant, year, method_name, num_sequences, start, end)
                    if result is not None:
                        all_results.append(result)

        # Combine all results into a single DataFrame
        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            return combined_df
        else:
            print("No results to combine.")
            return None
