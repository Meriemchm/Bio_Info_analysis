
import numpy as np
from Aligner import AlignmentHandler  
import json
import random
class GeneticAlgorithm():
    def __init__(self,alignment_handler=None,sequences = ["ACTGG", "ACTTGG", "CTTG", "ACTGC", "ATGCT"],kids=5,mutation=0.5):
        self.alignment_handler = alignment_handler() if callable(alignment_handler) else alignment_handler
        self.kids = kids
        self.mutation = mutation
        self.sequences = sequences
        
    def generate_gen(self, parents):

        first, second = parents[0], parents[1]

        m = len(first) 
        m2 = len(second)

        n_max = max(m, m2)
        align = []
        for i in range(n_max):
            n = np.random.choice([True, False])
            if i < m and n:
                align.append(first[i])
            elif i < m2 and not n:
                align.append(second[i])
            else:

                if i < m:
                    align.extend(first[i:])
                    break
                else:
                    align.extend(second[i:])
                    break
        return "".join(align)

          
    def mutation_process(self,seq):
        if np.random.choice([True,False],p=[self.mutation,1-self.mutation]):
            p_size = len(seq)
            start = np.random.randint(0,p_size)
            #print('start : ',start)
            choices = np.random.choice(['A','G',"T","C","-"],size=5)
            mutat = "".join(choices)
            #print(mutat)

            #print('seq',seq)
            seq = seq[:start]+mutat+seq[start+5:]
            #print('nvseq ',seq)
        return seq

    def best_process(self):
        aligned = self.alignment_handler.pair_par_pair(self.sequences,'k')
        resut = aligned[0]
        best_score = resut[-1]
        return resut[:-1],best_score
      
    def processing(self, iterations):
        scores = []  
        no_improvement = 0  
        patience = 5 

        for iter in range(iterations):
            print(f"Gen : {iter}")
            
            best_parents, best_score = self.best_process()
            print(f"meilleure score: {best_score}, meilleure Parents: {best_parents}")
            
            scores.append(best_score)
            
            new_gen = []
            print("nouvelle generation:")
            
            for gen in range(self.kids): 
                new_gen.append(self.mutation_process(self.generate_gen(best_parents)))
            
            self.sequences = new_gen + self.sequences  
            
            print(new_gen)

            if len(scores) > 1 and best_score <= scores[-2]: 
                no_improvement += 1
            else:
                no_improvement = 0 

            if no_improvement >= patience:
                print(f"Arrêt : Pas d'amélioration après {patience} générations")
                break

        return scores, best_parents, best_score
