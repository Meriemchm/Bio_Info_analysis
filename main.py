from dbHandler import Handler
import os
from Aligner import AlignmentHandler  
from GeneticAlgorithm import GeneticAlgorithm
from visualize import Visualize
json_path = os.path.join(os.path.dirname(__file__), "Sars-cov-2.json")

handler = Handler(json_path)

resultat =handler.choisir_sequences_random()

print("\nTableau des séquences choisies :")
print(resultat)
visualize = Visualize()
alignment_handler = AlignmentHandler(visualize)


#---------------------test genetic---------------------------
# model = GeneticAlgorithm(AlignmentHandler,resultat)
# aligned = model.processing(15)
# print('les scores obtenu : ',aligned[0])
# print('meilleure parents : ',aligned[1])
# print('meilleure score atteint : ',aligned[2])


#---------------------ClustalW---------------------------

# aligned = alignment_handler.progressive_alignment(resultat)
# print('aligned sequences : ')
# for seq in aligned:
#     print(seq)

#---------------------Upgma---------------------------

# aligned = alignment_handler.upgma_processing(resultat)
# print('aligned sequences : ')
# for seq in aligned:
#     print(seq)

#---------------------NJ---------------------------

alignment_handler.NJ_processing(resultat)
