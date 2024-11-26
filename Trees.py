import numpy as np
from Aligner import AlignmentHandler  
class GeneticAlgorithm():
    def __init__(self,alignment_handler=None,sequences = ["ACTGG", "ACTTGG", "CTTG", "ACTGC", "ATGCT"]    ,kids=5,mutation=0.5):
        self.alignment_handler = alignment_handler() if callable(alignment_handler) else alignment_handler
        self.kids = kids
        self.mutation = mutation
        self.sequences = sequences