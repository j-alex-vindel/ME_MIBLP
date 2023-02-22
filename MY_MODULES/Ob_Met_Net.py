from typing import List, NewType
from support_functions import set_constructor,wildtype_FBA
import numpy as np


S_Matrix = List[List[int]]
Lower_Bound = List[int]
Upper_Bound = List[int]
Reactions = List[str]
Metabolites = List[str]
Knockouts = List[int]
Index = NewType('Index',int)
Big_M = NewType('Big M', int)
Target = NewType('Prod Target %',float)

class FBAError(Exception):
    pass

class Metabolic_Network:

    def __init__(self, 
                 S:S_Matrix = None, 
                 LB:Lower_Bound = None, 
                 UB:Upper_Bound = None, 
                 Rxn: Reactions = None, 
                 Met: Metabolites = None,
                 KO: Knockouts = None, 
                 Name: str = None, 
                 biomass: Index = None, 
                 chemical: Index = None, 
                 infeas:float = 1e-6, 
                 time_limit: int = 1000, 
                 BM: Big_M = 1000, 
                 target: Target = .5):
        
        self.S = S
        self.LB = LB
        self.UB = UB
        self.Rxn = Rxn
        self.Met = Met
        self.KO = KO
        self.Name = Name
        self.biomass = biomass
        self.chemical = chemical
        self.infeas = infeas
        self.time_limit = time_limit
        self.BM = BM
        self.target = target  
        self.M = set_constructor(self.Rxn)
        self.N = set_constructor(self.Met)
        self.b = np.array([0 for i in self.N])
        self.c = np.array([1 if i == self.biomass else 0 for i in self.M])
        self.FBA = wildtype_FBA(self)
        self.minprod = target*self.FBA[self.biomass]



