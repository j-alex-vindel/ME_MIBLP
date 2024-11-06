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

class Met_Net:
    """
    A class to represent a Metabolic Network

    ...
    Attributes
    ----------
    S : numpy array (NxM)
        A numpy array which contains the Stochiometric coefficients of the Metabolic Network (NxM)
    LB : array (1xM) 
        An array or list which contains the reactions' lower bound values of the metabolic network
    UB : array (1xM)
        An array or list which contains the reactions' upper bound values of the metabolic network
    Rxn : list[str]
        A list of strings with the names of the reactions in the Metabolic Network size (M)
    Met : list[str]
        A list of strings with the names of the metabolites in the Metabolic Network size (N)
    KO : list[int]
        A list of indeces of non essential reactions ie reactions cadidate for knockouts
    Name : str
        The name of the Metabolic Network, usually the name of the strain
    biomass : int
        The index of the biomass reaction (growth)
    chemical : int
        The index of the chemical of interest
    infeas : float
        The infeasibility parameter (default at 1e-6)
    time_limit : int
        A large number to account for the maximum computing time allowed inside the solver (default is 1000)
    BigM : int
        A large integer number for linearization and constraint construction (default is 1000)
    target : float
        Fraction threshold of the maximal possible biomass production rate (default to 50%)
    minprod : float
        Biomass threshold given by target*FBA[biomass]
    FBA : list[float]
        Reaction production rates under undisturbed Metabolic Network, often called widltype 
    FVA : list[float]
        Reaction production rate under undisturbed Metabolic Network, where the f.o. is the chemical of interest.
        And, the biomass reaction has a minprod threshold
    M : list[int]
        Set of indeces size (M)
    N : list[int]
        Set of indices size (N)
    b : list[int]
        Array of zeroes for the rhs in Sv=b 
    c : list[int]
        Array of obj coeff for the wildtype
    """
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
                 BM: Big_M = 1000):
        
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
        self.M = set_constructor(self.Rxn)
        self.N = set_constructor(self.Met)
        self.b = np.array([0 for i in self.N])
        self.c = np.array([1 if i == self.biomass else 0 for i in self.M])
        self.target = .5
        
    @property
    def target(self):
        return self._target

    @target.setter
    def target(self,value:Target=.5):
        self._minprod = None
        self._FBA = None
        # self._FVA = None
        self._target = value
    
    @property
    def chemical(self):
        return self._chemical

    @chemical.setter
    def chemical(self,value:Index):
        self._FVA = None
        self._chemical = value

    @property
    def minprod(self):
        if self._minprod is None:
            self._minprod = self._target*self.FBA[self.biomass]
        return self._minprod

    @property
    def FBA(self):
        if self._FBA is None:
            self._FBA = wildtype_FBA(self)
        return self._FBA
    @property
    def FVA(self):
        if self._FVA is None:
            self._FVA = wildtype_FBA(self,wildtype=False,mutant=True)
        return self._FVA
