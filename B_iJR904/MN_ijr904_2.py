# Checking the MN iJR904 with the LB and UB values from the Reack Knock

# First check the FBA value if we ran it witouht changing the values in the LB,UB just the glc uptake

import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES')))

from Ob_Met_Net import Metabolic_Network,Met_Net_2
from pymatreader import read_mat
from typing import List, NewType
from support_functions import set_constructor,wildtype_FBA
import numpy as np

## ============================================================================================= ====================================================================================
ijr904 = 'iJR904'

data = read_mat(f"../B_iJR904/Data/iJR904.mat")

LB = data[ijr904]['lb'].tolist()
UB = data[ijr904]['ub'].tolist()
Met = data[ijr904]['mets']
Rxn = data[ijr904]['rxns']
S = data[ijr904]['S']

biomas = Rxn.index('BIOMASS_Ecoli')
chemical = Rxn.index('EX_succ_e')

LB[Rxn.index('EX_glc__D_e')] = -10
UB[Rxn.index('EX_glc__D_e')] = -10




mn = Met_Net_2(S=S,LB=LB,UB=UB,Rxn=Rxn,Met=Met,biomass=biomas,chemical=chemical)

print(f"FBA -> {mn.FBA}")

print(f"Target -> {mn.target}")

print(f"Min prod -> {mn.minprod}")
