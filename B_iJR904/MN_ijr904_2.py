# Checking the MN iJR904 with the LB and UB values from the Reack Knock

# First check the FBA value if we ran it witouht changing the values in the LB,UB just the glc uptake

import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES')))

from Ob_Met_Net import Metabolic_Network,Met_Net
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

# LB[Rxn.index('EX_glc__D_e')] = -10
# UB[Rxn.index('EX_glc__D_e')] = -10

UB[Rxn.index("EX_o2_e")] = -20
LB[Rxn.index("EX_o2_e")] = -20

LB[Rxn.index("ATPM")] = 7.6
UB[Rxn.index("ATPM")] = 7.6

MN = Metabolic_Network(S=S,LB=LB,UB=UB,Met=Met,Rxn=Rxn,biomass=biomas,chemical=chemical)
print(f"Conventional Object")
print(f"MN FBA {MN.FBA[MN.biomass]:.5}")

print(f"MN Target -> {MN.target}")

print(f"MN Min prod -> {MN.minprod:.5}")

print(f"Alternative Object")

mn  =Met_Net(S=S,LB=LB,UB=UB,Met=Met,Rxn=Rxn,biomass=biomas,chemical=chemical)

print(f"MN FBA {mn.FBA[MN.biomass]:.5}")

print(f"MN Target -> {mn.target}")

print(f"MN Min prod -> {mn.minprod:.5}")


print(f"Changing the values from the minprod and target")

mn.target = .1

print(f"MN Target -> {mn.target}")

print(f"MN Min prod -> {mn.minprod:.5}")