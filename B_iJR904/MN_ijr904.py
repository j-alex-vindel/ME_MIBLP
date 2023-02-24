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

# Identiying the index of the reactions, biomass and chemical of interest 
biomas = Rxn.index('BIOMASS_Ecoli')
chemical = Rxn.index('EX_succ_e')

# Biological Assumptions

LB[Rxn.index('EX_glc__D_e')] = -10
UB[Rxn.index('EX_glc__D_e')] = -10

UB[Rxn.index("EX_o2_e")] = -20
LB[Rxn.index("EX_o2_e")] = -20

LB[Rxn.index("ATPM")] = 7.6
UB[Rxn.index("ATPM")] = 7.6

# Defining non essentials and KO reactions

non_essentials = ['HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK',
'LDH_D', 'PFL','PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH','ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

knockout = [Rxn.index(i) for i in non_essentials]


# Creating Object 

MN_ijr904 = Met_Net(S=S,LB=LB,UB=UB,Met=Met,Rxn=Rxn,biomass=biomas,chemical=chemical,Name=ijr904,KO=knockout)

