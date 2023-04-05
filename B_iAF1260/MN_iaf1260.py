 
import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net import Met_Net
from pymatreader import read_mat

iaf1260 = 'iAF1260'

# =========== Access Data ====================
data_iaf1260 = read_mat(f"../B_iAF1260/Data/iAF1260.mat")

LB  = data_iaf1260[iaf1260]['lb'].tolist()
met = data_iaf1260[iaf1260]['mets']
UB  = data_iaf1260[iaf1260]['ub'].tolist()
rxn = data_iaf1260[iaf1260]['rxns']
S   = data_iaf1260[iaf1260]['S']

# =========== Biological Assumptions ====================
# from  iAF1260_FBA_Mu_Min.ipynb
# Change Biological Assumptions

# print("Before Changing Biological Params")

# print(f"->>Glucose LB: {LB[rxn.index('EX_glc__D_e')]:.5}; UB:{UB[rxn.index('EX_glc__D_e')]:.5}")
# print(f"->>Oxigen LB: {LB[rxn.index('EX_o2_e')]:.5}; UB:{UB[rxn.index('EX_o2_e')]:.5}")
# print(f"->>ATPM LB: {LB[rxn.index('ATPM')]:.5}; UB:{UB[rxn.index('ATPM')]:.5}")


LB[rxn.index('EX_glc__D_e')] = -10
UB[rxn.index('EX_glc__D_e')] = -10

chemical = rxn.index('EX_succ_e')
biomas = rxn.index('BIOMASS_Ec_iAF1260_core_59p81M')

non_essentials = ['GLCabcpp', 'GLCptspp', 'HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK',
'LDH_D', 'PFL', 'ALCD2x', 'PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH', 'ACONTa', 'ACONTb', 'ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

knockout = [rxn.index(i) for i in non_essentials]


# =========== Metabolic Network Object ====================

MN_iaf1260 = Met_Net(S=S,LB=LB,UB=UB,Rxn=rxn,Met=met,Name=iaf1260,KO=knockout,biomass= biomas, chemical=chemical)

