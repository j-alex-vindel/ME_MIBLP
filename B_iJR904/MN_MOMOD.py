import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES')))

from Ob_Met_Net import Met_Net
from pymatreader import read_mat
from support_functions import set_constructor,wildtype_FBA
import pandas as pd

ecoli = "../B_iJR904/Data"

S = pd.read_csv(f'{ecoli}/MOMO_SSM.txt',sep=" ",header=None).values
rxn = pd.read_csv(f'{ecoli}/MOMO_Rxn.csv',header=None)[0].tolist()
met = pd.read_csv(f'{ecoli}/MOMO_Met.csv',header=None)[0].tolist()
LB = pd.read_csv(f'{ecoli}/MOMO_SLB.txt',header=None)[0].tolist()
UB = pd.read_csv(f'{ecoli}/MOMO_SUB.txt',header=None)[0].tolist()

biomass = 12

chemical = 23


# Set of Reactions essentials from MOMO (indexes)
T_E = [2,6,9,11,12,13,14,16,17,19,42,43,45,47,50,52,56,58,63,67,69,70,78,84,87,88,92,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39]


T_p = [i-1 for i in T_E] # Reactions essentials with their proper index

knockout = [i for i in range(len(rxn)) if i not in T_p]


MN_MOMO_D = Met_Net(S=S,LB=LB,UB=UB,Rxn=rxn,Met=met,biomass=biomass,chemical=chemical,KO=knockout,Name='MoMD')

