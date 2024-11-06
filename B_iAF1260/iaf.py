import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES')))
from support_functions import Wildtpe_Check
from Ob_Met_Net import Met_Net
from CB_Sol_P import CB_P
from Ob_Met_Net_solmethods import CB_sol_OP, MILP_sol_OP
from pymatreader import read_mat
import pandas as pd

bacteria = 'iAF1260'

data = read_mat(f"../B_iAF1260/Data/iAF1260.mat")

LB = data[bacteria]['lb'].tolist()
UB = data[bacteria]['ub'].tolist()
Met = data[bacteria]['mets']
Rxn = data[bacteria]['rxns']
S = data[bacteria]['S']

LB[Rxn.index('EX_glc__D_e')] = -10
UB[Rxn.index('EX_glc__D_e')] = -10

chemical = Rxn.index('EX_succ_e')
biomas = Rxn.index('BIOMASS_Ec_iAF1260_core_59p81M')

non_essentials = ['GLCabcpp', 'GLCptspp', 'HEX1', 'PGI', 'PFK', 'FBA', 'TPI', 'GAPD','PGK', 'PGM', 'ENO', 'PYK',
'LDH_D', 'PFL', 'ALCD2x', 'PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH', 'ACONTa', 'ACONTb', 'ICDHyr', 'CS', 'MDH','MDH2', 'MDH3', 'ACALD']

knockout = [Rxn.index(i) for i in non_essentials]



# =========== Metabolic Network Object ====================

mn = Met_Net(S=S,LB=LB,UB=UB,Rxn=Rxn,KO=knockout,Met=Met,biomass=biomas,chemical=chemical,Name=bacteria)



# -------------------- FROM RObustKnock ----------------------------

Ys = [1 for _ in mn.M]

WT = Wildtpe_Check(obj=mn,wildtype=True,Ys=Ys)

# ------------------ FROM CB Pess ---------------------


tgt = int(sys.argv[2])/100
k = int(sys.argv[1])
mn.target = tgt 

p = CB_P(network=mn,k=k,log=False)
o = CB_sol_OP(network=mn,k=k,log=False)
m = MILP_sol_OP(network=mn,k=k,log=False)


# print(f"\nPessimistic")
# print(f"Strategy {p.Strategy}")
# for i in p.Strategy:
#     print(data[bacteria]['rxnNames'][Rxn.index(i)])

# print(f"v[b] {p.Vs[mn.biomass]} -> %{(p.Vs[mn.biomass]/WT[mn.biomass])*100:.5}")
# print(f"v[c] {p.Vs[mn.chemical]}")
# print(f"Time: {p.Time}\n")

# print(f"Optimistic")
# print(f"Strategy {o.Strategy}")
# for i in o.Strategy:
#     print(data[bacteria]['rxnNames'][Rxn.index(i)])
# print(f"v[b] {o.Vs[mn.biomass]} -> %{(o.Vs[mn.biomass]/WT[mn.biomass])*100:.5}")
# print(f"v[c] {o.Vs[mn.chemical]}")
# print(f"Time: {o.Time}\n")

# print(f"MILP")
# print(f"Strategy {m.Strategy}")
# for i in m.Strategy:
#     print(data[bacteria]['rxnNames'][Rxn.index(i)])
# print(f"v[b] {m.Vs[mn.biomass]} -> %{(m.Vs[mn.biomass]/WT[mn.biomass])*100:.5}")
# print(f"v[c] {m.Vs[mn.chemical]}")
# print(f"Time: {m.Time}\n")


r = {'v[b]':[m.Vs[mn.biomass],o.Vs[mn.biomass],p.Vs[mn.biomass]],
     'v[c]':[m.Vs[mn.chemical],o.Vs[mn.chemical],p.Vs[mn.chemical]],
     'Bobj':[mn.Rxn[mn.chemical],mn.Rxn[mn.chemical],mn.Rxn[mn.chemical]],
     'pctBio':[m.Vs[mn.biomass]/WT[mn.biomass],o.Vs[mn.biomass]/WT[mn.biomass],p.Vs[mn.biomass]/WT[mn.biomass]],
     'K':[k,k,k],
    'Strat':[m.Strategy,o.Strategy,p.Strategy],
    "Method":['M','O','P'],
    'Time':[m.Time,o.Time,p.Time],
    "minprod":[tgt,tgt,tgt],
    'Bacteria':[mn.Name[:3],mn.Name[:3],mn.Name[:3]]}

df = pd.DataFrame.from_dict(r)


os.makedirs('Results/OPM2',exist_ok=True)

df.to_csv(f"Results/OPM2/OPM_{bacteria[:3]}.csv",mode='a')