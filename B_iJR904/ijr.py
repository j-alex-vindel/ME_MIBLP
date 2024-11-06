import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES')))
from support_functions import Wildtpe_Check
from Ob_Met_Net import Met_Net
from CB_Sol_P import CB_P
from Ob_Met_Net_solmethods import CB_sol_OP, MILP_sol_OP
from pymatreader import read_mat
import pandas as pd

ijr904 = 'iJR904'

data = read_mat(f"../B_iJR904/Data/iJR904.mat")

LB = data[ijr904]['lb'].tolist()
UB = data[ijr904]['ub'].tolist()
Met = data[ijr904]['mets']
Rxn = data[ijr904]['rxns']
S = data[ijr904]['S']

LB[Rxn.index('EX_glc__D_e')] = -10
UB[Rxn.index('EX_glc__D_e')] = -10


LB[Rxn.index("EX_o2_e")] = -20

LB[Rxn.index("ATPM")] = 7.6
UB[Rxn.index("ATPM")] = 7.6

# print(data[ijr904]['rxnNames'].index('ATP synthase (four protons for one ATP)'))
# 1. ATP 222

# for index,rxn in enumerate(data[ijr904]['rxns']):
#     if rxn[:3] in ("ATP"):
#         print(index,'->',rxn,'->',data[ijr904]['rxnNames'][index])
        
# # 2. 784
# print(data[ijr904]['rxnNames'].index('Phosphoglycerate mutase'))
# print(data[ijr904]['rxns'][784])


# print(data[ijr904]['rxnNames'].index("Succinate dehydrogenase"))
# print(data[ijr904]['rxns'][927])

k_pes = [222,784,927]
# print("KS")
# for k in ks:
#     print(data[ijr904]['rxns'][k],'->',data[ijr904]['rxnNames'][k])

# print(" ")
# for index,rxn in enumerate(data[ijr904]['rxnNames']):
#     if rxn[:3] in ("Ace"):
#         print(index,'->',data[ijr904]['rxns'][index],'->',rxn)

# To find the index of the optimistic approach
# for index,rxn in enumerate(data[ijr904]['rxnNames']):
#     if rxn in ['2-dehydro-3-deoxy-phosphogluconate aldolase','Triose-phosphate isomerase']:
#         print(index,'->',rxn)

k_opt = [222,302,1033]


biomass = Rxn.index('BIOMASS_Ecoli')
chemical = Rxn.index('EX_ac_e') #acetate


mn = Met_Net(S=S,LB=LB,UB=UB,Rxn=Rxn,Met=Met,biomass=biomass,chemical=chemical,Name=ijr904)



# -------------------- FROM RObustKnock ----------------------------
y_pes = [0 if i in k_pes else 1 for i in mn.M]
Ys = [1 for _ in mn.M]

WT = Wildtpe_Check(obj=mn,wildtype=True,Ys=Ys)

# FBA = Wildtpe_Check(obj=mn,wildtype=True,Ys=y_pes)


# FVA = Wildtpe_Check(obj=mn,wildtype=False,mutant=True,Ys=y_pes)



strat_pes = [mn.Rxn[i] for i in k_pes]


# ------------------ FROM CB Pess ---------------------

non_essentials = ['PFL','PTAr', 'ACKr','G6PDH2r', 'PGL', 'GND', 'RPI', 'RPE', 'TKT1', 'TALA', 'TKT2', 'FUM',
'FRD2', 'SUCOAS', 'AKGDH']


ki = [Rxn.index(i) for i in non_essentials]

kn = k_pes.copy()
kn.extend(ki)
kn.extend(k_opt)

ko = list(set(kn))
tgt = int(sys.argv[2])/100
k = int(sys.argv[1])

mn.KO = ko
mn.target = tgt

p = CB_P(network=mn,k=k,log=False)
o = CB_sol_OP(network=mn,k=k,log=False)
m = MILP_sol_OP(network=mn,k=k,log=False)


# print(f"\nPessimistic")
# print(f"Strategy {p.Strategy}")
# for i in p.Strategy:
#     print(data[ijr904]['rxnNames'][Rxn.index(i)])

# print(f"v[b] {p.Vs[mn.biomass]} -> %{(p.Vs[mn.biomass]/WT[mn.biomass])*100:.5}")
# print(f"v[c] {p.Vs[mn.chemical]}")
# print(f"Time: {p.Time}\n")

# print(f"Optimistic")
# print(f"Strategy {o.Strategy}")
# for i in o.Strategy:
#     print(data[ijr904]['rxnNames'][Rxn.index(i)])
# print(f"v[b] {o.Vs[mn.biomass]} -> %{(o.Vs[mn.biomass]/WT[mn.biomass])*100:.5}")
# print(f"v[c] {o.Vs[mn.chemical]}")
# print(f"Time: {o.Time}\n")

# print(f"MILP")
# print(f"Strategy {m.Strategy}")
# for i in m.Strategy:
#     print(data[ijr904]['rxnNames'][Rxn.index(i)])
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

df.to_csv(f"Results/OPM2/OPM_{ijr904[:3]}.csv",mode='a')

