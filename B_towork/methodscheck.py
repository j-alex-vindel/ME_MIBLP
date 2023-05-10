import os
import sys
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iAF1260'))) 
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iJR904'))) 
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iJO1366'))) 

import Algorithms
from Ob_Met_Net_solmethods import CB_sol_OP,Inner_check_vs_ys_NOP
from CB_Sol_P import CB_P
from MN_iaf1260 import MN_iaf1260
from MN_ijr904 import MN_ijr904
from MN_ijo1366 import MN_ijo1366
import pandas as pd

def strainselector(strain:str=None):
    if strain == 'ijo':
        met = MN_ijo1366
    elif strain == 'ijr':
        met = MN_ijr904
    elif strain == 'iaf':
        met = MN_iaf1260
    return met

def checkalgori(df:pd.DataFrame=None):
    for i in range(len(df)):
        row = df.iloc[i]
        yield row.Tgt, row.K, row.Method

def dataframeselector(df:pd.DataFrame=None,check:str=None):
    if check == 'A':
        t_df = df[(df.Strain == metnet.Name)&(df.IC_Che ==1)&(df.IB_Bio == 1)]
    elif check == 'B':
        t_df = df[(df.Strain == metnet.Name)&(df.IB_Bio == 1)]
    elif check == 'C':
        t_df = df[(df.Strain == metnet.Name)&(df.IC_Che == 1)]
    return t_df



if len(sys.argv) <3:
    sys.exit('Insuficient arguments')
else:
    strain  = sys.argv[1]
    if sys.argv[2] not in ['A','B','C']:
        sys.exit('Write A or B')
    else:
        check = sys.argv[2]

metnet = strainselector(strain=strain)

df = pd.read_csv("../Rev_iAF_iJO_iJR_m.csv")

tdf = dataframeselector(df=df,check=check)

print(f"Name: {metnet.Name}")
print(f"Strain: {strain}")
print(f"Check: {check}")

if tdf.empty:
    sys.exit('Empty Data Frame')

else:
    solver = Algorithms.BilevelMethods(log=False)

    results = {'C_Bio':[],
          'C_Che':[],
          'Tgt':[],
          'K':[],
          'Method':[],
          'ICB_Bio':[],
          'ICC_Che':[],
          'ICC_Bio':[],
          'ICC_Che':[]}
    for pair in checkalgori(tdf):
        target,k,method = pair
        if method == 'O':
            print(f"optimistic")
            print(f"{target}->{k}")
            metnet.target = target/100
            print(f"FBA [b]:{metnet.FBA[metnet.biomass]:.2}")
            print(f"FVA [b]:{metnet.FVA[metnet.biomass]:.2}")
            optimistic = CB_sol_OP(network=metnet,k=k)
            ocb = Inner_check_vs_ys_NOP(network=metnet,result_cb=optimistic,criteria='ys',objective='biomass')
            occ = Inner_check_vs_ys_NOP(network=metnet,result_cb=optimistic,criteria='ys',objective='chemical')

            results['C_Bio'].append(optimistic.Vs[metnet.biomass])
            results['C_Che'].append(optimistic.Vs[metnet.chemical])
            results['Tgt'].append(target)
            results['K'].append(k)
            results['Method'].append(method)
            results['ICB_Bio'].append(ocb.Biomass)
            results['ICB_Che'].append(ocb.Chemical)
            results['ICC_Bio'].append(occ.Biomass)
            results['ICC_Che'].append(occ.Chemical)

        elif method == 'P':
            print(f"pessimistic")
            print(f"{target}->{k}")
            metnet.target = target/100
            print(f"FBA [b]:{metnet.FBA[metnet.biomass]:.2}")
            print(f"FVA [b]:{metnet.FVA[metnet.biomass]:.2}")
            pessimistic = CB_sol_OP(network=metnet,k=k)
            pcb = Inner_check_vs_ys_NOP(network=metnet,result_cb=pessimistic,criteria='ys',objective='biomass')
            pcc = Inner_check_vs_ys_NOP(network=metnet,result_cb=pessimistic,criteria='ys',objective='chemical')
            
            results['C_Bio'].append(pessimistic.Vs[metnet.biomass])
            results['C_Che'].append(pessimistic.Vs[metnet.chemical])
            results['Tgt'].append(target)
            results['K'].append(k)
            results['Method'].append(method)
            results['ICB_Bio'].append(pcb.Biomass)
            results['ICB_Che'].append(pcb.Chemical)
            results['ICC_Bio'].append(pcc.Biomass)
            results['ICC_Che'].append(pcc.Chemical)

    cdf = pd.DataFrame.from_dict(results)

    filename = f"../Results/Envelopes/Revised/DR/Revised_{metnet.Name[:3]}_C{check}.csv"

    cdf.to_csv(filename)

        



    


