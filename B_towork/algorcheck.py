import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P
from Ob_Met_Net_solmethods import Inner_check_vs_ys_NOP
from Algorithms import BilevelMethods
from strainselector import strainsele,strain_id,method_id
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from collections import namedtuple
from typing import List,Type
import math


Result = namedtuple('Result_cb',['MetNet','Strategy','Ys','Vs','Vij','Time','Soltype','Method'])
Tgt = List[float]
Knock = List[int]
MN = Type[object]

def Pess_Algo():
    new_r = {'C_Bio':[],
            'C_Che':[],
            'Method': [],
            'Time': [],
            'Tgt': [],
            'Strain':[],
            'K': [],
            'Strategy':[],
            'KIndex':[],
            'E_Bio':[],
            'E_Che':[],
            'ICB_Bio':[],
            'ICB_Che':[],
            'ICC_Bio':[],
            'ICC_Che':[],
            'ID':[],
            'Pct_Che':[],
            'Pct_Bio':[],
            'Soltype':[]}


    for index,pair in enumerate(product(tgts,KS)):
        target,k = pair
        metnet.target = target/100
        print(f"Strain:{metnet.Name[:3]} -> MinProd: {metnet.minprod:.2} -> K: {k} -> {index+1}/{len(list(product(tgts,KS)))}")
        p = solver.CB_P(network=metnet,K=k)

        icb_p = solver.FBA_check(network=metnet,solution=p,obj_v='biomass',c_params='ys')
        icc_p = solver.FBA_check(network=metnet,solution=p,obj_v='chemical',c_params='ys')

        bio = metnet.FVA[metnet.biomass]
        che = metnet.FVA[metnet.chemical]

        p_bio = p.Vs[metnet.biomass]
        p_che = p.Vs[metnet.chemical]
        
        
        if v_che != 0:
            pct_che = ((p_che - v_che)/v_che)*100
        else:
            pct_che = math.inf

        pct_bio = (p_bio/v_bio)*100


        pid = strain_id(metnet.Name[:3]) + method_id(p.Method) + target + k
        
        # c - p - m
        new_r['Tgt'].extend([target])
        new_r['K'].extend([k])
        new_r['Strain'].extend([name])

        new_r['C_Bio'].extend([p_bio])
        new_r['C_Che'].extend([p_che])
        new_r['Method'].extend([p.Method])
        new_r['Time'].extend([p.Time])
        new_r['Strategy'].extend([p.Strategy])
        new_r['KIndex'].extend([p.KO_index])
        new_r['ICB_Bio'].extend([icb_p.Biomass])
        new_r['ICB_Che'].extend([icb_p.Chemical])
        new_r['ICC_Bio'].extend([icc_p.Biomass])
        new_r['ICC_Che'].extend([icc_p.Chemical])
        
        new_r['E_Bio'].extend([bio])
        new_r['E_Che'].extend([che])
        new_r['ID'].extend([pid])
        new_r['Pct_Bio'].extend([pct_bio])
        new_r['Pct_Che'].extend([pct_che])
        new_r['Soltype'].extend([p.Soltype])

    for i,v in enumerate(new_r['C_Che']):
        print(f"tgt {tgts[i]} -> {new_r['C_Che'][i]:.3} -> {new_r['ICB_Che'][i]:.3} -> {new_r['Strategy'][i]} -> {new_r['Soltype'][i]}")
    
    
    if write in ['y',"Y"]:
        df = pd.DataFrame.from_dict(new_r)
        df.to_csv(f"../Results/Envelopes/Revised/DR/DEF_P_{name}.csv")
        sys.exit('File created - Run terminated!')
    else:
        sys.exit('Program Terminated - No csv File!')

def Pess_D(write:str=None,tgts:Tgt=None,KS:Knock=None,mn:MN=None):
   
    new_r = {'C_Bio':[],
            'C_Che':[],
            'Method': [],
            'Time': [],
            'Tgt': [],
            'Strain':[],
            'K': [],
            'Strategy':[],
            'KIndex':[],
            'E_Bio':[],
            'E_Che':[],
            'ICB_Bio':[],
            'ICB_Che':[],
            'ICC_Bio':[],
            'ICC_Che':[],
            'ID':[],
            'Pct_Che':[],
            'Pct_Bio':[]}
    for pair in product(tgts,KS):
        target,k = pair
        mn.target = target/100

        p = CB_P(network=mn,k=k,log=True)
        ibp = Inner_check_vs_ys_NOP(network=mn,result_cb=p,criteria='ys',objective='biomass')
        icp = Inner_check_vs_ys_NOP(network=mn,result_cb=p,criteria='ys',objective='chemical')
        
        bio = mn.FVA[mn.biomass]
        che = mn.FVA[mn.chemical]

        p_bio = p.Vs[mn.biomass]
        p_che = p.Vs[mn.chemical]
        
        
        if v_che != 0:
            pct_che = ((p_che - v_che)/v_che)*100
        else:
            pct_che = math.inf

        pct_bio = (p_bio/v_bio)*100


        pid = strain_id(mn.Name[:3]) + method_id(p.Method) + target + k
        
        # c - p - m
        new_r['Tgt'].extend([target])
        new_r['K'].extend([k])
        new_r['Strain'].extend([name])

        new_r['C_Bio'].extend([p_bio])
        new_r['C_Che'].extend([p_che])
        new_r['Method'].extend([p.Method])
        new_r['Time'].extend([p.Time])
        new_r['Strategy'].extend([p.Strategy])
        new_r['KIndex'].extend([p.KO_index])
        new_r['ICB_Bio'].extend([ibp.Biomass])
        new_r['ICB_Che'].extend([ibp.Chemical])
        new_r['ICC_Bio'].extend([icp.Biomass])
        new_r['ICC_Che'].extend([icp.Chemical])
        
        new_r['E_Bio'].extend([bio])
        new_r['E_Che'].extend([che])
        new_r['ID'].extend([pid])
        new_r['Pct_Bio'].extend([pct_bio])
        new_r['Pct_Che'].extend([pct_che])
        # new_r['Soltype'].extend([p.Soltype])
    
    if write in ['y',"Y"]:
        df = pd.DataFrame.from_dict(new_r)
        df.to_csv(f"../Results/Envelopes/Revised/DR/DEF_P3_{name}.csv")
        sys.exit('File created - Run terminated!')
    elif write in ['e','E']:
            sys.exit('Program Terminated - No csv File!')
    else:
        for i,v in enumerate(new_r['C_Che']):
            print(f"tgt {tgts[i]} -> {new_r['C_Che'][i]:.3} -> {new_r['ICB_Che'][i]:.3} -> {new_r['Strategy'][i]} -> {new_r['Soltype'][i]}")
    

def yvector(strain=None,index=None):
    yv = [0 if i in index else 1 for i in strain.M]
    return yv

strain = sys.argv[1]

write = sys.argv[2]

debug = sys.argv[3]
# index = [int(i) for i in sys.argv[3:]]

metnet = strainsele(strain)
name = metnet.Name[:3]

v_che = metnet.FBA[metnet.chemical]
v_bio = metnet.FBA[metnet.biomass]

print(name)
print(v_che)
print(v_bio)
print(metnet.Rxn[metnet.chemical])



solver = BilevelMethods(log=False)

obj = metnet.Rxn[metnet.chemical]

if sys.argv[4] in ['a','all','A','All']:
    KS = [1,2,3]
else:
    KS = [int(sys.argv[4])]

if sys.argv[5] in ['All','A','a','all']:
    tgts = [_ for _ in range(10,100,10)]
elif sys.argv[5] in ['-']:
    tgts = None
else:
    tgts = [int(i) for i in sys.argv[5:]]


if debug not in ['D','d']:
    Pess_Algo()

elif debug in ['D','d']:
    Pess_D(write=write,tgts=tgts,KS=KS,mn=metnet)


else:
    sys.exit()

