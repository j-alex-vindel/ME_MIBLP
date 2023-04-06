import os
import sys
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P
from Ob_Met_Net_solmethods import MILP_sol_OP,CB_sol_OP, Inner_check_vs_ys_NOP
from Ob_Met_Net import Met_Net
import pandas as pd

metnet = Met_Net

target = 1.0

points = {'MB':[],
          'MC':[],
          'CB':[],
          'CC':[],
          'PB':[],
          'PC':[],
          'tgt':[],
          'Bio':[],
          'Che':[],
          'Mys':[],
          'Cys':[],
          'Pys':[],
          'ICB_mb':[],
          'ICC_mb':[],
          'ICB_mc':[],
          'ICC_mc':[],

          'ICB_cb':[],
          'ICC_cb':[],
          'ICB_cc':[],
          'ICC_cc':[],
          
          'ICB_pb':[],
          'ICC_pb':[],
          'ICB_pc':[],
          'ICC_pc':[]}

while target > 0:
    metnet.target = target
    print(f"Minprod: {metnet.minprod}")
    m = MILP_sol_OP(network=metnet,k=2)
    c = CB_sol_OP(network=metnet,k=2)
    p = CB_P(network=metnet,k=2)
    m_bio = m.Vs[metnet.biomass]
    m_che = m.Vs[metnet.chemical]
    c_bio = c.Vs[metnet.biomass]
    p_bio = p.Vs[metnet.biomass]
    p_che = p.Vs[metnet.chemical]
    c_che = c.Vs[metnet.chemical]
    bio = metnet.FVA[metnet.biomass]
    che = metnet.FVA[metnet.chemical]
    
    icb_m = Inner_check_vs_ys_NOP(network=metnet,result_milp=m,criteria='ys',objective='biomass')
    icc_m = Inner_check_vs_ys_NOP(network=metnet,result_milp=m,criteria='ys',objective='chemical')
    
    icb_c = Inner_check_vs_ys_NOP(network=metnet,result_cb=c,criteria='ys',objective='biomass')
    icc_c = Inner_check_vs_ys_NOP(network=metnet,result_cb=c,criteria='ys',objective='chemical')

    icb_p = Inner_check_vs_ys_NOP(network=metnet,result_cb=p,criteria='ys',objective='biomass')
    icc_p = Inner_check_vs_ys_NOP(network=metnet,result_cb=p,criteria='ys',objective='chemical')


    points['MB'].append(m_bio)
    points['MC'].append(m_che)
    points['CB'].append(c_bio)
    points['CC'].append(c_che)
    points['Bio'].append(bio)
    points['Che'].append(che)
    points['tgt'].append(f"{target:.2}")
    points['PB'].append(p_bio)
    points['PC'].append(p_che)
    points['Mys'].append(m.Strategy)
    points['Cys'].append(c.Strategy)
    points['Pys'].append(p.Strategy)

    points['ICB_mb'].append(icb_m.Biomass)
    points['ICC_mb'].append(icc_m.Biomass)
    points['ICB_mc'].append(icb_m.Chemical)
    points['ICC_mc'].append(icc_m.Chemical)

    points['ICB_cb'].append(icb_c.Biomass)
    points['ICC_cb'].append(icc_c.Biomass)
    points['ICB_cc'].append(icb_c.Chemical)
    points['ICC_cc'].append(icc_c.Chemical)

    points['ICB_pb'].append(icb_p.Biomass)
    points['ICC_pb'].append(icc_p.Biomass)
    points['ICB_pc'].append(icb_p.Chemical)
    points['ICC_pc'].append(icc_p.Chemical)

    target -= .1


rdf = pd.DataFrame.from_dict(points)
rdf.round(decimals=5)


file_name_r = f"../Results/Envelopes/Full_Detail_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:   
    rdf.to_csv(file_name_r)



