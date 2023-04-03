import os
import sys
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 


from Ob_Met_Net_solmethods import MILP_sol_OP,CB_sol_OP
from Ob_Met_Net import Met_Net
import pandas as pd

metnet = Met_Net

target = 1.0

points = {'MB':[],
          'MC':[],
          'CB':[],
          'CC':[],
          'tgt':[],
          'Bio':[],
          'Che':[]}

while target > 0:
    metnet.target = target
    print(f"Minprod: {metnet.minprod}")
    m = MILP_sol_OP(network=metnet,k=2)
    c = CB_sol_OP(network=metnet,k=2)
    m_bio = m.Vs[metnet.biomass]
    m_che = m.Vs[metnet.chemical]
    c_bio = c.Vs[metnet.biomass]
    c_che = c.Vs[metnet.chemical]
    bio = metnet.FVA[metnet.biomass]
    che = metnet.FVA[metnet.chemical]
    points['MB'].append(m_bio)
    points['MC'].append(m_che)
    points['CB'].append(c_bio)
    points['CC'].append(c_che)
    points['Bio'].append(bio)
    points['Che'].append(che)
    points['tgt'].append(f"{target:.2}")
    target -= .1


rdf = pd.DataFrame.from_dict(points)
rdf.round(decimals=5)


file_name_r = f"../Results/Envelopes/Full_Detail_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:   
    rdf.to_csv(file_name_r)



