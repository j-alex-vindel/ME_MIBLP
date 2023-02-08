import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P
from Ob_Met_Net_solmethods import CB_sol_2_NOP
from Ob_Met_Net_ijr904 import MetNet_ijr904

import pandas as pd

metnet = MetNet_ijr904
K = 1

cb1_p = CB_P(network=metnet,k=K,log=True)

cb1_o = CB_sol_2_NOP(network=metnet,k=K,log=True)


cb2_p = CB_P(network=metnet,k=K+1,log=True)
cb2_o = CB_sol_2_NOP(network=metnet,k=K+1,log=True)


cb3_p = CB_P(network=metnet,k=K+2,log=True)
cb3_o = CB_sol_2_NOP(network=metnet,k=K+2,log=True)


r  = {
            'Approach':['O'                     ,'P'                      ,'O'                      ,'P'                      ,'O'                      ,'P'],
            'Biomass': [cb1_o.Vs[metnet.biomass], cb1_p.Vs[metnet.biomass], cb2_o.Vs[metnet.biomass], cb2_p.Vs[metnet.biomass], cb3_o.Vs[metnet.biomass],cb3_p.Vs[metnet.biomass]],
            'Chemical':[cb1_o.Vs[metnet.chemical], cb1_p.Vs[metnet.chemical], cb2_o.Vs[metnet.chemical], cb2_p.Vs[metnet.chemical], cb3_o.Vs[metnet.chemical],cb3_p.Vs[metnet.chemical]],
            'Time':    [cb1_o.Time, cb1_p.Time, cb2_o.Time, cb2_p.Time, cb3_o.Time,cb3_p.Time],
            'Strategy':[cb1_o.Strategy, cb1_p.Strategy, cb2_o.Strategy, cb2_p.Strategy, cb3_o.Strategy,cb3_p.Strategy]
}

file_name_r = f"../Results/Methods_OP_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:
    rdf = pd.DataFrame.from_dict(r)
    rdf.to_csv(file_name_r)