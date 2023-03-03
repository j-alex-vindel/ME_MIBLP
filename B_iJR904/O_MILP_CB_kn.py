import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net_solmethods import MILP_sol_OP,CB_sol_OP
from support_functions import save_results, do_them_graphs
from MN_ijr904 import MN_ijr904

import pandas as pd

metnet = MN_ijr904
metnet.target = .1
K = 1

c1 = CB_sol_OP(network=metnet,k=K,log=False)
m1 = MILP_sol_OP(network=metnet,k=K,log=False)


c2 = CB_sol_OP(network=metnet,k=K+1,log=False)
m2 = MILP_sol_OP(network=metnet,k=K+1,log=False)


c3 = CB_sol_OP(network=metnet,k=K+2,log=False)
m3 = MILP_sol_OP(network=metnet,k=K+2,log=False)


results = [c1,c2,c3,m1,m2,m3]

for result in results:
    save_results(result)

# --------------- Compile Results ----------------------------------------------------------
r  = {
            'Method':[m1.Method,c1.Method,m2.Method,c2.Method,m3.Method,c3.Method],
            'Biomass':[m1.Vs[metnet.biomass],c1.Vs[metnet.biomass],m2.Vs[metnet.biomass],c2.Vs[metnet.biomass],m3.Vs[metnet.biomass],c3.Vs[metnet.biomass]],
            'Chemical':[m1.Vs[metnet.chemical],c1.Vs[metnet.chemical],m2.Vs[metnet.chemical],c2.Vs[metnet.chemical],m3.Vs[metnet.chemical],c3.Vs[metnet.chemical]],
            'Time':[m1.Time,c1.Time,m2.Time,c2.Time,m3.Time,c3.Time],
            'Strategy':[m1.Strategy,c1.Strategy,m2.Strategy,c2.Strategy,m3.Strategy,c3.Strategy],
            'K':[len(m1.Strategy),len(c1.Strategy),len(m2.Strategy),len(c2.Strategy),len(m3.Strategy),len(c3.Strategy)]
}

rdf = pd.DataFrame.from_dict(r)
rdf.round(decimals=5)


file_name_r = f"../Results/Optimistic/Op_K{sys.argv[0][-4]}_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:   
    rdf.to_csv(file_name_r)

# -----------------  Graphs ---------------------------------------------------------------- 
do_them_graphs(metnet,rdf,'Opt')




        