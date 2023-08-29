import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net_solmethods import MILP_sol_OP,CB_sol_OP
from strainselector import strainsele
# from MN_MOMO import MN_MOMO_ijr

import pandas as pd

metnet = strainsele('momo')
metnet.target = .10
print(f"Biomass: {metnet.FBA[metnet.biomass]}")
print(f"Minprod {metnet.minprod}")

K = 2

# YS = [0 if i in [2,6] else 1 for i in metnet.M]


m2 = MILP_sol_OP(network=metnet,k=K,log=False)

print(m2.Strategy)
print(m2.Vs[metnet.biomass],"->",metnet.Rxn[metnet.biomass])
print(m2.Vs[metnet.chemical],"->",metnet.Rxn[metnet.chemical])

c2 = CB_sol_OP(network=metnet,k=K,log=False)

print(c2.Strategy)
print(c2.Vs[metnet.biomass],"->",metnet.Rxn[metnet.biomass])
print(c2.Vs[metnet.chemical],"->",metnet.Rxn[metnet.chemical])