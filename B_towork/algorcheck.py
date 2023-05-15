import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Algorithms import BilevelMethods

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'B_iJO1366')))
from MN_ijo1366 import MN_ijo1366
import pandas as pd

solver = BilevelMethods(log=False)

metnet = MN_ijo1366

tgt = float(sys.argv[1])
metnet.target = tgt
k = int(sys.argv[2])
print(f"FVA v[b] ={metnet.FVA[metnet.biomass]}")
print(f"FVA v[c] ={metnet.FVA[metnet.chemical]}")

o = solver.CB_O(network=metnet,K=k)

print(len(o.Ys))
print(f"{o.Vs[metnet.biomass]}-> {o.Vs[metnet.chemical]}")

ocb = solver.FBA_check(network=metnet,solution=o,obj_v='biomass',c_params='ys')

print(f"{ocb.Biomass} -> {ocb.Chemical}")