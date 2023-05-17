import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Algorithms import BilevelMethods
from strainselector import strainsele

strain = 'ijo'
metnet = strainsele(strain=strain)

solver = BilevelMethods(log=False)

target = .5
k =2

m = solver.MILP(network=metnet,K=k)

print(f"Biomass {m.Vs[metnet.biomass]}")
print(f"Strat: {m.Strategy}")
print(f"Chemical {m.Vs[metnet.chemical]}")