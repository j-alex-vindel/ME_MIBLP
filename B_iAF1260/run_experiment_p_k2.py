import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net_solmethods import CB_sol_2_P
from Ob_Met_Net_iaf1260 import MetNet_iaf1260

import pandas as pd

metnet = MetNet_iaf1260
K = int(sys.argv[0][-4])

cb = CB_sol_2_P(network=metnet,k=K,log=True)

print(f"Solutions")

print(f"VS biomas: {cb.Vs[metnet.biomass]}")
print(f"VS chemical: {cb.Vs[metnet.chemical]}")
print(f"Strategy:{cb.Strategy}")