import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P
from Ob_Met_Net_ijo1366 import MetNet_ijo1366

import pandas as pd

metnet = MetNet_ijo1366
K = int(sys.argv[0][-4])

cb = CB_P(network=metnet,k=K,log=True)

print(f"Solutions")

print(f"VS biomas: {cb.Vs[metnet.biomass]}")
print(f"VS chemical: {cb.Vs[metnet.chemical]}")
print(f"Strategy:{cb.Strategy}")
