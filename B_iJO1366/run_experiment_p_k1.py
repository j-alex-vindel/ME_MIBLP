import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net_solmethods import CB_sol_2_P
from Ob_Met_Net import Metabolic_Network

import pandas as pd

metnet = Metabolic_Network
K = sys.argv[0][-4]

cb = CB_sol_2_P(network=metnet,k=K,log=True)

