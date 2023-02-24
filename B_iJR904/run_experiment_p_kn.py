import sys
import os

sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P
from MN_ijr904 import MN_ijr904

import pandas as pd

metnet = MN_ijr904
metnet.target = .1

K = 1

cb1 = CB_P(network=metnet,k=K,log=True)

print(f"Solutions")

print(f"VS biomas: {cb1.Vs[metnet.biomass]}")
print(f"VS chemical: {cb1.Vs[metnet.chemical]}")
print(f"Strategy:{cb1.Strategy}")

cb2 = CB_P(network=metnet,k=K+1,log=True)

print(f"Solutions")

print(f"VS biomas: {cb2.Vs[metnet.biomass]}")
print(f"VS chemical: {cb2.Vs[metnet.chemical]}")
print(f"Strategy:{cb2.Strategy}")

cb3 = CB_P(network=metnet,k=K+2,log=True)

print(f"Solutions")

print(f"VS biomas: {cb3.Vs[metnet.biomass]}")
print(f"VS chemical: {cb3.Vs[metnet.chemical]}")
print(f"Strategy:{cb3.Strategy}")

r  = {
            'K':[1,2,3],
            'Biomass':[cb1.Vs[metnet.biomass],cb2.Vs[metnet.biomass],cb3.Vs[metnet.biomass]],
            'Chemical':[cb1.Vs[metnet.chemical],cb2.Vs[metnet.chemical],cb3.Vs[metnet.chemical]],
            'Time':[cb1.Time,cb2.Time,cb3.Time],
            'Strategy':[cb1.Strategy,cb2.Strategy,cb3.Strategy]
}

file_name_r = f"../Results/Pessimistic/Methods_K{sys.argv[0][-4]}_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:
    rdf = pd.DataFrame.from_dict(r)
    rdf.to_csv(file_name_r)