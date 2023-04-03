import os
import sys
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from MN_ijr904 import MN_ijr904
import matplotlib.pyplot as plt
import pandas as pd

metnet = MN_ijr904

target = 1
points = []
while target > 0:
    print(target)
    metnet.target = target
    b = metnet.FVA[metnet.biomass]
    c = metnet.FVA[metnet.chemical]
    points.append((b,c))
    target -= .1

# Data Normalization

bp = [point[0] for point in points]
cp = [point[1] for point in points]

max_b = max(point[0] for point in points)
max_c = max(point[1] for point in points)

normb = [point[0]/max_b for point in points]
normc = [point[1]/max_c for point in points]


r = {'B':bp,
     'C':cp,
     'N_B':normb,
     'N_C':normc}

file_name_r = f"../Results/Envelopes/FE_{metnet.Name}.csv"

isfile_r = os.path.exists(file_name_r)

if not isfile_r:
    rdf = pd.DataFrame.from_dict(r)
    rdf.to_csv(file_name_r)


fig,ax = plt.subplots()
plt.plot(bp,cp)

ax.set(xlim=(0,max(bp)+.01),ylim=(0,max(cp)+1))
ax.set_title(f"Strain ({metnet.Name}) Production Envelope")
plt.xlabel(r"""Biomass production $mmol/g(Dw)h$""")
plt.ylabel(r"""Chemical production $mmol/g(Dw)h$""")

plt.show()
plt.savefig(f"../Results/Graphs/NON_FE_{metnet.Name}.png")

fig,ax = plt.subplots()
plt.plot(normb,normc)

ax.set(xlim=(0,1.01),ylim=(0,1.01))
ax.set_title(f"Strain ({metnet.Name}) Normalized Production Envelope")
plt.xlabel('% Biomass')
plt.ylabel('% Chemical')
plt.show()
plt.savefig(f"../Results/Graphs/N_FE_{metnet.Name}.png")