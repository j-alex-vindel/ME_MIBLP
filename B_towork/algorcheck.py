import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from CB_Sol_P import CB_P
from Algorithms import BilevelMethods
from strainselector import strainsele
import pandas as pd
import matplotlib.pyplot as plt


strain = sys.argv[1]

metnet = strainsele(strain)
name = metnet.Name[:3]

obj = metnet.Rxn[metnet.chemical]
K = 2

tgts = [_ for _ in range(10,100,10)]
solver= BilevelMethods(log=False)
chems_solver = []
chems_function = []
e_chems = []
for target in tgts:
    metnet.target = target/100
    e_chems.append(metnet.FVA[metnet.chemical])
    dp = CB_P(network=metnet,k=K,log=False)
    ap = solver.CB_P(network=metnet,K=K)
    chems_function.append(dp.Vs[metnet.chemical])
    chems_solver.append(ap.Vs[metnet.chemical])

fig, ax = plt.subplots()
plt.plot(tgts,e_chems,linestyle='dotted',c='grey')
ax.set_title(f"{obj} Production K={K} - {name}")
plt.xlabel(r""" % Biomass min production $mmol/g(Dw)h$""")
plt.ylabel(r""" Chemical production $mmol/g(Dw)h$""")
plt.scatter(tgts,chems_function,marker='x',label="Function",s=65,c='black')
plt.scatter(tgts,chems_solver,marker='d',label="Algorithms",s=60,c='teal')
ax.legend()
plt.show()

