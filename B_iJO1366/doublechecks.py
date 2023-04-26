import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','MY_MODULES')))

from Ob_Met_Net_solmethods import CB_sol_OP

from MN_ijo1366 import MN_ijo1366
import matplotlib.pyplot as plt

metnet = MN_ijo1366

k = 2

tgt = 90
nocuts = []
cuts = []
pct = []

while tgt/100 > 0:
    pct.append(f"{int(tgt)}%")
    metnet.target=tgt/100
    print(f"minprod: {metnet.minprod}")
    c = CB_sol_OP(network=metnet,k=k,log=False)
    c1 = CB_sol_OP(network=metnet,k=k,log=False,extra=True)
    nocuts.append(c.Time)
    cuts.append(c1.Time)
    tgt -= 10

# print(f"Biomass: {c.Vs[metnet.biomass]}")
# print(f"Chemical: {c.Vs[metnet.chemical]}")
# print(f"Strategy: {c.Strategy}")
# print(f"Time: {c.Time}")

# print(f"\n Extra cuts:")
# print(f"Biomass: {c1.Vs[metnet.biomass]}")
# print(f"Chemical: {c1.Vs[metnet.chemical]}")
# print(f"Strategy: {c1.Strategy}")
# print(f"Time: {c1.Time}")
plt.plot(pct,nocuts,c='blue')
plt.plot(pct,cuts,c='green')
plt.show()