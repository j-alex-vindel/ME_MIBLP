import sys
import gurobipy as gp
from gurobipy import GRB 

model = gp.read(sys.argv[1])
filename = sys.argv[1][-5:-3]
model.optimize()

if model.Status == GRB.INF_OR_UNBD:
    # Turn presolve off to determine whether model is infeasible
    # or unbounded
    model.setParam(GRB.Param.Presolve, 0)
    model.optimize()

if model.Status == GRB.OPTIMAL:
    print('Optimal objective: %g' % model.ObjVal)
    model.write(f"model_{filename}.sol")
    sys.exit(0)
elif model.Status != GRB.INFEASIBLE:
    print('Optimization was stopped with status %d' % model.Status)
    sys.exit(0)


# Model is infeasible - compute an Irreducible Inconsistent Subsystem (IIS)

print('')
print('Model is infeasible')
model.computeIIS()
model.write(f"pess_{filename}.ilp")
print("IIS written to file 'model.ilp'")