import numpy as np
from scipy.optimize import minimize

# Initialize all variables here
R = 2; K = 2; 
N = [506, 664]; 
q_star = [1000, 2000];
rho = [ [0.035, 0.045], [0.020, 0.060] ]

i = 0; 
ci = [0.0] * K # target vector c
# objective function: cost minimization
def min_cost(ci):
  E = 0.0
  for k in range(K):
      E = E + ci[k]
  return E

# set boundary conditions for c
bnds = ((0.0, 2.0/rho[i][0]),)
for k in range(1, K):
  bnds = bnds + ((0.0, 2.0/rho[i][k]),)

# equality constraint
def broker_q_cons(ci):
  q = 0.0
  for x in range(K):
      for y in range(x, K):
          r = 1.0 if (x == y) else 0.5
          q = q + N[i] * r \
                       * np.sqrt((1-np.exp(-1.0*ci[x]*rho[i][x])) \
                                *(1-np.exp(-1.0*ci[y]*rho[i][y])))
  return q_star[i] - q
    
con1 = {'type': 'eq', 'fun': broker_q_cons}
cons = [con1]

solution = minimize(min_cost, ci, method='SLSQP', \
                    bounds=bnds, constraints=cons)

print (solution)