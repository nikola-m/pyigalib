# p14 - solve nonlinear BVP u_xx = -lambda*exp(u), u(-1)=u(1)=0 by iteration.
import numpy as np
from pyigalib import *
from nurbs_factory import *
from scipy.linalg import solve
from math import log10
import matplotlib.pyplot as plt

### DEFINE APPROXIMATING B-SPLINE AND COLLOCATION MATRICES ###
# B-spline order:
p=8
# Index of the highest control point:
n=100
# Generate (at the moment only uniform clamped) knot vector (it must be clamped):
uvec = make_knot_vector(p,n+1,"clamped")
#Scale knot vector to match the domain size:
a=0.; b=1 # domain size in 1d
uvec = np.array(uvec)
#print uvec
# scale to [-1,1]
#uvec = 2*uvec/float(uvec[p+n+1])-1. 
uvec = uvec/float(uvec[p+n+1])
# Generate collocation points in 1D - we choose them to be Greville apscissae:
xcp = greville(uvec,p,n+1)
#print xcp
# B-spline basis function and derivatives upt to 2nd order:
ders = basis_funs_and_derivatives_cpt(p,n,uvec,xcp,2)
# Matrix N_{i,j}==N_{j,p}(xcp(i)), (i=0,n+1,j=0,n+1) of basis functions, each row for each collocation point, each column for each basis fun.
N = ders[:,0,:]
# Matrix D2_{i,j}==d^2 N_{j,p}(xcp(i)) / dx^2, (i=0,n+1,j=0,n+1) of 2nd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D2 = ders[:,2,:]

### BOUNDARY CONDITIONS ###
# Apply boundary conditions (Homogenous Dirichlet):
D2 = D2[1:n,1:n]
N = N[1:n,1:n]

u=np.zeros(n-1)
err=np.zeros(n-1)
change = 0.1 # arbitrary small number
it = 0

while True:
  u=np.dot(N,u)
  unew=solve(D2,np.exp(u))
  err = np.abs(np.dot(N,unew)-u)
  change = err.max(); print change
  u = unew
  it += 1
  if(change < 1e-15):
    break

Px = xcp
Py = np.zeros(n+1)
# Find y-coordinates of control points from PDE:
Py[1:n] = u
# Create list of tuples (Px,Py) for B-spline plotting via nurbs_factory code
P = [(Px[i], Py[i]) for i in range(n+1)]

### POST PROCESSING ###

# Denser sample for plotting:
xx=np.linspace(a,b,100)
# Basis function for this vector of knot values
N = basis_funs_cpt(p,n,uvec,xx)
# Interpolate grid data
uu = np.dot(N,Py)    

# Analytical solution
cc=1.3360557
usol = -np.log(2)+2*np.log(cc*1./(np.cos(cc*(xx-0.5)/2.)))
maxerr = np.max(abs(uu-usol))

# Create a matplotlib figure
fig = plt.figure()
plt.title('p = %d, n = %d , maxerr = %e, no.steps = %d' % (p,n,maxerr,it))
fig.set_figwidth(22)
ax  = fig.add_subplot(111)
# Draw the curve points
ax.scatter( xx,uu, marker="o", c=xx, cmap="jet", alpha=0.5, label = "IgA collocation" )
# Draw the control cage.
ax.plot( xx, usol, alpha=0.7, label = "Analytical solution")
plt.legend()
plt.show()

