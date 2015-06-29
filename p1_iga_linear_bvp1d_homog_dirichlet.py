#
# Solving two-point BVP: -u''+u'+u = f, f= (1+4\pi^2)*sin(2\pi x)-2*\pi*cos(2\pi x), x \in [0,1], 
# Exact solution: u  = sin(2\pi x)
# Boundary conditions are homogenous Dirichlet u[0]=0, u[1]=0.
#
import numpy as np
from pyigalib import *
from nurbs_factory import *
from scipy.linalg import solve
from math import log10
import matplotlib.pyplot as plt

### DEFINE APPROXIMATING B-SPLINE AND COLLOCATION MATRICES ###
# B-spline order:
p=2
# Index of the highest control point:
n=100
# Generate (at the moment only uniform clamped) knot vector (it must be clamped):
uvec = make_knot_vector(p,n+1,"clamped")
#print uvec
#Scale knot vectro to match the domain size:
uvec = np.array(uvec)
#print uvec
uvec = uvec/float(uvec[p+n+1])
# Generate collocation points in 1D - we choose them to be Greville apscissae:
xcp = greville(uvec,p,n+1)
#print xcp
# B-spline basis function and derivatives upt to 2nd order:
ders = basis_funs_and_derivatives_cpt(p,n,uvec,xcp,2)
# Matrix N_{i,j}==N_{j,p}(xcp(i)), (i=0,n+1,j=0,n+1) of basis functions, each row for each collocation point, each column for each basis fun.
N = ders[:,0,:]
# Matrix D2_{i,j}==d N_{j,p}(xcp(i)) / dx, (i=0,n+1,j=0,n+1) of first derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D = ders[:,1,:]
# Matrix D2_{i,j}==d^2 N_{j,p}(xcp(i)) / dx^2, (i=0,n+1,j=0,n+1) of 2nd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D2 = ders[:,2,:]

### BOUNDARY CONDITIONS ###
# Apply boundary conditions (Homogenous Dirichlet):
N = N[1:n,1:n]
D = D[1:n,1:n]
D2 = D2[1:n,1:n]

# Create RHS vector for the given problem:
f=(1+4*np.pi**2)*np.sin(2*np.pi*xcp[1:n])+2*np.pi*np.cos(2*np.pi*xcp[1:n])

### SOLUTION ###
#Initialize solution
Px = xcp
Py = np.zeros(n+1)
# Find y-coordinates of control points from PDE:
Py[1:n] = solve((-D2+D+N),f) 
# Create list of tuples (Px,Py) for B-spline plotting via nurbs_factory code
P = [(Px[i], Py[i]) for i in range(n+1)]


### POST PROCESSING ###
# Create the Curve function
C = C_factory(P, uvec, p)

# Regularly spaced samples
xx = np.linspace(C.min, C.max, 100, endpoint=C.endpoint)

# Analytical solution
exact = np.sin(2*np.pi*xx)

# Sample the curve and make comparation with analytical solution
sampling = [t for t in xx]
curvepts = [ C(s) for s in sampling ]

# Max error
#print zip(*curvepts)[0] # samo x coordinatre tacaka
maxerr = np.max(np.abs(exact-zip(*curvepts)[1])) # samo y coordinatre tacaka
print p, n, maxerr

# Create a matplotlib figure
fig = plt.figure()
plt.title('B-Spline order p = %d, number of control points: %d , maxerr = %e' % (p,n+1,maxerr))
#plt.label('IgA','Analytical','Control polygon')
fig.set_figwidth(16)
ax  = fig.add_subplot(111)
# Draw the curve points
ax.scatter( *zip(*curvepts), marker="o", c=sampling, cmap="jet", alpha=0.5, label = "IgA collocation" )
# Draw analytical solution
#xx=np.linspace(0.,1.,100)
ax.plot( xx,exact, color="blue", alpha=0.7, label = "Analytical sol." )
# Draw the control cage.
ax.plot(*zip(*P), alpha=0.3, label = "Control polygon")
plt.legend()
plt.show()
