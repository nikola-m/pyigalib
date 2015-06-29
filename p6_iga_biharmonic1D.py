# p38 - solve u_xxxx = exp(x), u(-1)=u(1)=u'(-1)=u'(1)=0
#       compare with p13
#
import numpy as np
from pyigalib import *
from nurbs_factory import *
from scipy.linalg import solve
from math import log10
import matplotlib.pyplot as plt

### DEFINE APPROXIMATING B-SPLINE AND COLLOCATION MATRICES ###
# B-spline order:
p=9
# Index of the highest control point:
n=64
# Generate (at the moment only uniform clamped) knot vector (it must be clamped):
uvec = make_knot_vector(p,n+1,"clamped")
#Scale knot vector to match the domain size:
a=-1; b=1 # domain size in 1d
uvec = np.array(uvec)
#print uvec
# scale to [-1,1]
uvec = 2*uvec/float(uvec[p+n+1])-1. 
# Generate collocation points in 1D - we choose them to be Greville apscissae:
xcp = greville(uvec,p,n+1)
#print xcp
# B-spline basis function and derivatives upt to 2nd order:
ders = basis_funs_and_derivatives_cpt(p,n,uvec,xcp,4)
# Matrix N_{i,j}==N_{j,p}(xcp(i)), (i=0,n+1,j=0,n+1) of basis functions, each row for each collocation point, each column for each basis fun.
N = ders[:,0,:]
# Matrix D2_{i,j}==d N_{j,p}(xcp(i)) / dx, (i=0,n+1,j=0,n+1) of first derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D = ders[:,1,:]
# Matrix D2_{i,j}==d^2 N_{j,p}(xcp(i)) / dx^2, (i=0,n+1,j=0,n+1) of 2nd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D2 = ders[:,2,:]
# Matrix D3_{i,j}==d^3 N_{j,p}(xcp(i)) / dx^3, (i=0,n+1,j=0,n+1) of 3rd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D3 = ders[:,3,:]
# Matrix D4_{i,j}==d^4 N_{j,p}(xcp(i)) / dx^4, (i=0,n+1,j=0,n+1) of 4th derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D4 = ders[:,4,:]

D3=np.dot(D2,D)
D4=np.dot(D2,D2)

#v = np.zeros(n+1) # Work array
#v[1:n] = 1./(1.-xcp[1:n]**2)
#S = np.diag(v)
#D4 = np.dot((np.dot(np.diag(1.-xcp**2),D4) - 8*np.dot(np.diag(xcp),D3) - 12*D2),S)

D4 = (np.dot(np.diag(1.-xcp**2),D4) - 8*np.dot(np.diag(xcp),D3) - 12*D2)

D4 = D4[1:n,1:n]

# RHS vector
f=np.exp(xcp[1:n])

# Solve boundary-value problem
Py=np.zeros(n+1)
Py[0]=0.0
Py[n]=0.0
Py[1:n]=solve(D4,f)

# Interpolant
uu = np.zeros(n+1)
#uu = (1.-xcp**2)*np.dot(N,np.dot(S,Py))
uu = (1.-xcp**2)*np.dot(N,Py)             

# Exact solution and max err:
A = np.array([ [1,-1,1,-1],
               [0,1,-2,3],
               [1,1,1,1],
               [0,1,2,3] ])
V = np.vander(xcp)
V = V[:,n-3:n+1]
V = (V.T[::-1]).T # transpose  - reverse up-down - transpose
v = np.array([-1,-1,1,1])
c = solve(A,np.exp(v))
exact=(np.exp(xcp)-np.dot(V,c))

# maxerr
maxerr=max(abs(uu-exact))

# Interpolate to finer grid for plotting
nsamples = 100
xnd = np.linspace(a,b,nsamples)
# Basis functions matrix for sampled points
N = basis_funs_cpt(p,n,uvec,xnd)

#uuu = (1.-xnd**2)*np.dot(N,np.dot(S,Py))
uuu = (1.-xnd**2)*np.dot(N,Py)

plt.title('B-Spline order p = %d, number of control points: %d , maxerr = %e' % (p,n+1,maxerr))
plt.plot(xnd,uuu,'b')
plt.show()

