# p16 - Poisson eq. on [-1,1]x[-1,1] with u=0 on boundary
from pyigalib import *
from nurbs_factory import *
from math import log10
import numpy as np
from scipy.linalg import solve,inv
from scipy.interpolate import interp2d
from matplotlib import pyplot as plt

import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import cm

from time import time

start = time()

### DEFINE APPROXIMATING B-SPLINE AND COLLOCATION MATRICES ###
# B-spline order:
p=2
# Index of the highest control point:
n=7
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

# Tensor product grid
ycp=xcp  
x,y = np.meshgrid(xcp[1:n], ycp[1:n])
x=x.flatten(1)
y=y.flatten(1)

# B-spline basis function and derivatives upt to 2nd order:
ders = basis_funs_and_derivatives_cpt(p,n,uvec,xcp,2)
# Matrix N_{i,j}==N_{j,p}(xcp(i)), (i=0,n+1,j=0,n+1) of basis functions, each row for each collocation point, each column for each basis fun.
N = ders[:,0,:]
# Matrix D2_{i,j}==d^2 N_{j,p}(xcp(i)) / dx^2, (i=0,n+1,j=0,n+1) of 2nd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D2 = ders[:,2,:]

### BOUNDARY CONDITIONS ###
# Apply boundary conditions (Homogenous Dirichlet):
N = N[1:n,1:n]
D2 = D2[1:n,1:n]

# RHS vector
pi=np.pi
f=2.*pi**2*np.sin(pi*x)*np.sin(pi*y)

# PDE operator-tensor product Laplacian
L=-(np.kron(N,D2)+np.kron(D2,N))
# condition number
#print p,n+1,np.linalg.cond(L) 
#sparsity pattern of L
#plt.title('Sparsity pattern of discretization matrix for Laplace operator, p = %d, n = %d'  % (p,n))
#plt.spy(L, precision=1e-50, marker='s', markersize=2)
#plt.show()

# Solve system
u=solve(L,f) 
elapsed = (time() - start)

# Reshape long 1D results to 2D grid:
Py=np.zeros((n+1,n+1))
Py[1:n,1:n] = u.reshape(n-1,n-1)

# Interpolate to finer grid for plotting
nsamples = 50
xnd = np.linspace(a,b,nsamples)  # equidistant grid
ynd=xnd
xx,yy = np.meshgrid(xnd,ynd)
xx=xx.flatten(1)
yy=yy.flatten(1)

# Basis functions matrix for sampled points
N = basis_funs_cpt(p,n,uvec,xnd)

# Napravi nurbs povrs:
uu=np.zeros((nsamples,nsamples))
# Sustina B-spline interpolacije u jednoj komandi
uu = np.dot(np.kron(N,N),Py.flatten(1))
exact=np.sin(pi*xx)*np.sin(pi*yy)

# Exact solution and Error
maxerr=max(abs(uu-exact))
#print p,n,maxerr,elapsed

# Prepare for plotting
uu=np.reshape(uu,(nsamples,nsamples))
xx=np.reshape(xx,(nsamples,nsamples))
yy=np.reshape(yy,(nsamples,nsamples))
# For control polygon:
x,y = np.meshgrid(xcp, ycp)

# Plot results
fig=plt.figure()
ax = p3.Axes3D(fig)
ax.plot_wireframe(x,y,Py,color='black', alpha=0.4)
ax.plot_surface(xx, yy, uu, rstride=1, cstride=1, cmap=cm.cool, linewidth=0.3, antialiased=True, alpha=0.7)
ax.text(-2.3,1.58,0.5,'p = %d, n = %d, (max err) = %e, time = %f [s]' % (p,n+1,maxerr,elapsed))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u')
plt.show()

# If you have it installed:
#from mayavi import mlab
#mlab.surf(xnd,ynd,uu)
#mlab.show()
#mlab.savefig('poisson-mlab.png')
