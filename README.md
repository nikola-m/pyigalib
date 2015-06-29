# pyigalib
Isogeometric collocation in Python.
==================================================================

Basis of the library is Fortran code found in pyigalib.f95 file, which is turned into Python extension module using f2py:

``` 
f2py -c -m --fcompiler=gnu95 pyigalib pyigalib.f95
``` 

The result is pyigalib.so file. I'm using gfortran, therefore I use --fcompiler=gnu95 switch. You can no use pyigalib by importing it in Python:  

```python
from pyigalib imort *
from nurbs_factory import *
# B-spline order:
p=4
# Index of the highest control point:
n=20
# Generate (at the moment only uniform clamped) knot vector (it must be clamped):
uvec =  (p,n+1,"clamped")
#print uvec
#Scale knot vector to match the domain size:
uvec = np.array(uvec)
uvec = uvec/float(uvec[p+n+1])

# Generate collocation points in 1D - we choose them to be Greville apscissae:
xcp = greville(uvec,p,n+1)

# B-spline basis function and derivatives upt to 2nd order:
ders = basis_funs_and_derivatives_cpt(p,n,uvec,xcp,2)
# Matrix N_{i,j}==N_{j,p}(xcp(i)), (i=0,n+1,j=0,n+1) of basis functions, each row for each collocation point, each column for each basis fun.
N = ders[:,0,:]
# Matrix D_{i,j}==d N_{j,p}(xcp(i)) / dx, (i=0,n+1,j=0,n+1) of 1st derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D = ders[:,1,:]
# Matrix D2_{i,j}==d^2 N_{j,p}(xcp(i)) / dx^2, (i=0,n+1,j=0,n+1) of 2nd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
D2 = ders[:,2,:]
```

The 'nurbs_factory.py' file has some useful function such as 'make_knot_vector' or 'greville' and uses the code found [here](http://nbviewer.ipython.org/gist/dbarbeau/8b5ae150a65ce144a1bb).


A short instruction on how to use functions from pyigalib is found in pyigalib_doc file.

Nikola Mirkov
largeddysimulation@gmail.com
