#
# Solving 1D diffusion equation IVP: u_t = u_xx, x \in [-5,5]
# Exact solution: u0  = 
# Boundary conditions are homogenous Dirichlet u[0]=0, u[1]=0.
#
import numpy as np
from pyigalib import *
from nurbs_factory import *
from scipy.linalg import solve, inv
from math import log10
import matplotlib.pyplot as plt


def rhs(A,u):
    """  Diffusion operator acting on points of control polygon 'py'
    """ 
    rhs = np.dot(A,u)

    return rhs

def ssprk54(t,dt,A,u):
    """SSPRK(5,4) The five-stage fourth-order SSPRK.
       SSPRK - Strong Stability Preserving Runge-Kutta
    """
    u1 = u + 0.391752226571890*dt*rhs(A,u)
    u2 = 0.444370493651235*u + 0.555629506348765*u1+0.368410593050371*dt*rhs(A,u1)
    u3 = 0.620101851488403*u + 0.379898148511597*u2+0.251891774271694*dt*rhs(A,u2)
    u4 = 0.178079954393132*u + 0.821920045606868*u3+0.544974750228521*dt*rhs(A,u3)
    u  = 0.517231671970585*u2 +0.096059710526147*u3 + 0.063692468666290*dt*rhs(A,u3) \
                              +0.386708617503269*u4 + 0.226007483236906*dt*rhs(A,u4)
    return t+dt,u

def rk2(t, dt, A, u):
    """ 2nd Order Runge-Kutta
        x - current time, like t
        h - like dt, a timestep
        u - initial value, here represented by beta coeff.
        f - derivative function u' = f(t, u)
    """
    up = rhs(A,u)

    u1 = u + 0.5*up         # odmah pripremam argument za sledeci korak, vektorka operacija

    up = rhs(A,u1)

    u = u + up*dt           # vektorska operacija 

    return t+dt, u

def rk4(t, dt, A, u):
    """ Forth order Runge-Kutta
        x - current time, like t
        h - like dt, a timestep
        u - initial value, here represented by beta coeff.
        f - derivative function u' = f(t, u)
    """

    k1 = dt * rhs(A,u)

    arg = u + 0.5*k1
    k2 = dt * rhs(A,arg)

    arg = u + 0.5*k2
    k3 = dt * rhs(A,arg)

    arg = u + 0.5*k3
    k4 = dt * rhs(A,arg)

    u = u + (k1 + 2*(k2 + k3) + k4)/6.0

    return t+dt, u


def main():

# Parameters for the problem:
    Dpar=1.; t0=0.1

# Solution domain:
    a = -5.
    b = 5.
# B-spline order:
    p=4
# Index of the highest control point:
    n=100
# Generate (at the moment only uniform clamped) knot vector (it must be clamped):
    uvec = make_knot_vector(p,n+1,"clamped")

# Scale knot vectro to match the domain size:
    uvec = np.array(uvec)
    uvec = (abs(a)+b)*uvec/float(uvec[p+n+1])-b
# Generate collocation points in 1D - we choose them to be Greville apscissae:
    x = greville(uvec,p,n+1)

# B-spline basis function and derivatives upt to 2nd order:
    ders = basis_funs_and_derivatives_cpt(p,n,uvec,x,2)
# Matrix N_{i,j}==N_{j,p}(xcp(i)), (i=0,n+1,j=0,n+1) of basis functions, each row for each collocation point, each column for each basis fun.
    N = ders[:,0,:]
# Matrix D2_{i,j}==d^2 N_{j,p}(xcp(i)) / dx^2, (i=0,n+1,j=0,n+1) of 2nd derivatives of basis functions, each row for each collocation point, each column for each basis fun.
    D2 = ders[:,2,:]

# Apply boundary conditions (Homogenous Dirichlet):
    D2 = D2[1:n,1:n]
   
# Solution parameters
    t = 0.1 # start time
    dt = 1e-3 # timestep-size
    tfinal = 1.0 #t + nsteps*dt # final time

# Set initial value of function
    u = np.exp(-x**2/(4*Dpar*t0))  
    py=solve(N,u);

    plt.plot(x, u,'mo');# plt.show()
    plt.hold(True)

# Matrix operators for 1-d diffusion IVP
    Ni=inv(N[1:n,1:n]); A=np.dot(Ni,D2)

    xx=np.linspace(a,b,200)
    uin = np.exp(-xx**2/(4*Dpar*t0)) 


# Do Runge-Kutta integration till final time
    while t < tfinal:
        #py[1:n] += dt*np.dot(A,py[1:n]); t += dt
#        print 'Calculating values for the new time step...'
#   Second order Runge-Kutta:
#        t, py[1:n] = rk2(t, dt, A, py[1:n])
#   Forth order Runge-Kutta:
        #t, py[1:n] = rk4(t, dt, A, py[1:n])
#   Five stage, forth order SSP Runge-Kutta
        t, py[1:n] = ssprk54(t, dt, A, py[1:n]); #py[0]=py[n] #<-periodic BC

    uex = np.zeros(n+1)
    usol = np.zeros(n+1)

    usol = np.dot(N,py)
    uex = np.sqrt(t0/t)*np.exp(-x**2/(4*Dpar*t))
    maxerr = abs(max(usol-uex))

# Basis functions matrix for sampled points
    N = basis_funs_cpt(p,n,uvec,xx)
    uu = np.dot(N,py)

### Plot solutions #####################
    plt.plot(xx, uin,'m', xx,uu,'c')
    #@plt.plot(x, usol,'cx-', x,uex,'b*-')
    plt.plot(x, usol,'c')
    plt.title('ssprk54, p=%d, n=%d, dt=%f, max err = %e' % (p,n,dt,maxerr))
    plt.legend(('initial','interpolant', 'final'))
    #plt.legend(('initial','final'))
    plt.grid()
    plt.axis('tight')
    plt.show()
 
if __name__ == "__main__":
    main()
