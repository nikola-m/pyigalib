import numpy as np

def C_factory(P, V=None, n=2):
    """ Returns a b-spline curve C(t) configured with P, V and n.
    
    Parameters
    ==========
    - P (list of D-tuples of reals) : List of de Boor points of dimension D.
    - n (int) : degree of the curve
    - V (list of reals) : list of knots in increasing order (by definition).
    
    Returns
    =======
    A D-dimensionnal B-Spline curve.
    """
    
    # TODO: check that p_len is ok with the degree and > 0
    m = len(P)    # the number of points in P    
    D = len(P[0]) # the dimension of a point (2D, 3D)
    
    # TODO: check the validity of the input knot vector.
    # TODO: create an initial Vector Point.
    
    #############################################################################
    # The following line will be detailed later.                                #
    # We create the highest degree basis spline function, aka. our entry point. #
    # Using the recursive formulation of b-splines, this b_n will call          #
    # lower degree basis_functions. b_n is a function.                          #
    #############################################################################
    b_n = basis_factory(n)
    
    def S(t, d):
        """ The b-spline funtion, as defined in eq. 3. """
        out = 0.
        for i in range(m): #: Iterate over 0-indexed point indices
            out += P[i][d]*b_n(t, i, V)
        return out
    
    def C(t):
        """ The b-spline curve, as defined in eq. 4. """
        out = [0.]*D           #: For each t we return a list of D coordinates
        for d in range(D):     #: Iterate over 0-indexed dimension indices
            out[d] = S(t,d)
        return out
    
    ####################################################################
    # "Enrich" the function with information about its "configuration" #
    ####################################################################    
    C.V = V                   #: The knot vector used by the function
    C.spline = S              #: The spline function.
    C.basis = b_n             #: The highest degree basis function. Useful to do some plotting.
    C.min = V[0]              #: The domain of definition of the function, lower bound for t
    C.max = V[-1]             #: The domain of definition of the function, upper bound for t
    C.endpoint = C.max!=V[-1] #: Is the upper bound included in the domain.
    return C

def basis_factory(degree):
    """ Returns a basis_function for the given degree """
    if degree == 0:
        
        def basis_function(t, i, knots):
            """The basis function for degree = 0 as per eq. 7"""
            t_this = knots[i]
            t_next = knots[i+1]
            out = 1. if (t>=t_this and t< t_next) else 0.         
            return out
    else:
        
        def basis_function(t, i, knots):
            """The basis function for degree > 0 as per eq. 8"""
            out = 0.
            t_this = knots[i]
            t_next = knots[i+1]
            t_precog  = knots[i+degree]
            t_horizon = knots[i+degree+1]            

            top = (t-t_this)
            bottom = (t_precog-t_this)
     
            if bottom != 0:
                out  = top/bottom * basis_factory(degree-1)(t, i, knots)
                
            top = (t_horizon-t)
            bottom = (t_horizon-t_next)
            if bottom != 0:
                out += top/bottom * basis_factory(degree-1)(t, i+1, knots)
         
            return out
        
    ####################################################################
    # "Enrich" the function with information about its "configuration" #
    ####################################################################         
    basis_function.lower = None if degree==0 else basis_factory(degree-1)
    basis_function.degree = degree
    return basis_function


def make_knot_vector(n, m, style="clamped"):
    """
    Create knot vectors for the requested vector type.
    
    Parameters
    ==========
    - n (int) : degree of the bspline curve that will use this knot vector
    - m (int) : number of vertices in the control polygone
    - style (str) : type of knot vector to output
    
    Returns
    =======
    - A knot vector (tuple)
    """
    if style != "clamped":
        raise NotImplementedError
        
    total_knots = m+n+2                           # length of the knot vector, this is exactly eq.6
    outer_knots = n+1                             # number of outer knots at each of the vector.
    inner_knots = total_knots - 2*(outer_knots)   # number of inner knots
    # Now we translate eq. 5:
    knots  = [0]*(outer_knots)
    knots += [i for i in range(1, inner_knots)]
    knots += [inner_knots]*(outer_knots)
    
    return tuple(knots) # We convert to a tuple. Tuples are hashable, required later for memoization

def greville(V,n,m):
    """
    Create Greville collocation points from the knot vector V.
    - n (int) : degree of the bspline curve
    - m (int) : number of vertices in control polygon  
    """
    out = np.zeros(m)
    for i in range(m):
        out[i] = np.sum(V[i+1:i+n+1])/float(n)
    return out

