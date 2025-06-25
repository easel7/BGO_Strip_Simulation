#!/usr/bin/python3

'''
Implementation of two and three exponent cases from Jean Jacquelin's
REGRESSION et EQUATION INTEGRALE, available at
https://www.scribd.com/doc/14674814/Regressions-et-equations-integrales
'''    

__author__    = "M Nelson"
__date__      = "January 19, 2020"
__version__   = "0.1"


import sys
import argparse

import numpy as np
from math import sqrt
from scipy.linalg import lstsq
from scipy.optimize import curve_fit

def biexpfit( xdata, ydata, debug=False ):

    x = np.array(xdata)
    y = np.array(ydata)
    if debug:
        print( 'x', x )
        print( 'y', y )

    S = np.empty_like(y)
    S[0] = 0
    S[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    if debug:
        print('S', S )

    SS = np.empty_like(y)
    SS[0] = 0
    SS[1:] = np.cumsum(0.5 * (S[1:] + S[:-1]) * np.diff(x))
    if debug:
        print('SS', SS )

    x2 = x * x
    x3 = x2 * x
    x4 = x2 * x2

    M = [ [sum(SS*SS),  sum(SS*S), sum( SS*x2 ), sum(SS*x), sum(SS)],
          [sum(SS*S),   sum(S*S),  sum(S*x2),    sum(S*x), sum(S) ],
          [sum(SS*x2),  sum(S*x2), sum(x4),      sum(x3),  sum(x2) ],
          [sum(SS*x),   sum(S*x),  sum(x3),      sum(x2),  sum(x) ],
          [sum(SS),     sum(S),    sum(x2),      sum(x),   len(xdata) ] ]

    if debug:
        print( 'M' )
        for n in range(5):
            print( M[n] )

    Ycol = np.array( [ sum(SS*y), sum(S*y), sum(x2*y), sum(x*y), sum(y) ] )

    (A,B,C,D,E),residues,rank,singulars = list( lstsq( M, Ycol ) )
    if debug:
        print( 'A-E', A, B, C, D, E )

    '''
    Minv = np.linalg.inv(M)    
    A,B,C,D,E = list( np.matmul(Minv,Ycol) )
    '''

    p = (1/2.)*(B + sqrt(B*B+4*A))
    q = (1/2.)*(B - sqrt(B*B+4*A))
    if debug:
        print( 'p,q', p, q )

    beta = np.exp(p*x)
    eta = np.exp(q*x)

    betaeta = beta * eta

    L = [ [ len(xdata), sum(beta), sum(eta) ],
          [ sum(beta),  sum(beta*beta), sum(betaeta) ],
          [ sum(eta),   sum(betaeta), sum(eta*eta)] ]

    Ycol = np.array( [ sum(y), sum(beta*y), sum(eta*y) ] )

    (a,b,c),residues,rank,singulars = list( lstsq( L, Ycol ) )    
    if debug:
        print( 'a,b,c', a,b,c )

    '''
    Linv = np.linalg.inv(L)
    a,b,c = list( np.matmul( Linv, Ycol ) )
    '''

    # sort in ascending order (fastest negative rate first)
    (b,p),(c,q) = sorted( [[b,p],[c,q]], key=lambda x: x[1])

    return a,b,c,p,q


def triexpfit( xdata, ydata, debug=False ):

    x = np.array(xdata)
    y = np.array(ydata)
    if debug:
        print( 'x', x )
        print( 'y', y )

    S = np.empty_like(y)
    S[0] = 0
    S[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    if debug:
        print('S', S )

    SS = np.empty_like(y)
    SS[0] = 0
    SS[1:] = np.cumsum(0.5 * (S[1:] + S[:-1]) * np.diff(x))
    if debug:
        print('SS', SS )

    SSS = np.empty_like(y)
    SSS[0] = 0
    SSS[1:] = np.cumsum(0.5 * (SS[1:] + SS[:-1]) * np.diff(x))
    if debug:
        print('SSS', SSS )

    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x
    x6 = x5 * x

    M = [ [ sum(SSS*SSS), sum(SSS*SS), sum(SSS*S), sum(SSS*x3), sum(SSS*x2), sum(SSS*x), sum(SSS) ],
          [ sum(SSS*SS), sum(SS*SS), sum(SS*S), sum(SS*x3), sum(SS*x2), sum(SS*x), sum(SS) ],
          [ sum(SSS*S), sum(SS*S), sum(S*S), sum(S*x3), sum(S*x2), sum(S*x), sum(S) ],
          [ sum(SSS*x3), sum(SS*x3), sum(S*x3), sum(x6), sum(x5), sum(x4), sum(x3) ],
          [ sum(SSS*x2), sum(SS*x2), sum(S*x2), sum(x5), sum(x4), sum(x3), sum(x2) ],
          [ sum(SSS*x), sum(SS*x), sum(S*x), sum(x4), sum(x3), sum(x2), sum(x) ],
          [ sum(SSS), sum(SS), sum(S), sum(x3), sum(x2), sum(x), len(xdata) ] ]

    if debug:
        print( 'M' )
        for n in range(7):
            print( M[n] )

    Ycol = [ sum(y*SSS), sum(y*SS), sum(y*S), sum(y*x3), sum(y*x2), sum(y*x), sum(y) ]

    if debug:
        print( 'Y', Ycol )

    '''
    Minv = np.linalg.inv(M)
    A,B,C,D,E,F,G = list( np.matmul( Minv, Ycol ) )
    print( 'A-G linalg', A, B, C, D, E, F, G )
    '''

    (A,B,C,D,E,F,G),residues,rank,singulars = list( lstsq( M, Ycol ) )
    if debug:
        print( 'A-G', A, B, C, D, E, F, G )

    p,q,r = np.roots( [1.,-C,-B,-A] )
    p,q,r = sorted( [p,q,r] )
    if debug:
        print( 'p,q,r', p, q, r )

    L = [ [ len(xdata), sum(np.exp(p*x)),sum(np.exp(q*x)),sum(np.exp(r*x)) ],
          [ sum(np.exp(p*x)), sum(np.exp(2*p*x)),sum(np.exp((p+q)*x)),sum(np.exp((p+r)*x)) ],
          [ sum(np.exp(q*x)), sum(np.exp((p+q)*x)),sum(np.exp(2*q*x)),sum(np.exp((q+r)*x)) ],
          [ sum(np.exp(r*x)), sum(np.exp((p+r)*x)),sum(np.exp((q+r)*x)),sum(np.exp(2*r*x)) ] ]

    Ycol = [ sum(y), sum(y*np.exp(p*x)),sum(y*np.exp(q*x)),sum(y*np.exp(r*x)) ]

    (a,b,c,d),residues,rank,singulars = list( lstsq( L, Ycol ) )
    if debug:
        print( 'a,b,c,d', a,b,c,d )


    # sort in ascending order (fastest negative rate first)
    (b,p),(c,q),(d,r) = sorted( [[b,p],[c,q],[d,r]], key=lambda x: x[1])

    return a,b,c,d, p,q,r


def __regression_tests__(debug=True):

    # x = np.linspace( 0, 2., 100 )

    # a, b = 1., 2.
    # p = -5.
    # print( 'inputs', a, b, p )

    # y = a + b * np.exp( p*x )

    # a, b, p = expfit( x, y )
    # print( 'fits', a, b, p )

    x = np.array([12.75, 38.25, 63.75, 89.25, 114.75, 140.25, 165.75, 191.25,
              216.75, 242.25, 267.75, 293.25, 318.75, 344.25])
    y = np.array([8721.5, 7603.5, 6627, 5819.5, 5068, 4432.5, 3876.5, 3402.5,
                2979.5, 2611, 2287, 2015, 1752, 1523])

    a, b, c = 1.0, 5000., 5000.
    p,q = -1/170., -1/200.
    print( 'inputs', a,b,c,p,q )

    y = a + b * np.exp( p*x ) + c * np.exp( q * x )

    a,b,c,p,q =  biexpfit( x, y, debug )
    print( 'fits Result', a,b,c, -1/p,-1/q )


    # a, b, c, d = -1., 5., 4., 2.
    # p, q, r = .5, -3., -2.
    # print( 'inputs', a,b,c,d, p,q,r )

    # x = np.linspace( .1, 1.5, 15 )
    # print( x )    

    # y = a + b * np.exp( p*x ) + c * np.exp( q* x ) + d * np.exp( r* x )
    # print( y )

    # a,b,c,d, p,q,r = triexpfit( x, y, debug )
    # print( 'fits', a,b,c,d, p,q,r )

if __name__ == "__main__":
    __regression_tests__(debug=True)