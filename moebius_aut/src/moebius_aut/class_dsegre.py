'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''
from sage.all import *
from class_ma_ring import ring
from class_ma_ring import MARing


class DSegre( object ):
    '''
    This class represents the Segre embedding of P^2xP^2 
    restricted to CxC where C is the Veronese embedding 
    of P^1 into P^2. We refer to such a surface as the 
    "double Segre surface" and it lives in P^8.
    The double Segre surface is the anticanonical model of a 
    Del Pezzo surface of degree 8.    
    '''

    @staticmethod
    def get_ideal_lst( exc_idx_lst = [], varname = 'x' ):
        '''
        We consider a parametrization of double Segre surface
        as provided by "get_pmz_lst". 
        
        (1:s:s^{-1}:u:u^{-1}:s*u:s^{-1}*u^{-1}:s*u^{-1}:s^{-1}*u)
        =
        (1:x1:x2:x3:x4:x5:x6:x7:x8)
        
        We can put the exponents of the monomials in a lattice
        where x0 corresponds to coordinate (0,0), x6 to (-1,-1)
        x5 to (1,1) and x8 to (-1,1). 
                
              x8 x3 x5
              x2 x0 x1
              x6 x4 x7  
        
        If we compose with the projection 
            (1:x1:x2:x3:x4:x5:x6:x7:x8)
             |->
            (1:x1:x2:x3:x4:x7:x8)
        then we obtain a degree 6 Del Pezzo surface. 
        The corresponding ideal is obtained by setting 
        the input parameter "exc_idx_lst" to "[5,6]". 

        INPUT:
            - "exc_idx_lst" -- A list of integers in [0,8].
            - "varname"     -- A character in [ "x", "y", "z" ].
        OUTPUT:
            - Returns a (sub-)list of generators for the ideal 
              of the double Segre surface.
              The ideal lives in a subring 
                  QQ[x0,...,x8] 
              of the ring represented by "MARing".
              For each index i in "exc_idx_lst" the 
              generators that contain xi are omitted.  
        '''

        s_lst = []
        s_lst += ['x0^2-x1*x2']
        s_lst += ['x0^2-x3*x4']
        s_lst += ['x0^2-x5*x6']
        s_lst += ['x0^2-x7*x8']

        s_lst += ['x1^2-x5*x7']
        s_lst += ['x2^2-x6*x8']
        s_lst += ['x3^2-x5*x8']
        s_lst += ['x4^2-x6*x7']

        s_lst += ['x0*x1-x4*x5']
        s_lst += ['x0*x2-x3*x6']
        s_lst += ['x0*x3-x2*x5']
        s_lst += ['x0*x4-x1*x6']

        s_lst += ['x0*x1-x3*x7']
        s_lst += ['x0*x2-x4*x8']
        s_lst += ['x0*x3-x1*x8']
        s_lst += ['x0*x4-x2*x7']

        s_lst += ['x0*x5-x1*x3']
        s_lst += ['x0*x6-x2*x4']
        s_lst += ['x0*x7-x1*x4']
        s_lst += ['x0*x8-x2*x3']

        s_lst = [ s.replace( 'x', varname ) for s in s_lst ]

        # exclude generators depending on exc_lst
        exc_lst = [ varname + str( i ) for i in exc_idx_lst ]
        for exc in exc_lst:
            s_lst = [ s for s in s_lst if str( exc ) not in s ]

        e_lst = [ ring( s ) for s in s_lst ]

        return e_lst


    @staticmethod
    def get_pmz_lst():
        '''
        OUTPUT:
            - Returns a list of polynomials of bidegree (2,2)
              in the subring QQ[s,t;u,w] of the ring "MARing".
        '''

        s_lst = []
        s_lst += [ 's*t*u*w' ]
        s_lst += [ 's*s*u*w' ]
        s_lst += [ 't*t*u*w' ]
        s_lst += [ 's*t*u*u' ]
        s_lst += [ 's*t*w*w' ]
        s_lst += [ 's*s*u*u' ]
        s_lst += [ 't*t*w*w' ]
        s_lst += [ 's*s*w*w' ]
        s_lst += [ 't*t*u*u' ]

        p_lst = [ ring( s ) for s in s_lst ]

        return p_lst


    @staticmethod
    def get_aut_P8():
        '''
        The double Segre surface S is isomorphic to P^1xP^1.
        The following pair of 2x2 matrices denotes an automorphism 
        of P^1xP^1:
              ( [ a b ]   [ e f ] ) 
              ( [ c d ] , [ g h ] )
        We compute the representation of this automorphism in P^8
        by using the parametrization as provided by ".get_pmz_lst".
        Since we consider the 2x2 matrix up to multiplication
        by a constant, the automorphism group is 6-dimensional.
        
        OUTPUT:
            - Returns a matrix with entries in the fraction field
              QQ(a,b,c,d,e,f,g,h). This parametrized matrix represents
              automorphisms of P^8 that preserve the double Segre 
              surface.  
        '''

        x0, x1, y0, y1 = ring( 'x0,x1,y0,y1' )
        a, b, c, d, e, f, g, h = ring( 'a,b,c,d,e,f,g,h' )
        s, t, u, w = ring( 's,t,u,w' )  # coordinates of P^1xP^1

        pmz_lst = DSegre.get_pmz_lst()

        # compute automorphisms double Segre surface
        #
        dct1a = {}
        dct1a.update( {s:a * x0 + b * x1} )
        dct1a.update( {t:c * x0 + d * x1} )
        dct1a.update( {u:e * y0 + f * y1} )
        dct1a.update( {w:g * y0 + h * y1} )

        dct1b = {x0:s, x1:t, y0:u, y1:w}
        spmz_lst = [ pmz.subs( dct1a ) for pmz in pmz_lst]
        spmz_lst = [ spmz.subs( dct1b ) for spmz in spmz_lst]

        mat = []
        for spmz in spmz_lst:
            row = []
            for pmz in pmz_lst:
                row += [spmz.coefficient( pmz )]
            mat += [row]
        mat = matrix( MARing.FF, mat )

        return mat


    @staticmethod
    def get_qmat_qpol():
        '''
        OUTPUT:
            - Returns a matrix and a quadric polynomial.
              
              * The matrix consists 
              
              
        '''
        x = MARing.x()
        q = MARing.q()

        g_lst = DASegre.get_ideal_lst()
        qpol = 0
        for i in range( len( g_lst ) ):
            qpol += q[i] * g_lst[i]

        qmat = invariant_theory.quadratic_form( qpol, x ).as_QuadraticForm().matrix()
        qmat = Matrix( MARing.R, qmat )

        return qmat, qpol


    def get_sig( pol ):

        x = AutRing.x()

        M = invariant_theory.quadratic_form( pol, x ).as_QuadraticForm().matrix()
        M = matrix( QQ, M )
        D, V = M.eigenmatrix_right()  # D has first all negative values on diagonal

        # determine signature of quadric
        #
        num_neg = len( [ d for d in D.diagonal() if d < 0 ] )
        num_pos = len( [ d for d in D.diagonal() if d > 0 ] )

        return ( num_neg, num_pos )



