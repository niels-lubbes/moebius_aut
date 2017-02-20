'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''
from sage.all import *

from class_ma_tools import MATools
from class_ma_ring import ring
from class_ma_ring import MARing



class DSegre( object ):
    '''
    This class represents the Segre embedding of P^2xP^2 
    restricted to CxC where C is the Veronese embedding 
    of P^1 into P^2. We refer to such a surface as the 
    "double Segre surface" and it lives in projective 8-space P^8.
    The double Segre surface is the anticanonical model of a 
    Del Pezzo surface of degree 8.        
    '''

    @staticmethod
    def get_ideal_lst( exc_idx_lst = [], varname = 'x' ):
        '''
        We consider a toric parametrization of double Segre surface
        whose completion to P^1xP^1 is provided by "get_pmz_lst": 
        
        (s,u)
        |-->
        (1:s:s^{-1}:u:u^{-1}:s*u:s^{-1}*u^{-1}:s*u^{-1}:s^{-1}*u)
        =
        (1:x1:x2:x3:x4:x5:x6:x7:x8)
        
        We can put the exponents of the monomials in a lattice
        where x0 corresponds to coordinate (0,0), x6 to (-1,-1)
        x5 to (1,1) and x8 to (-1,1): 
                
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
              of the ring represented by "MARing.R".
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
              in the subring QQ[s,t;u,w] of the ring "MARing.R".
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
              ( [ a b ]   [ e f ] ) = (A,B)
              ( [ c d ] , [ g h ] )
        We compute the representation of this automorphism in P^8
        by using the parametrization as provided by ".get_pmz_lst".
        Since we consider the 2x2 matrix up to multiplication
        by a constant, the automorphism group is 6-dimensional.
        
        Formally, this method computes Sym^2(A)@Sym^2(B) 
        where @ denotes the tensor product (otimes in tex).
    
        
        OUTPUT:
            - Returns a matrix with entries in the fraction field
              QQ(a,b,c,d,e,f,g,h) (see MARing.FF). 
              This parametrized matrix represents automorphisms 
              of P^8 that preserve the double Segre surface S.  
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
    def get_qmat():
        '''
        OUTPUT:
            - Returns a symmetric 9x9 matrix with entries
              in the ring QQ[q0,...,q19] which is a subring 
              of "MARing.R". It represents the Gramm matrix 
              of a quadratic form in the ideal of the 
              double Segre surface.              
        '''
        x = MARing.x()
        q = MARing.q()

        g_lst = DSegre.get_ideal_lst()
        qpol = 0
        for i in range( len( g_lst ) ):
            qpol += q[i] * g_lst[i]

        qmat = invariant_theory.quadratic_form( qpol, x ).as_QuadraticForm().matrix()
        qmat = Matrix( MARing.R, qmat )

        return qmat


    @staticmethod
    def get_invariant_q_lst( c_lst ):
        '''                        
        INPUT: 
            -- "c_lst" - A list of length 8 with 
                         elements c0,...,c7 in QQ(k), 
                         where QQ(k) is a subfield of "MARing.FF".
                         If we substitute k:=0 in the entries of 
                         "c_lst" then we should obtain the list:
                             [1,0,0,1,1,0,0,1].
                                     
        OUTPUT:
            -- Let H be the representation of the pair
               of matrices  
               
                    ( [ c0 c1 ]   [ c4 c5 ] ) 
                    ( [ c2 c3 ] , [ c6 c7 ] )
               
               into P^8 (see also ".get_aut_P8()").
               We assume here that H is an element
               in Aut(P^1xP^1) and normalized so that
               each 2x2 matrix has determinant 1.
                
               Thus H corresponds to a 1-parameter subgroup
               of Aut(P^8), such that each automorphism preserves 
               the double Segre surface S in projective 8-space P^8.
               
               This method returns a list of generators 
               of an ideal J in the subring QQ[q0,...,q19] 
               of "MARing.R". 
               
               Each point p in the zeroset V(J),
               when substituted in the matrix ".get_qmat()",
               defines a quadratic form in the ideal ".get_ideal_lst()"
               that is preserved by the 1-parameter subgroup H.                                        
        '''
        # create a dictionary for substitution
        #
        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )
        k = ring( 'k' )
        c0, c1, c2, c3, c4, c5, c6, c7 = c_lst
        dct = {a:c0, b:c1, c:c2, d:c3, e:c4, f:c5, g:c6, h:c7}

        # get representation of 1-parameter subgroup in Aut(P^8)
        #
        H = DSegre.get_aut_P8().subs( dct )

        # consider the tangent vector of the curve H at the identity
        #
        D = MARing.diff_mat( H, k ).subs( {k:0} )

        # Note that if we differentiate the condition
        # A=H.T*A*H on both sides, evalute k=0, then
        # we obtain the condition D.T * A + A * D=0.
        # Here A denotes the matrix of a quadratic form
        # in the ideal of the double Segre surface S.
        #
        A = DSegre.get_qmat()
        Z = D.T * A + A * D

        return Z.list()


    @staticmethod
    def get_invariant_ideal( c_lst_lst ):
        '''
        INPUT:
            -- "c_lst_lst" - A list of lists c_lst, such that c_lst is 
                             a list of length 8 with 
                             elements c0,...,c7 in QQ(k), 
                             where QQ(k) is a subfield of "MARing.FF".
                             If we substitute k:=0 in the entries of 
                             "c_lst" then we should obtain the list:
                             [1,0,0,1,1,0,0,1].
        OUTPUT:
            -- A list of quadratic forms in the ideal of the double Segre 
               surface S, such that the quadratic forms are invariant 
               under the automorphisms of S as defined by "c_lst_lst"
               and such that the quadratic forms generate the ideal of  
               all invariant quadratic forms.   
               
               For the ideal of all quadratic forms see ".get_ideal_lst()" 
        '''
        # initialize vectors for indeterminates of "MARing.R"
        #
        x = MARing.x()
        q = MARing.q()
        r = MARing.r()

        # obtain algebraic conditions on q0,...,q19
        # so that the associated quadratic form is invariant
        # wrt. the automorphism defined by input "c_lst_lst"
        #
        iq_lst = []
        for c_lst in c_lst_lst:
            iq_lst += DSegre.get_invariant_q_lst( c_lst )
        iq_lst = list( MARing.R.ideal( iq_lst ).groebner_basis() )

        # solve the ideal defined by "iq_lst"
        # In order to use solve of Sage we cast to the symbolic ring "SR".
        #
        sol_dct = solve( [SR( str( iq ) ) for iq in iq_lst], [var( str( qi ) ) for qi in q], solution_dict = True )
        sol_dct = ring( sol_dct )[0]

        # substitute the solution in the quadratic form
        # associated to "get_q_mat()".
        #
        qmat = DSegre.get_qmat()
        qpol = list( vector( x ).row() * qmat * vector( x ).column() )[0][0]
        sqpol = qpol.subs( sol_dct )
        iqf_lst = []  # iqf=invariant quadratic form
        for i in range( len( r ) ):
            coef = sqpol.coefficient( r[i] )
            if coef != 0:
                iqf_lst += [ coef ]

        # verbose output
        #
        mt = MATools()
        mt.p( 'sqpol   =', sqpol )
        mt.p( 'iqf_lst =', iqf_lst )

        return iqf_lst







