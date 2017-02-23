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

    # private static dictionary for translating
    # c_lst_lst to a string representation (see ".to_str"())
    #
    __str_dct = {}


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
    def change_basis( iqf_lst, involution = 'identity' ):
        '''
        INPUT:
            - "iqf_lst"     -- A list of elements in the subring NF[x0,...,x8] 
                               of "MARing.R" where NF denotes the Guassian 
                               rationals QQ(I) with I^2=-1.
                               
            - "involution"  -- Either one of the following strings:
                               'identity', 'leftright', 'rotate', 'diagonal'.
                                
        OUTPUT:
            - We consider a toric parametrization of double Segre surface
              whose completion to P^1xP^1 is provided by "get_pmz_lst": 
            
                (s,u) |-->
                (1:s:s^{-1}:u:u^{-1}:s*u:s^{-1}*u^{-1}:s*u^{-1}:s^{-1}*u)
                =
                (x0:x1:x2:x3:x4:x5:x6:x7:x8)
            
              We can put the exponents of the monomials in a lattice
              where x0 corresponds to coordinate (0,0), x6 to (-1,-1)
              x5 to (1,1) and x8 to (-1,1): 
                    
                    x8 x3 x5
                    x2 x0 x1
                    x6 x4 x7
                  
              An antiholomorphic involution which preserves the toric structure
              acts on the above lattice as a unimodular involution:
                  
                    identity_*:    ( a, b ) |--> ( a, b)
                    leftright_*:   ( a, b ) |--> (-a, b)
                    rotate_*:      ( a, b ) |--> (-a,-b)  
                    diagonal_*:    ( a, b ) |--> ( b, a)
        
              These unimodular lattice involutions induce an involution on P^8: 

                    identity: (x0:...:x8) |--> (x0:...:x8) 

                    leftright:  x0 |--> x0,
                                x1 |--> x1 + I*x2, 
                                x2 |--> x1 - I*x2, 
                                x3 |--> x3,  
                                x4 |--> x4, 
                                x5 |--> x5 + I*x8, 
                                x6 |--> x7 - I*x6, 
                                x7 |--> x7 + I*x6, 
                                x8 |--> x5 - I*x8 

                    rotate: x0 |--> x0,
                            x1 |--> x1 + I*x2, x2 |--> x1 - I*x2, 
                            x3 |--> x3 + I*x4, x4 |--> x3 - I*x4, 
                            x5 |--> x5 + I*x6, x6 |--> x5 - I*x6, 
                            x7 |--> x7 + I*x8, x8 |--> x7 - I*x8      
                            
                            
                    diagonal:   x5 |--> x5,
                                x0 |--> x0, 
                                x6 |--> x6, 
                                x3 |--> x3 + I*x1,  
                                x1 |--> x3 - I*x1,
                                x8 |--> x8 + I*x7,  
                                x7 |--> x8 - I*x7,
                                x2 |--> x2 + I*x4,  
                                x4 |--> x2 - I*x4
                                
        '''

        I = ring( 'I' )
        x = x0, x1, x2, x3, x4, x5, x6, x7, x8 = MARing.x()
        z = z0, z1, z2, z3, z4, z5, z6, z7, z8 = MARing.z()

        dct = {}

        dct['identity'] = { x[i]:z[i] for i in range( 9 ) }

        dct['rotate'] = {
            x0:z0,
            x1:z1 + I * z2, x2:z1 - I * z2,
            x3:z3 + I * z4, x4:z3 - I * z4,
            x5:z5 + I * z6, x6:z5 - I * z6,
            x7:z7 + I * z8, x8:z7 - I * z8 }

        dct['leftright'] = {
            x0:z0,
            x1:z1 + I * z2,
            x2:z1 - I * z2,
            x3:z3, x4:z4,
            x5:z5 + I * z8,
            x6:z7 - I * z6,
            x7:z7 + I * z6,
            x8:z5 - I * z8  }

        dct['diagonal'] = {
                 x0:z0,
                 x6:z6,
                 x5:z5,
                 x3:z3 + I * z1,
                 x1:z3 - I * z1,
                 x8:z8 + I * z7,
                 x7:z8 - I * z7,
                 x2:z2 + I * z4,
                 x4:z2 - I * z4
                 }

        zx_dct = { z[i]:x[i] for i in range( 9 ) }

        new_lst = [ iqf.subs( dct[involution] ).subs( zx_dct ) for iqf in iqf_lst ]

        return new_lst


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
    def get_aut_P8( c_lst = None ):
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

        INPUT: 
            -- "c_lst" - Either None, or a list of length 8 with 
                         elements c0,...,c7 in "MARing.FF". We 
                         assume that the pair of matrices                                                                                 
                            ( [ c0 c1 ]   [ c4 c5 ] ) = (C,D) 
                            ( [ c2 c3 ] , [ c6 c7 ] )
                         represent an automorphism of P^1xP^1
                         such that both matrices are normalized 
                         to have determinant 1. 
                                                             
        OUTPUT:
            - This method returns a 9x9 matrix defined over "MARing.FF",
              which represents a (parametrized) automorphism of P^8
              that preserves the double Segre surface S. 
            
              Returns Sym^2(A)@Sym^2(B) if "c_lst==None".
              Returns Sym^2(C)@Sym^2(D) if "c_lst!=None".                                
        '''
        # obtain parametrization in order to compute Sym^2(?)@Sym^2(?)
        #
        pmz_lst = DSegre.get_pmz_lst()

        # compute automorphisms double Segre surface
        #
        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )
        x0, x1, y0, y1 = ring( 'x0,x1,y0,y1' )
        s, t, u, w = ring( 's,t,u,w' )  # coordinates of P^1xP^1
        dct1a = {}
        dct1a.update( {s:a * x0 + b * x1} )
        dct1a.update( {t:c * x0 + d * x1} )
        dct1a.update( {u:e * y0 + f * y1} )
        dct1a.update( {w:g * y0 + h * y1} )
        spmz_lst = [ pmz.subs( dct1a ) for pmz in pmz_lst]
        dct1b = {x0:s, x1:t, y0:u, y1:w}
        spmz_lst = [ spmz.subs( dct1b ) for spmz in spmz_lst]

        # compute matrix from reparametrization "spmz_lst"
        # this is a representation of element in Aut(P^1xP^1)
        #
        mat = []
        for spmz in spmz_lst:
            row = []
            for pmz in pmz_lst:
                row += [spmz.coefficient( pmz )]
            mat += [row]
        mat = matrix( MARing.FF, mat )

        # substitute
        #
        if c_lst != None:
            c0, c1, c2, c3, c4, c5, c6, c7 = c_lst
            subs_dct = {a:c0, b:c1, c:c2, d:c3, e:c4, f:c5, g:c6, h:c7}
            mat = mat.subs( subs_dct )

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
            -- "c_lst" - A list of length 8 with elements 
                         c0,...,c7 in QQ(k), 
                         where QQ(k) is a subfield of "MARing.FF".
                         If we substitute k:=0 in the entries of 
                         "c_lst" then we should obtain the list:
                                [1,0,0,1,1,0,0,1].                                                                                  
                         A c_lst represents a pair of two matrices:                                
                            ( [ c0 c1 ]   [ c4 c5 ] ) 
                            ( [ c2 c3 ] , [ c6 c7 ] )                                   
                         with the property that 
                             c0*c3-c1*c2=c4*c7-c5*c6=1.   
                                     
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
        # get representation of 1-parameter subgroup in Aut(P^8)
        #
        H = DSegre.get_aut_P8( c_lst )

        # consider the tangent vector of the curve H at the identity
        #
        k = ring( 'k' )
        D = MARing.diff_mat( H, k ).subs( {k:0} )

        # Note that if we differentiate the condition
        # A=H.T*A*H on both sides, evalute k=0, then
        # we obtain the condition D.T * A + A * D=0.
        # Here A denotes the matrix of a quadratic form
        # in the ideal of the double Segre surface S.
        #
        A = DSegre.get_qmat()
        Z = D.T * A + A * D
        iq_lst = [iq for iq in Z.list() if iq != 0 ]

        return iq_lst


    @staticmethod
    def get_invariant_ideal( c_lst_lst ):
        '''
        INPUT:
            -- "c_lst_lst"    - A list of "c_lst"-lists.
                                A c_lst is a list of length 8 with elements 
                                c0,...,c7 in QQ(k), 
                                where QQ(k) is a subfield of "MARing.FF".
                                If we substitute k:=0 in the entries of 
                                "c_lst" then we should obtain the list:
                                   [1,0,0,1,1,0,0,1].                                                                      
                                A c_lst represents a pair of two matrices:                                
                                   ( [ c0 c1 ]   [ c4 c5 ] ) 
                                   ( [ c2 c3 ] , [ c6 c7 ] )                                   
                                with the property that 
                                    c0*c3-c1*c2=c4*c7-c5*c6=1.            
                                    
        OUTPUT:
            -- A list of quadratic forms in the ideal of the double Segre 
               surface S, such that the quadratic forms are invariant 
               under the automorphisms of S as defined by "c_lst_lst"
               and such that the quadratic forms generate the ideal of  
               all invariant quadratic forms.   
                           
               For the ideal of all quadratic forms see ".get_ideal_lst()".              
        '''

        # for verbose output
        #
        mt = MATools()

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
        #
        sol_dct = MARing.solve( iq_lst, q )

        # substitute the solution in the quadratic form
        # associated to "get_q_mat()".
        #
        qmat = DSegre.get_qmat()
        qpol = list( vector( x ).row() * qmat * vector( x ).column() )[0][0]
        sqpol = qpol.subs( sol_dct )
        mt.p( 'sqpol   =', sqpol )
        mt.p( 'r       =', r )
        assert sqpol.subs( {ri:0 for ri in r} ) == 0
        iqf_lst = []  # iqf=invariant quadratic form
        for i in range( len( r ) ):
            coef = sqpol.coefficient( r[i] )
            if coef != 0:
                iqf_lst += [ coef ]
        mt.p( 'iqf_lst =', iqf_lst )

        return iqf_lst


    @staticmethod
    def get_c_lst_lst_dct():
        '''
        OUTPUT:
              
            - A c_lst is a list of length 8 with elements c0,...,c7 in QQ(k), 
              where QQ(k) is a subfield of "MARing.FF".
              If we substitute k:=0 in the entries of "c_lst" then we obtain 
              the list:
               
                  [1,0,0,1,1,0,0,1].                                          
              
              A c_lst represents a pair of two matrices:
              
                 ( [ c0 c1 ]   [ c4 c5 ] ) 
                 ( [ c2 c3 ] , [ c6 c7 ] )
                 
              with the property that c0*c3-c1*c2=c4*c7-c5*c6=1.
              
              This pair represent a 1-parameter subgroup G 
              (with parameter k) in Aut(P^1xP^1).
              Both matrices are normalized to have determinant 1 
              and is the identity automorphism at k=0. 
              
              For example SO(2)xSO(2) in Aut(P^1xP^1) has two generators:
              
                c_lst_0 = [k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1]
                c_lst_1 = [1, 0, 0, 1, k + 1, 0, 0, 1 / ( k + 1 )]

              We now define c_lst_lst = [c_lst_0, c_lst_1].
              
              The tangent vector of these two 1-parameter subgroups at the 
              identity determines an element in the Lie algebra sl2+sl2 
              of Aut(P^1xP^1):
               
                      [1, 0, 0, -1, 1, 0, 0, 1]
                      [1, 0, 0, 1, 1, 0, 0, -1]
              
              It is possible to classify all Lie subalgebra's
              of sl2+sl2 up to conjugacy.
                
              Each Lie subalgebra with n>0 generators determines a c_lst_lst of
              length n. Now "c_lst_lst_dct[n]" is a list of "c_lst_lst"-lists
              of lenght n.   
        
              Thus this methods returns a dictionary whose keys are positive 
              integers (1,2,3...) and each value is a list of c_lst_lst lists 
              of length equal to the corresponding key. Each c_lst_lst is a list 
              of "c_lst"-list and c_lst is as described above.                         
        '''
        # define c_list of all possible of 1-parameter
        # subgroups of Aut(P^1xP^1) up to conjugacy.
        #
        k = ring( 'k' )

        a = k + 1
        b = 1 / ( k + 1 )

        g1 = [1, k, 0, 1] + [1, 0, 0, 1]
        g2 = [1, 0, k, 1] + [1, 0, 0, 1]
        g3 = [a, 0, 0, b] + [1, 0, 0, 1]

        t1 = [1, 0, 0, 1] + [1, k, 0, 1]
        t2 = [1, 0, 0, 1] + [1, 0, k, 1]
        t3 = [1, 0, 0, 1] + [a, 0, 0, b]

        g1xt1 = [1, k, 0, 1] + [1, k, 0, 1]
        g2xt2 = [1, 0, k, 1] + [1, 0, k, 1]
        g3xt3 = [a, 0, 0, b] + [a, 0, 0, b]
        g1xt3 = [1, k, 0, 1] + [a, 0, 0, b]

        # construct dictionary for string representation
        #
        DSegre.__str_dct[str( g1 )] = 'g1'
        DSegre.__str_dct[str( g2 )] = 'g2'
        DSegre.__str_dct[str( g3 )] = 'g3'

        DSegre.__str_dct[str( t1 )] = 't1'
        DSegre.__str_dct[str( t2 )] = 't2'
        DSegre.__str_dct[str( t3 )] = 't3'

        DSegre.__str_dct[str( g1xt1 )] = 'g1xt1'
        DSegre.__str_dct[str( g2xt2 )] = 'g2xt2'
        DSegre.__str_dct[str( g3xt3 )] = 'g3xt3'
        DSegre.__str_dct[str( g1xt3 )] = 'g1xt3'


        # construct dictionary for output
        #
        c_lst_lst_dct = {}

        c_lst_lst_dct[1] = [ [g1], [g3], [g1xt1], [g3xt3], [g1xt3] ]

        c_lst_lst_dct[2] = [ [g1, g3], [g1, t1],
                             [g1, t3], [g3, t3],
                             [g1xt3, t1], [g3xt3, t1],
                             [g1xt1, g3xt3] ]

        c_lst_lst_dct[3] = [ [g1, g2, g3],
                             [g1, t1, t3],
                             [g3, t1, t3],
                             [g1xt1, g2xt2, g3xt3],
                             [g1, t1, g3xt3] ]

        c_lst_lst_dct[4] = [ [g1, g3, t1, t3],
                             [g1, g2, g3, t1],
                             [g1, g2, g3, t3] ]

        c_lst_lst_dct[5] = [ [g1, g2, g3, t1, t3] ]

        c_lst_lst_dct[6] = [ [g1, g2, g3, t1, t2, t3]]


        return c_lst_lst_dct


    @staticmethod
    def to_str( c_lst_lst ):
        '''
        INPUT: 
            - "c_lst_lst" -- See output "get_c_lst_lst_dct".
        OUTPUT:
            - A string representation for "c_lst_lst".
        '''
        if DSegre.__str_dct == {}:
            DSegre.get_c_lst_lst_dct()

        s = '< '
        for c_lst in c_lst_lst:
            s += DSegre.__str_dct[ str( c_lst ) ] + ', '
        s = s[:-2]
        s += ' >'

        return s











