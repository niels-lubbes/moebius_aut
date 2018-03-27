'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''

from moebius_aut.class_ma_tools import MATools
from moebius_aut.class_ma_ring import ring
from moebius_aut.class_ma_ring import MARing

from moebius_aut.sage_interface import sage_matrix
from moebius_aut.sage_interface import sage_vector
from moebius_aut.sage_interface import sage_invariant_theory

class DSegre( object ):
    '''
    This class represents the Veronese-Segre embedding of P^1xP^1 
    We refer to such a surface as the "double Segre surface" and 
    it lives in projective 8-space P^8.
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
                
        If we compose this parametrization with for example 
        the projection 
            (1:x1:x2:x3:x4:x5:x6:x7:x8)
             |->
            (1:x1:x2:x3:x4:x7:x8)
        then we obtain a degree 6 Del Pezzo surface. 
        The ideal of this surface is obtained by setting 
        the input parameter "exc_idx_lst" to "[5,6]". 

        Parameters
        ----------
        exc_idx_lst: list<int> 
            A list of integers in [0,8].
        
        varname: ch     
            A character in [ 'x', 'y', 'z' ].
        
        Returns
        -------
        list
            A (sub-)list of generators for the ideal of the 
            double Segre surface. The ideal lives in a subring             
                QQ[x0,...,x8] (or QQ[y0,...,y8] or QQ[z0,...,z8])
            of the ring represented by "MARing.R". 
            The variable names are determined by "varname". 
            For each index i in "exc_idx_lst" the generators that contain 
            xi (or yi or zi) are omitted.                                                     
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
        The double Segre surface comes together with an antiholomorphic 
        involution. This method allows us to change coordinates so that 
        this antiholomorphic involution becomes complex conjugation.
                        
        A double Segre surface in projective 8-space P^8 is represented 
        by a toric parametrization whose completion to P^1xP^1 is provided 
        by "get_pmz_lst": 
            
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
                  
        An antiholomorphic involution, that preserves the toric structure,
        acts on the above lattice as a unimodular involution:
                  
            identity_*:    ( a, b ) |--> ( a, b)
            leftright_*:   ( a, b ) |--> (-a, b)
            rotate_*:      ( a, b ) |--> (-a,-b)  
            diagonal_*:    ( a, b ) |--> ( b, a)
        
        These unimodular lattice involutions induce an involution on P^8.
        
        We compose the toric parametrization of the double Segre surface 
        with the one of the following maps in order to make the antiholomorphic
        involution that acts on P^8, equal to complex conjugation.  

            identity: (x0:...:x8) |--> (x0:...:x8) 

            leftright:  x3 |--> x3,
                        x0 |--> x0,  
                        x4 |--> x4,                     
                        x1 |--> x1 + I*x2, 
                        x2 |--> x1 - I*x2, 
                        x5 |--> x5 + I*x8, 
                        x8 |--> x5 - I*x8,                                
                        x7 |--> x7 + I*x6,
                        x6 |--> x7 - I*x6 
                                 
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

    
        Parameters
        ----------
        iqf_lst: list<MARing>     
            A list of elements in the subring NF[x0,...,x8] 
            of "MARing.R" where NF denotes the Gaussian 
            rationals QQ(I) with I^2=-1.
                               
        involution: str  
            Either one of the following strings:
            'identity', 'leftright', 'rotate', 'diagonal'.

        Returns
        -------
            We consider the input "iqf_lst" as a map. 
            We compose this map composed with the map corresponding to 
            <identity>, <leftright>, <rotate> or <diagonal>.
            We return a list of elements in NF[x0,...,x8] that
            represents the composition. 
                                
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
            x3:z3,
            x4:z4,
            x5:z5 + I * z8,
            x8:z5 - I * z8,
            x1:z1 + I * z2,
            x2:z1 - I * z2,
            x7:z7 + I * z6,
            x6:z7 - I * z6 }

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
        Returns
        -------
        list<MARing>
            Returns a list of polynomials of bidegree (2,2)
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
    def get_aut_P8( c_lst ):
        '''
        The double Segre surface S is isomorphic to P^1xP^1.
        The pair (A,B) of 2x2 matrices denotes an automorphism 
        of P^1xP^1. We compute the representation of this 
        automorphism in P^8 by using the parametrization as 
        provided by ".get_pmz_lst". Since we consider the 
        2x2 matrices up to multiplication by a constant, 
        it follows that the automorphism group is 6-dimensional.        
        Formally, this method computes Sym^2(A)@Sym^2(B) 
        where @ denotes the tensor product (otimes in tex).        
        
        Parameters
        ---------- 
        c_lst: list<MARing.FF>
            A list of length 8 with elements c0,...,c7 in "MARing.FF". 
            We assume that the pair of matrices 
                ( [ c0 c1 ]   [ c4 c5 ] ) = (A,B) 
                ( [ c2 c3 ] , [ c6 c7 ] )                                                                                                                                                         
            represent an automorphism of P^1xP^1.
                                                             
        Returns
        -------
        sage_matrix
            A 9x9 matrix defined over "MARing.FF", which represents a 
            (parametrized) automorphism of P^8 that preserves the 
            double Segre surface S.                                                 
        '''
        # obtain parametrization in order to compute Sym^2(?)@Sym^2(?)
        #
        pmz_lst = DSegre.get_pmz_lst()

        # compute automorphisms double Segre surface
        #
        c0, c1, c2, c3, c4, c5, c6, c7 = c_lst
        x0, x1, y0, y1 = ring( 'x0,x1,y0,y1' )
        s, t, u, w = ring( 's,t,u,w' )  # coordinates of P^1xP^1
        dct1 = {}
        dct1.update( {s:c0 * x0 + c1 * x1} )
        dct1.update( {t:c2 * x0 + c3 * x1} )
        dct1.update( {u:c4 * y0 + c5 * y1} )
        dct1.update( {w:c6 * y0 + c7 * y1} )
        dct2 = {x0:s, x1:t, y0:u, y1:w}
        spmz_lst = [ pmz.subs( dct1 ).subs( dct2 ) for pmz in pmz_lst]

        # compute matrix from reparametrization "spmz_lst"
        # this is a representation of element in Aut(P^1xP^1)
        #
        mat = []
        for spmz in spmz_lst:
            row = []
            for pmz in pmz_lst:
                row += [spmz.coefficient( pmz )]
            mat += [row]
        mat = sage_matrix( MARing.FF, mat )

        MATools.p( 'c_lst =', c_lst )
        MATools.p( 'mat =\n' + str( mat ) )

        return mat


    @staticmethod
    def get_qmat( exc_idx_lst = [] ):
        '''
        We obtain generators of the ideal 
        of the (projection of) the double Segre surface
        with the method "get_ideal_lst( exc_idx_lst )".
        If the ideal is of a projection of the double Segre
        surface, then the returned matrix with parameters
        q0,...,q19 is not of full rank.
        
        Parameters
        ----------
        exc_idx_lst : list<int>
            A list of integers in [0,8].      
              
        Returns
        -------
        sage_matrix<MARing.R>
            A symmetric 9x9 matrix with entries
            in the ring QQ[q0,...,q19] which is a subring 
            of "MARing.R". It represents the Gramm matrix 
            of a quadratic form in the ideal of the              
            double Segre surface or a projection of 
            the double Segre surface with ideal defined
            by "get_ideal_lst( exc_idx_lst )".                        
        '''
        x = MARing.x()
        q = MARing.q()

        g_lst = DSegre.get_ideal_lst( exc_idx_lst )
        qpol = 0
        for i in range( len( g_lst ) ):
            qpol += q[i] * g_lst[i]

        qmat = sage_invariant_theory.quadratic_form( qpol, x ).as_QuadraticForm().matrix()
        qmat = sage_matrix( MARing.R, qmat )

        return qmat


    @staticmethod
    def get_invariant_q_lst( c_lst, exc_idx_lst = [] ):
        '''                        
        Parameters
        ----------
        c_lst : list<MARing.FF>
            A list of length 8 with elements 
                c0,...,c7 in QQ(k), 
            where QQ(k) is a subfield of "MARing.FF".
            If we substitute k:=0 in the entries of 
            "c_lst" then we should obtain the list:
                [1,0,0,1,1,0,0,1].                                                                                  
            A c_lst represents a pair of two matrices:                                
                ( [ c0 c1 ]   [ c4 c5 ] ) 
                ( [ c2 c3 ] , [ c6 c7 ] )                                   
            with the property that 
                c0*c3-c1*c2!=0 and c4*c7-c5*c6!=0. 
            If the two matrices are not normalized
            to have determinant 1, then the method should be 
            taken with care (it should be checked that the
            tangent vectors at the identity generate the 
            correct Lie algebra).                                                    
                            
        exc_idx_lst : list<int> 
            A list of integers in [0,8].                              
                                     
        Returns
        -------
        list<MARing.R>
            Let H be the representation of the pair of matrices  
               
                ( [ c0 c1 ]   [ c4 c5 ] ) 
                ( [ c2 c3 ] , [ c6 c7 ] )
               
            into P^8 (see also ".get_aut_P8()"). We assume here 
            that H is an element in Aut(P^1xP^1) and normalized 
            so that each 2x2 matrix has determinant 1.
                
            Thus H corresponds to a 1-parameter subgroup of Aut(P^8), 
            such that each automorphism preserves the double Segre 
            surface S in projective 8-space P^8.
               
            This method returns a list of generators of an ideal J 
            in the subring QQ[q0,...,q19] of "MARing.R". 
               
            Each point p in the zeroset V(J), when substituted in 
            the matrix 
                ".get_qmat(exc_idx_lst)",
            defines a quadratic form in the ideal 
                ".get_ideal_lst(exc_idx_lst)"
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
        A = DSegre.get_qmat( exc_idx_lst )
        Z = D.T * A + A * D
        iq_lst = [iq for iq in Z.list() if iq != 0 ]

        return iq_lst


    @staticmethod
    def get_invariant_qf( c_lst_lst, exc_idx_lst = [] ):
        '''
        Parameters
        ----------        
        c_lst_lst : list<list<MARing.FF>>  
            A list of "c_lst"-lists.
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
                c0*c3-c1*c2!=0 and c4*c7-c5*c6!=0.
            If the two matrices are not normalized
            to have determinant 1 then the method should be 
            taken with care (it should be checked that the
            tangent vectors at the identity generate the 
            correct Lie algebra).                                             
        
        exc_idx_lst : list<int>
            A list of integers in [0,8].                                     
        
        Returns
        -------
        A list of quadratic forms in the ideal of (a projection of) 
        the double Segre surface S:
            ".get_ideal_lst( exc_idx_lst )"
        such that the quadratic forms are invariant 
        under the automorphisms of S as defined by "c_lst_lst"
        and such that the quadratic forms generate the module of  
        all invariant quadratic forms. Note that Aut(S)=Aut(P^1xP^1).   
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
            iq_lst += DSegre.get_invariant_q_lst( c_lst, exc_idx_lst )
        iq_lst = list( MARing.R.ideal( iq_lst ).groebner_basis() )

        # solve the ideal defined by "iq_lst"
        #
        sol_dct = MARing.solve( iq_lst, q )

        # substitute the solution in the quadratic form
        # associated to the symmetric matrix qmat.
        #
        qmat = DSegre.get_qmat( exc_idx_lst )
        qpol = list( sage_vector( x ).row() * qmat * sage_vector( x ).column() )[0][0]
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
    def get_gens_sl2():
        '''
        Returns
        -------
        list<MARing.FF>
            A list of lists L of length 4 with elements c0,...,c3 in QQ(k)
            where the function field QQ(k) is a subfield of "MARing.FF".
            If we substitute k:=0 in the entries of L then we obtain 
            the list [1,0,0,1]. A list L represent a 2x2 matrix
                [ c0 c1 ]
                [ c2 c3 ]
            with the property that c0*c3-c1*c2!=0.
              
            The antiholomorphic involution coming from the real
            structure induces---up to conjugacy---two possible 
            antiholomorphic involutions acting on the 2x2-matrices:
                R0: [c0,c1,c2,c3] |--> [c0,c1,c2,c3] followed by complex conjugation.
                R1: [c0,c1,c2,c3] |--> [c3,c2,c1,c0] followed by complex conjugation.              
              
            A list L represents a real 1-parameter subgroup of PSL(2) wrt R0 or R1. 
            We assume that we are only interested in the tangent vector of this 
            subgroup at the identity. This tangent vector is an element in the 
            Lie algebra sl(2). 
            
            For example for the rotations wrt. R0 we have that  
                [ cos(k) -sin(k) ]
                [ sin(k)  cos(k) ]
            but we represent this element by "r": 
                [ 1  -k ]
                [ k   1 ] 
            because their tangent vectors at the identity coincide. 
            
            The 1-parameter subgroup of translations with real structure R0 
            corresponds to "t":
                [ 1  k ]
                [ 0  1 ]
            The corresponding element in the Lie algebra is 
                [ 0  1 ]
                [ 0  0 ]
            and is conjugate to 
                [  I   I ]
                [ -I  -I ]
            where I is the imaginary unit. This latter element is real wrt R1.
            A 1-parameter subgroup which has---up to scalar multiplication---this 
            element as tangent vector is "T":
                [  k+1   k         ]
                [ -k    (k+1)^(-1) ]
                
            See the code for an overview of generators L.
            The tangent vectors of t, q and s generate the Lie algebra sl(2).  
                                        
        Examples
        --------
            t, q, s, r, e, T = DSegre.get_gens_sl2()
            c_lst_lst = [s+e, e+s]
            iq_lst = DSegre.get_invariant_qf( c_lst_lst, [] )
            iq_lst = DSegre.change_basis( iq_lst, "rotate" )
            iq_lst = MARing.replace_conj_pairs( iq_lst )
            sig_lst = MARing.get_rand_sigs( iq_lst, 10 )            
        '''
        k = ring( 'k' )

        a = k + 1
        b = 1 / ( k + 1 )

        t = [1, k, 0, 1]  # translation for R0
        q = [1, 0, k, 1]
        s = [a, 0, 0, b]  # scalings for R0/rotations for R1
        r = [1, -k, k, 1]  # rotations for R0/scalings for R1
        e = [1, 0, 0, 1]  # identity

        T = [a, k, -k, b]  # translations for R1

        return t, q, s, r, e, T


    @staticmethod
    def get_c_lst_lst_lst():
        '''
        Returns
        -------
            A c_lst is a list of length 8 with elements c0,...,c7 in QQ(k), 
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
            Both matrices are the identity automorphism when k=0. 
            
            For example SO(2)xSO(2) in Aut(P^1xP^1) has two generators:
            
              c_lst_0 = [1, -k, k, 1, 1, 0, 0, 1]
              c_lst_1 = [1, 0, 0, 1, 1, -k, k, 1]

            We now define c_lst_lst = [c_lst_0, c_lst_1].
            
            The tangent vector of these two 1-parameter subgroups at the 
            identity determines an element in the Lie algebra sl2+sl2 
            of Aut(P^1xP^1):
             
                    [0, -1, 1, 0, 1, 0, 0, 1]
                    [1, 0, 0, 1, 0, -1, 1, 0]
            
            We assume that the antiholomorphic involution of the 
            real structure acts as complex conjugation on the matrices.
                          
            The output is a list of "c_lst_lst" elements. 
            Each "c_lst_lst" represent Lie subalgebra's of sl2+sl2 up to conjugacy.
            Up to flipping the left- and right-factor, exactly one representative 
            for each conjugacy class is contained in the output list.                                                                         
        '''
        #
        # obtain 1-parameter subgroups whose tangent vectors at the
        # identity generates the Lie algebra sl(2)
        #
        t, q, s, r, e, T = DSegre.get_gens_sl2()

        # shorthand notation
        #
        t1 = t + e; t2 = e + t
        q1 = q + e; q2 = e + q
        s1 = s + e; s2 = e + s
        r1 = r + e; r2 = e + r

        # construct dictionary for string representation
        #
        DSegre.__str_dct[str( t1 )] = 't1'
        DSegre.__str_dct[str( q1 )] = 'q1'
        DSegre.__str_dct[str( s1 )] = 's1'
        DSegre.__str_dct[str( r1 )] = 'r1'
        DSegre.__str_dct[str( t2 )] = 't2'
        DSegre.__str_dct[str( q2 )] = 'q2'
        DSegre.__str_dct[str( s2 )] = 's2'
        DSegre.__str_dct[str( r2 )] = 'r2'

        DSegre.__str_dct[str( t + t )] = 't1+t2'
        DSegre.__str_dct[str( q + q )] = 'g1+q2'
        DSegre.__str_dct[str( s + s )] = 's1+s2'
        DSegre.__str_dct[str( r + r )] = 'r1+r2'
        DSegre.__str_dct[str( r + s )] = 'r1+s2'
        DSegre.__str_dct[str( s + r )] = 's1+r2'
        DSegre.__str_dct[str( r + t )] = 'r1+t2'
        DSegre.__str_dct[str( t + r )] = 't1+r2'
        DSegre.__str_dct[str( s + t )] = 's1+t2'
        DSegre.__str_dct[str( t + s )] = 't1+s2'

        # construct classification of real Lie subalgebras
        #
        c_lst_lst_lst = []

        # Second projection Lie algebra is trivial
        c_lst_lst_lst += [[ t1, q1, s1 ]]
        c_lst_lst_lst += [[ t1, s1 ]]
        c_lst_lst_lst += [[ t1 ]]
        c_lst_lst_lst += [[ s1 ]]
        c_lst_lst_lst += [[ r1 ]]

        # Lie algebra is a product
        c_lst_lst_lst += [[ t1, q1, s1, t2, q2, s2 ]]
        c_lst_lst_lst += [[ t1, q1, s1, t2, s2 ]]
        c_lst_lst_lst += [[ t1, q1, s1, t2 ]]
        c_lst_lst_lst += [[ t1, q1, s1, s2 ]]
        c_lst_lst_lst += [[ t1, q1, s1, r2 ]]
        c_lst_lst_lst += [[ t1, s1, t2, s2 ]]
        c_lst_lst_lst += [[ t1, s1, t2 ]]
        c_lst_lst_lst += [[ t1, s1, s2 ]]
        c_lst_lst_lst += [[ t1, s1, r2 ]]
        c_lst_lst_lst += [[ t1, t2 ]]
        c_lst_lst_lst += [[ t1, s2 ]]
        c_lst_lst_lst += [[ t1, r2 ]]
        c_lst_lst_lst += [[ s1, s2 ]]
        c_lst_lst_lst += [[ s1, r2 ]]
        c_lst_lst_lst += [[ r1, r2 ]]

        c_lst_lst_lst += [[ t + t, q + q, s + s ]]
        c_lst_lst_lst += [[ t + t, s + s ]]
        c_lst_lst_lst += [[ t + t ]]
        c_lst_lst_lst += [[ s + s ]]
        c_lst_lst_lst += [[ r + r ]]
        c_lst_lst_lst += [[ t + s ]]
        c_lst_lst_lst += [[ t + r ]]
        c_lst_lst_lst += [[ s + r ]]

        c_lst_lst_lst += [[ s + s, t1, t2 ]]
        c_lst_lst_lst += [[ s + t, t1 ]]
        c_lst_lst_lst += [[ s + s, t1 ]]
        c_lst_lst_lst += [[ s + r, t1 ]]

        return c_lst_lst_lst


    @staticmethod
    def to_str( c_lst_lst ):
        '''
        Parameters
        ---------- 
            c_lst_lst : list 
                See output "get_c_lst_lst_lst()".
        
        Returns
        -------
        str
            A string representation for "c_lst_lst".
        '''
        if DSegre.__str_dct == {}:
            DSegre.get_c_lst_lst_lst()

        s = '< '
        for c_lst in c_lst_lst:
            s += DSegre.__str_dct[ str( c_lst ) ] + ', '
        s = s[:-2]
        s += ' >'

        return s











