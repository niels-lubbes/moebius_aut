'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 17, 2017
@author: Niels Lubbes
'''
from sage.all import *

from class_ma_tools import MATools
from class_ma_ring import ring
from class_ma_ring import MARing
from class_dsegre import DSegre

mt = MATools()


def usecase__invariant_quadratic_forms( case ):
    '''
    INPUT:
        - "case" -- Three numerical characters, which characterize
                    the projection of a double Segre surface in S^n.
                    For example 
                        '078', 
                    means 
                        #families=0, #n=7, #degree=8.
                    If several examples are considered for the same
                    triple then we extend with a character in [a-z].
                    These are currently implemented cases:
                    
                    ['087','287','365','265s','265t','443','243ss','243st']
                    
    OUTPUT:
      - Let G be a subgroup of Aut(P^1xP^1)
        Compute the vectors space of real G-invariant quadratic forms in
        the ideal of a (projection of) the double Segre surface
        obtained by the method:
            "DSegre.get_ideal_lst( exc_idx_lst )"
                   
        The real structure is specified in terms of an involution. See 
            "DSegre.change_basis()" 
        for the specification of involutions.
        
        See in the code for a description of each of the cases of
        projections of double Segre surfaces that we consider.                
    '''

    #
    # declare generators of 1-parameter subgroups in k
    # such that k=0 is the identity.
    #
    k = ring( 'k' )
    I = ring( 'I' )
    U, V = k + 1, 1 / ( k + 1 )

    #
    # generators for real 1-parameter subgroups of Aut(P^1).
    #
    # For our algorithms only the tangent vectors of these
    # 1-parameter subgroups at the identity are relevant.
    #
    # For example if the involution is the identity then
    # the 1-parameter subgroup r has the same tangent vector
    # as the rotations [cos(k),-sin(k), sin(k), cos(k)].
    #
    # The antiholomorphic involution coming from the real
    # structure induces---up to conjugacy---two possible involutions
    # acting on the 2x2matrices:
    #     R1: complex conjugation
    #     R2: [a,b,c,d] --> [d,c,b,a] followed by complex conjugation.
    #
    # The 1-parameter subgroup i*T is conjugate to t and invariant
    # under involution R2.
    #
    # Note that t are scaling wrt. R1 but rotations wrt. R2.
    #
    t = [1, k, 0, 1]  # translation wrt. R1
    q = [1, 0, k, 1]  #
    s = [U, 0, 0, V]  # scalings (or rotations wrt. R2)
    r = [1, -k, k, 1]  # rotations (or scalings wrt. R2)
    e = [1, 0, 0, 1]  # identity
    T = [U, k, -k, V]  # translations wrt. R2

    #
    # "involution" is in { 'identity', 'leftright', 'rotate', 'diagonal' }
    #
    # "c_lst_lst" represents a subgroup G of Aut(P^1xP^1) as a list of
    # 1-parameter subgroups that generate G.
    #
    if case == '087':
        descr = '''
                This example shows that there exists a 
                quadric of signature (4,5) that is invariant
                under the full automorphism group of S.
                '''
        infoG = 'SL(2) x SL(2):  [t + e, q + e, s + e, e + t, e + q, e + s ]'
        exc_idx_lst = []
        c_lst_lst = [t + e, q + e, s + e, e + t, e + q, e + s ]
        involution = 'identity'


    elif case == '287':
        descr = '''
                This example shows quadrics of signatures (1,8), 
                (1,6) and (1,4) in the ideal of the double Segre
                surface, that are invariant under a 2-dimensional 
                group of toric Moebius automorphisms.         
                '''
        infoG = 'SO(2) x SO(2): [ s + e, e + s ]'
        exc_idx_lst = []
        c_lst_lst = [ s + e, e + s ]
        involution = 'rotate'


    elif case == '365':
        descr = '''
                This example shows that there exists a smooth sextic Del Pezzo 
                surface in S^5 that is invariant under a 2-dimensional group 
                of toric Moebius automorphisms.
                '''
        infoG = 'SO(2) x SO(2): [ s + e, e + s ]'
        exc_idx_lst = [5, 6]
        c_lst_lst = [ s + e, e + s ]
        involution = 'rotate'


    elif case == '265s':
        descr = ''' 
                This example shows that there exists a singular sextic 
                Del Pezzo surface in S^5 that is invariant under a 1-dimensional 
                group of toric Moebius automorphisms. Random signatures of type 
                [1,n+1] with n>=3 found so far are: [1,4], [1,6].

                If we extend the group with [e+s] then we obtain a quadric of 
                signature [1,4]; this is the horn cyclide case with 
                generators: x0^2 - x3*x4,  x0^2 - x1^2 - x2^2,  
                
                If we extend the group with [e+t] then we obtain a quadric of 
                signature [1,4]; this is the spindle cyclide case with 
                generators: x0^2 - x3*x4,  x4^2 - x6^2 - x7^2.                                  
                '''
        infoG = 'SO(2): [ s + e ]'
        exc_idx_lst = [5, 8]
        c_lst_lst = [ s + e ]
        involution = 'leftright'


    elif case == '265t':
        descr = ''' 
                This example indicates that there does not exists a singular 
                sextic Del Pezzo surface in S^5 that is invariant under a 
                1-dimensional group of Moebius translations. Random signatures 
                of type [1,n+1] with n>=3 found so far are: [1,4].                                
                '''
        infoG = 'SO(2): [ e + t ]'
        exc_idx_lst = [5, 8]
        c_lst_lst = [ e + t ]
        involution = 'leftright'


    elif case == '443':
        descr = '''
                This example shows that the Clifford torus in S^4 admits a 
                2-dimensional family of toric Moebius automorphisms. 
                '''
        infoG = 'SO(2) x SO(2): [ s + e, e + s ]'
        exc_idx_lst = [5, 6, 7, 8]
        c_lst_lst = [ s + e, e + s ]
        involution = 'rotate'

    elif case == '243ss':
        descr = ''' 
                This example shows that the horn cyclide in S^4 admits a 
                2-dimensional family of Moebius automorphisms.
                '''
        infoG = 'SO(2) x SO(2): [ s + e, e + s ]'
        exc_idx_lst = [5, 6, 7, 8]
        c_lst_lst = [ s + e, e + s ]
        involution = 'leftright'

    elif case == '243st':
        descr = ''' 
                This example shows that the spindle cyclide in S^4 admits a 
                2-dimensional family of Moebius automorphisms.
                '''
        infoG = 'SO(2) x SE(1): [ s + e, e + t ]'
        exc_idx_lst = [1, 2, 5, 8]
        c_lst_lst = [ s + e, e + t ]
        involution = 'leftright'

    else:
        raise ValueError( 'Unknown case: ', case )

    #
    # compute vector space of invariant quadratic forms
    #
    iq_lst = DSegre.get_invariant_qf( c_lst_lst, exc_idx_lst )
    iq_lst = DSegre.change_basis( iq_lst, involution )
    iq_lst = MARing.replace_conj_pairs( iq_lst )

    #
    # computes signatures of random quadrics in
    # the vector space of invariant quadratic forms.
    #
    sig_lst = MARing.get_rand_sigs( iq_lst, 10 )

    #
    # output results
    #
    while descr != descr.replace( '  ', ' ' ):
        descr = descr.replace( '  ', ' ' )
    mat_str = str( matrix( [[8, 3, 5], [2, 0, 1], [6, 4, 7]] ) )
    new_mat_str = mat_str
    for ei in exc_idx_lst:
        new_mat_str = new_mat_str.replace( str( ei ), ' ' )
    for ei in range( 0, 9 ):
        new_mat_str = new_mat_str.replace( str( ei ), '*' )

    mt.p( '\n' + 80 * '-' )
    mt.p( 'case        :', case )
    mt.p( 'description :\n', descr + '\n' + mat_str + '\n\n' + new_mat_str )
    mt.p( 'G           :', infoG )
    mt.p( 'exc_idx_lst :', exc_idx_lst )
    mt.p( 'involution  :', involution )
    mt.p( 'G-invariant quadratic forms:' )
    for iq in iq_lst:
        mt.p( '\t', iq )
    mt.p( 'random signatures:' )
    mt.p( '\t', sig_lst )
    for sig in sig_lst:
        if 1 in sig:
            mt.p( '\t', sig )
    mt.p( '\n' + 80 * '-' )


def usecase__toric_invariant_celestials():
    '''
    OUTPUT:
      - Suppose that G is the group of toric automorphisms
        of the double Serge surface X in P^8. The double 
        Segre surface admits a monomial parametrization. 
        The exponents corresponds to lattice points on 
        a square. We number these lattice points as follows:
        
            8  3  5
            2  0  1  
            6  4  7
        
        These exponents correspond to following parametrization:  
        
        (s,u) |-->
            (1:s:s^{-1}:u:u^{-1}:s*u:s^{-1}*u^{-1}:s*u^{-1}:s^{-1}*u)
            =
            (x0:x1:x2:x3:x4:x5:x6:x7:x8)
        
        The identity component of the toric automorphisms G are 
        represented by 9x9 diagonal matrices and are isomorphic to the 
        group SO(2)xSO(2). The vector space of G-invariant quadratic 
        forms in the ideal of the double Segre surface in P^8 is 
        generated by:
    
          V = < x0^2 - x7^2 - x8^2, 
                x0^2 - x5^2 - x6^2, 
                x0^2 - x3^2 - x4^2, 
                x0^2 - x1^2 - x2^2 >
        
        A quadric Q in V of signature (1,n+1) corresponds to 
        a celestial Y in the unit sphere S^n as follows: S^n
        is up to real projective equivalence defined by the
        projection of Q from its singular locus. The celestial
        Y is the image of the double Segre surface X under this 
        projection. Note that the singular locus of Q can intersect
        X such that Y is of lower degree. Another possibility is
        that the projection restricted to X is 2:1 such that the 
        degree of Y is 4. It turns out that projection of X with
        center singular locus of Q is toric in the sense we obtain
        ---up to projective equivalence---a parametrization by
        omitting monomials in the parametrization above. Using the 
        labels 0-8 of the lattice points in the square we consider
        the following projections, where P denotes the list of 
        monomials which are omitted:
        
        Q = (x0^2-x7^2-x8^2)+(x0^2-x5^2-x6^2)  
        P = [1,2,3,4]
        
        Q = (x0^2-x7^2-x8^2)+(x0^2-x3^2-x4^2)  
        P = [1,2,5,6]
        
        Q = (x0^2-x7^2-x8^2)+(x0^2-x1^2-x2^2)  
        P = [1,2,3,4]
        
        Q = (x0^2-x7^2-x8^2)+(x0^2-x3^2-x4^2)  
        P = [1,2,5,6]        
                                 
    '''

    # ideal of the double Segre surface in
    # coordinate ring R of P^8
    #
    IX = DSegre.get_ideal_lst()
    R = PolynomialRing( QQ, ['x' + str( i ) for i in range( 9 )] )
    IX = R.ideal( sage_eval( str( IX ), R.gens_dict() ) )

    # we consider all possible projections from the singular
    # locus of Q (see OUPUT documention).
    #
    prevP_lst = []
    for P1 in [[1, 2], [3, 4], [5, 6], [7, 8]]:
        for P2 in [[1, 2], [3, 4], [5, 6], [7, 8]]:

            # construct indices for projection map
            #
            P = sorted( list( set( P1 + P2 ) ) )
            if P not in prevP_lst:
                prevP_lst += [P]
            else:
                continue

            # project the double Segre surface wrt. P
            #
            J = IX.elimination_ideal( [sage_eval( 'x' + str( i ), R.gens_dict() ) for i in P] )

            # in order to get the right Hilbert polynomial we coerce into smaller ring
            #
            RP = PolynomialRing( QQ, ['x' + str( i ) for i in range( 9 ) if i not in P] )
            JP = RP.ideal( sage_eval( str( J.gens() ), RP.gens_dict() ) )
            hpol = JP.hilbert_polynomial()

            # output info
            #
            if not JP.is_zero():
                mt.p( 10 * '-' )
                mt.p( 'set of i such that <xi for i in P> is the ideal of the center of projection:' )
                mt.p( '\t', P )
                mt.p( 'projection of double Segre surface wrt. P:' )
                for gen in JP.gens():
                    mt.p( '\t', gen )
                    for involution in ['leftright', 'rotate']:
                        mt.p( '\t\t', involution, '\t', DSegre.change_basis( [ring( gen )], involution ) )
                mt.p( 'hilbert polynomial ((deg/dim!)*t^(dim)+...):' )
                mt.p( '\t', hpol )


def usecase__horn_and_spindle_cyclides():
    '''
    OUTPUT:
        We consider the following cyclides with monomial
        parametrization determined by the following lattice 
        polygons:  
        
            [  *  ]        [  *  ]
            [* * *]        [  *  ]
            [  *  ]        [* * *]
            
         horn cyclide    spindle cyclide   
         leftright       leftright  
         (2,4,3)         (2,4,3)
           
        and 'leftright' involution act as a modular involution
        with vertical symmetry axis. We use the following numbering 
        (see DSegre.get_ideal_lst() and DSegre.change_basis()):
                   
            [8 3 5]
            [2 0 1]
            [6 4 7]                
    '''
    if False:
        #
        # G-invariant quadratic forms for horn cyclide
        # (see usecase__invariant_quadratic_forms('243ss')):
        #
        # x0^2 - x3*x4
        # x0^2 - x1^2 - x2^2
        #
        # We send x3 |--> a*(x4-x3)
        #         x4 |--> a*(x4+x3)
        #         where a=sqrt(1/2)
        # followed by sending: x4 |--> x0, x0 |--> x4

        a = ring( 'a' )
        xv = x0, x1, x2, x3, x4 = ring( 'x0,x1,x2,x3,x4' )
        y0, y1, y2, y3 = ring( 'y0,y1,y2,y3' )
        m34 = matrix( MARing.R, [
            [1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, -a, a],
            [0, 0, 0, a, a]] )
        m04 = matrix( MARing.R, [
            [0, 0, 0, 0, 1],
            [0, 1, 0, 0, 0],
            [0, 0, 1, 0, 0],
            [0, 0, 0, 1, 0],
            [1, 0, 0, 0, 0]] )
        m = m34 * m04
        v = vector( xv )
        g1 = ring( 'x0^2 - x3*x4' )
        g2 = ring( 'x0^2 - x1^2 - x2^2' )
        mat1 = invariant_theory.quadratic_form( g1, xv ).as_QuadraticForm().matrix()
        mat2 = invariant_theory.quadratic_form( g2, xv ).as_QuadraticForm().matrix()
        ng1 = v.row() * m.T * mat1 * m * v.column()
        ng2 = v.row() * m.T * mat2 * m * v.column()
        G1 = SR( str( ng1[0] ) ).subs( {SR( 'a' ):1 / sqrt( 2 )} )
        G2 = SR( str( ng2[0] ) ).subs( {SR( 'a' ):1 / sqrt( 2 )} )

        S = PolynomialRing( QQ, var( 'x0,x1,x2,x3,x4,y0,y1,y2,y3' ) )
        x0, x1, x2, x3, x4, y0, y1, y2, y3 = S.gens()
        G1 = sage_eval( str( G1 ), S.gens_dict() )
        G2 = sage_eval( str( G2 ), S.gens_dict() )
        assert G1 in S
        assert G2 in S
        smap = [y0 - ( x0 - x3 ), y1 - x1, y2 - x2, y3 - x4]
        prj_lst = ideal( [G1, G2] + smap ).elimination_ideal( [x0, x1, x2, x3, x4] ).gens()
        eqn = str( prj_lst[0].subs( {y0:1} ) ).replace( 'y1', 'x' ).replace( 'y2', 'y' ).replace( 'y3', 'z' )

        mt.p( 80 * '-' )
        mt.p( 'HORN CYCLIDE' )
        mt.p( 80 * '-' )
        mt.p( 'We define a projective automorphism m:P^4--->P^4' )
        mt.p( '\t a   = 1/sqrt(2)' )
        mt.p( '\t m34 =', list( m34 ) )
        mt.p( '\t m04 =', list( m04 ) )
        mt.p( '\t m = m34 * m04 =', list( m ) )
        mt.p( '\t det(m) =', det( m ) )
        mt.p( '\t v |--> m*v =', v, '|-->', m * v )
        mt.p( 'Generators of ideal of cyclide in quadric of signature [1,4]:' )
        mt.p( '\t g1 =', g1 )
        mt.p( '\t g2 =', g2 )
        mt.p( 'Generators of ideal of cyclide in S^3 after applying m:' )
        mt.p( '\t G1 =', G1 )
        mt.p( '\t G2 =', G2 )
        mt.p( '\t G2-2*G1 =', G2 - 2 * G1 )
        mt.p( 'Stereographic projection to circular quadratic cone:' )
        mt.p( '\t smap    =', smap, '(stereographic projection map)' )
        mt.p( '\t prj_lst =', prj_lst )
        mt.p( '\t eqn     =', eqn )
        mt.p( 80 * '-' + 2 * '\n' )

    if True:
        #
        # G-invariant quadratic forms for spindle cyclides
        # (see usecase__invariant_quadratic_forms('243st')):
        #
        # x0^2 - x3*x4
        # x4^2 - x6^2 - x7^2
        #

        a = ring( 'a' )
        xv = x0, x3, x4, x6, x7 = ring( 'x0,x3,x4,x6,x7' )
        y0, y1, y2, y3 = ring( 'y0,y1,y2,y3' )
        m = matrix( MARing.R, [
            [a, 0, 0, 0, 0],
            [0, 1, 0, 0, 0],
            [0, -1, -1, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1]] )
        v = vector( xv )
        g1 = ring( 'x0^2 - x3*x4' )
        g2 = ring( 'x4^2 - x6^2 - x7^2' )
        mat1 = invariant_theory.quadratic_form( g1, xv ).as_QuadraticForm().matrix()
        mat2 = invariant_theory.quadratic_form( g2, xv ).as_QuadraticForm().matrix()
        ng1 = v.row() * m.T * mat1 * m * v.column()
        ng2 = v.row() * m.T * mat2 * m * v.column()
        G1 = SR( str( ng1[0] ) ).subs( {SR( 'a' ):1 / sqrt( 2 )} )
        G2 = SR( str( ng2[0] ) ).subs( {SR( 'a' ):1 / sqrt( 2 )} )

        S = PolynomialRing( QQ, var( 'x0,x3,x4,x6,x7,y0,y1,y2,y3' ) )
        x0, x3, x4, x6, x7, y0, y1, y2, y3 = S.gens()
        G1 = sage_eval( str( G1 ), S.gens_dict() )
        G2 = sage_eval( str( G2 ), S.gens_dict() )
        assert G1 in S
        assert G2 in S
        smap = [y0 - ( x4 + x3 ), y1 - x0, y2 - x7, y3 - x6]  # +x4^2-x0^2-x3^2-x6^2-x7^2
        prj_lst = S.ideal( [G1, G2] + smap ).elimination_ideal( [x0, x3, x4, x6, x7] ).gens()
        eqn = str( prj_lst[0].subs( {y0:1} ) ).replace( 'y1', 'x' ).replace( 'y2', 'y' ).replace( 'y3', 'z' )


        mt.p( 80 * '-' )
        mt.p( 'SPINDLE CYCLIDE' )
        mt.p( 80 * '-' )
        mt.p( 'We define a projective automorphism m:P^4--->P^4' )
        mt.p( '\t a   = 1/sqrt(2)' )
        mt.p( '\t m =', list( m ) )
        mt.p( '\t det(m) =', det( m ) )
        mt.p( '\t v |--> m*v =', v, '|-->', m * v )
        mt.p( 'Generators of ideal of cyclide in quadric of signature [1,4]:' )
        mt.p( '\t g1 =', g1 )
        mt.p( '\t g2 =', g2 )
        mt.p( 'Generators of ideal of cyclide in S^3 after applying m:' )
        mt.p( '\t G1 =', G1 )
        mt.p( '\t G2 =', G2 )
        mt.p( '\t G2-2*G1 =', G2 - 2 * G1 )
        mt.p( 'Stereographic projection to circular cylinder:' )
        mt.p( '\t smap    =', smap, '(stereographic projection map)' )
        mt.p( '\t prj_lst =', prj_lst )
        mt.p( '\t eqn     =', eqn )
        mt.p( 80 * '-' + 2 * '\n' )



def usecase__complex_classification():
    '''
    OUTPUT:
        - Prints a classification of quadratic forms
          that contain some fixed double Segre surface
          in projective 8-space, such that the quadratic
          form is invariant under subgroup of Aut(S).
          We consider subgroups equivalent, if they are
          complex conjugate in Aut(P^1xP^1).  
    '''


    c_lst_lst_dct = DSegre.get_c_lst_lst_dct()
    for key in c_lst_lst_dct:

        mt.p( 10 * '=' )
        mt.p( 'dimension =', key )
        mt.p( 10 * '=' )

        for c_lst_lst in c_lst_lst_dct[key]:

            mt.p( '\t', 10 * '-' )
            mt.p( '\t', 'group           =', DSegre.to_str( c_lst_lst ) )
            mt.p( '\t', 10 * '-' )

            for involution in ['identity', 'leftright', 'rotate']:

                J = DSegre.get_invariant_qf( c_lst_lst )
                J = DSegre.change_basis( J, involution )

                mt.p( '\t', 'involution      =', involution )
                mt.p( '\t', 'invariant ideal = <' )
                for gen in J:
                    sig = ''
                    if 'I' not in str( gen ):
                        sig = MARing.get_sig( gen )
                    mt.p( '\t\t', gen, '\t\t', sig )
                mt.p( '\t\t', '>' )
                mt.p( '\t', 5 * '.-' )

            mt.p( '\t', 10 * '-' )


if __name__ == '__main__':

    mt.start_timer()
    mt.filter( '__main__.py' )  # output only from this module
    set_verbose( -1 )  # surpresses warning message for slow for Groebner basis.

    ###############################################
    # (un)comment usecases for this package below #
    ###############################################

    # for case in ['087', '287', '365', '265s', '265t', '443', '243ss', '243st']:
    #    usecase__invariant_quadratic_forms( case )
    # usecase__toric_invariant_celestials()
    usecase__horn_and_spindle_cyclides()
    # usecase__complex_classification() # takes some time

    ###############################################

    mt.stop_timer()
    print
    print( 'The End' )

