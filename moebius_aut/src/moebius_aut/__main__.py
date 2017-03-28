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


def usecase__invariant_quadratic_forms( case = 2 ):
    '''
    OUTPUT:
        - Let G be a subgroup of Aut(P^1xP^1)
          Compute the vectors space of G-invariant quadratic forms in
          the ideal of the double Segre surface.
          For "DSegre.change_basis()" for the specification
          of involutions.
          
            case==1: 
                G = Aut(P^1xP^1)
                involution = 'identity' 
        
            case==2:
                G = toric automorphisms of P^1xP^1
                involution = 'rotate'
                  
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
    if case == 1:
        c_lst_lst = [t + e, q + e, s + e, e + t, e + q, e + s ]
        involution = 'identity'
        info = 'SL2+SL2'

    elif case == 2:
        c_lst_lst = [ s + e, e + s ]
        involution = 'rotate'
        info = 'toric'

    else:
        raise ValueError( 'Unknown case: ', case )

    #
    # compute vector space of invariant quadratic forms
    #
    iq_lst = DSegre.get_invariant_qf( c_lst_lst )
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
    mt.p( 'info about subgroup G of Aut(P^1xP^1):', info )
    mt.p( 'involution:', involution )
    mt.p( 'vector space of G-invariant quadratic forms:' )
    for iq in iq_lst:
        mt.p( '\t', iq )
    mt.p( 'signatures of random invariant quadratic forms:' )
    mt.p( '\t', sig_lst )
    mt.p( 'signatures of the form [1,n+1]:' )
    for sig in sig_lst:
        if 1 in sig:
            mt.p( '\t', sig )


def usecase__toric_invariant_celestials():
    '''
    OUTPUT:
        - Compute celestials that are invariant under 
          toric automorphisms
    '''

    R = PolynomialRing( QQ, ['x' + str( i ) for i in range( 9 )] )
    x0, x1, x2, x3, x4, x5, x6, x7, x8 = R.gens()
    lst = sage_eval( str( DSegre.get_ideal_lst() ), R.gens_dict() )
    J = R.ideal( lst )

    elim_lst = [1, 2]
    xelim_lst = [sage_eval( 'x' + str( i ), R.gens_dict() ) for i in elim_lst]
    J = J.elimination_ideal( xelim_lst )
    R2 = PolynomialRing( QQ, ['x' + str( i ) for i in range( 9 ) if i not in elim_lst] )
    J2 = R2.ideal( sage_eval( str( J.gens() ), R2.gens_dict() ) )
    for gen in J2.gens():
        mt.p( gen )
    hpol = J2.hilbert_polynomial()
    mt.p( hpol )



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

                J = DSegre.get_invariant_ideal( c_lst_lst )
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
    mt.filter( '__main__.py' )
    set_verbose( -1 )  # surpresses warning message for slow for Groebner basis.

    # usecase__invariant_quadratic_forms()
    usecase__toric_invariant_celestials()

    # usecase__complex_classification()



    mt.stop_timer()
    print
    print( 'The End' )

