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


def usecase__represent():
    '''
    OUTPUT:
        - Compute representations of 1-parameter subgroups of
          Aut(P^1xP^1) into End(P^8). 
    '''



    if False:
        w0 = ring( '2*x3*x4 - x6*x7 + x5*x8' )
        w1 = ring( 'x0*x4 - x2*x7 + x1*x8 ' )
        w2 = ring( 'x0*x4 + x2*x5 - x1*x6' )
        w3 = ring( 'x3^2 - x4^2 - x5*x7 - x6*x8' )
        w4 = ring( 'x0*x3 - x1*x7 - x2*x8' )
        w5 = ring( 'x0*x3 - x1*x5 - x2*x6' )
        w6 = ring( 'x0^2 - x7^2 - x8^2' )
        w7 = ring( 'x0^2 - x5^2 - x6^2 ' )
        w8 = ring( 'x0^2 - x3^2 - x4^2 ' )
        w9 = ring( 'x0^2 - x1^2 - x2^2' )

        mt.p( MARing.get_sig( w9 + w8 + w6 + w0 ) )
        sys.exit()

    if True:
        w0 = ring( 'x5*x6 - x7*x8' )
        w1 = ring( 'x1*x6 - x2*x7 ' )
        w2 = ring( 'x2*x5 - x1*x8 ' )
        w3 = ring( 'x4^2 - x6^2 - x7^2 ' )
        w4 = ring( 'x0*x4 - x2*x6 - x1*x7' )
        w5 = ring( 'x3^2 - x5^2 - x8^2 ' )
        w6 = ring( 'x0*x3 - x1*x5 - x2*x8' )
        w7 = ring( 'x0^2 - x5*x7 - x6*x8' )
        w8 = ring( 'x0^2 - x3*x4 ' )
        w9 = ring( 'x0^2 - x1^2 - x2^2' )

        # for w in [w0, w1, w2, w3, w4, w5, w6, w7, w8, w9]:
        #    mt.p( MARing.get_sig( w9 + w3 - w8 - w ) )

        # pol = w9 + w3 - w8
        pol = w3 + w5
        mt.p( pol, MARing.get_sig( pol ) )
        sys.exit()


    k = ring( 'k' )
    I = ring( 'I' )
    a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )
    U, V = k + 1, 1 / ( k + 1 )

    a, b, c, d = 1, 0, 0, 1
    A = matrix( MARing.R, [( a, b ), ( c, d )] )
    B = matrix( MARing.R, [( d, -b ), ( -c, a )] )
    M = matrix( MARing.R, [( U * I, 0 ), ( 0, -I * V )] )
    N = B * M * A
    c_lst = N.list() + [1, 0, 0, 1]

    # c_lst = [I * U, 0, 0, -I * V] + [1, 0, 0, 1]

    H = DSegre.get_aut_P8( c_lst )
    D = MARing.diff_mat( H, k ).subs( {k:0} )
    Q = DSegre.get_qmat()

    mt.p( 'M =\n' + str( M ) )
    mt.p( 'N = B*M*A =\n' + str( N ) )
    mt.p( 'H =\n' + str( H ) )
    mt.p( 'D =\n' + str( D ) )

    Z = D.T * Q + Q * D
    J = MARing.R.ideal( ring( Z.list() ) )
    gb_lst = list( J.groebner_basis() )
    mt.p( 'gb =' )
    for gb in gb_lst:
        mt.p( '\t', gb )


    iq_lst = DSegre.get_invariant_ideal( [c_lst] )
    iq_lst = DSegre.change_basis( iq_lst, 'leftright' )
    iq_lst = MARing.replace_conj_pairs( iq_lst )

    mt.p( 'iq_lst =' )
    for iq in iq_lst:
        mt.p( '\t', iq )
        abcd = False
        for char in ring( '[a, b, c, d]' ):
            if str( char ) in str( iq ):
                abcd = True
        if not abcd:
            mt.p( '\t\t', MARing.get_sig( iq ) )


    w0, w1, w2, w3, w4, w5, w6, w7, w8, w9 = iq_lst
    mt.p( MARing.get_sig( w9 ) )



def usecase__double_segre():
    '''
    OUTPUT:
        - Prints a classification of quadratic forms
          that contain some fixed double Segre surface
          in projective 8-space, such that the quadratic
          form is invariant under subgroup of Aut(S).
    '''
    # surpresses warning message for slow toy implementation
    # for Groebner basis.
    #
    set_verbose( -1 )

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

    usecase__represent()
    # usecase__double_segre()



    mt.stop_timer()
    print
    print( 'The End' )

