'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''

from moebius_aut.class_ma_ring import ring
from moebius_aut.class_ma_ring import MARing

from moebius_aut.sage_interface import sage_matrix


class TestMARing:
    '''
    We denote b:=1/a in method naming.
    '''

    def test__ring__I2ba( self ):
        assert ring( 'I^2*(1/a)*a' ) == -1


    def test__get_sig__14( self ):
        pol = MARing.ring( '-x0^2+x1^2+x2^2+x3^2+x4^2' )
        sig = MARing.get_sig( pol )
        print( pol )
        print( sig )
        assert sig == [1, 4]


    def test__get_sig__13( self ):
        pol = MARing.ring( 'x4^2-x6*x7 + x2^2-x6*x8 - x2*x4+x0*x6' )
        sig = MARing.get_sig( pol )
        print( pol )
        print( sig )
        assert sig == [1, 3]

    def test__get_rand_sigs__( self ):

        pol_lst = []
        pol_lst += [ring( 'x4^2 - x6^2 - x7^2 ' )]
        pol_lst += [ring( 'x0^2 - x3*x4 ' )]

        sig_lst = MARing.get_rand_sigs( pol_lst, 1 )
        print( sig_lst )
        assert [ 1, 2 ] in sig_lst


    def test__diff_mat__a00b( self ):
        mat = ring( 'matrix([(a,0),(0,1/a)])' )
        a = ring( 'a' )
        dmat = MARing.diff_mat( mat, a )
        print( mat )
        print( dmat )
        assert dmat == sage_matrix( MARing.FF, [( 1, 0 ), ( 0, -1 / a ** 2 )] )


    def test__solve( self ):

        q = MARing.q()
        pol_lst = [ q[0] + 2 * q[1] + 3 * q[2] - 4 ]
        var_lst = [ q[i] for i in [0, 1, 2]]

        sol_dct = MARing.solve( pol_lst, var_lst )
        print( sol_dct )
        assert 'r0' in str( sol_dct )

        # We call a second time so that sage internally
        # introduces new variable names. We undo this
        # in MARing.solve(), since we work in a
        # polynomial ring where we want to minimize
        # the number of variables.
        #
        sol_dct = MARing.solve( pol_lst, var_lst )
        print( sol_dct )
        assert 'r0' in str( sol_dct )


    def test__solve__empty_list( self ):

        q = MARing.q()
        pol_lst = []
        var_lst = [ q[i] for i in [0, 1, 2]]

        sol_dct = MARing.solve( pol_lst, var_lst )
        print( sol_dct )
        assert 'r0' in str( sol_dct )


    def test__replace_conj_pairs( self ):

        q_lst = ring( '[ x0-I*x1, x0+I*x1, x0+x2, x3*I, x3+I*x4]' )
        new_lst = MARing.replace_conj_pairs( q_lst )
        chk_lst = sorted( ring( '[ x0, x1, x0+x2, x3, x3+I*x4]' ) )

        print( q_lst )
        print( new_lst )
        print( chk_lst )
        assert chk_lst == new_lst


