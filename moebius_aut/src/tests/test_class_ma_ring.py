'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''
from sage.all import *
from moebius_aut.class_ma_ring import ring
from moebius_aut.class_ma_ring import MARing


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


    def test__diff_mat__a00b( self ):
        mat = ring( 'matrix([(a,0),(0,1/a)])' )
        a = ring( 'a' )
        dmat = MARing.diff_mat( mat, a )
        print( mat )
        print( dmat )
        assert dmat == matrix( MARing.FF, [( 1, 0 ), ( 0, -1 / a ** 2 )] )


    def test__solve( self ):

        q = MARing.q()
        pol_lst = [ q[0] + 2 * q[1] + 3 * q[2] - 4 ]
        var_lst = [ q[i] for i in [0, 1, 2]]

        sol_dct = MARing.solve( pol_lst, var_lst )
        print( sol_dct )
        assert 'r0' in str( sol_dct )

        sol_dct = MARing.solve( pol_lst, var_lst )
        print( sol_dct )
        assert 'r0' in str( sol_dct )


    def test__replace_conj_pairs( self ):

        q_lst = ring( '[ x0-I*x1, x0+I*x1, x0+x2, x3*I]' )
        new_lst = MARing.replace_conj_pairs( q_lst )
        chk_lst = sorted( ring( '[ x0, x1, x0+x2, x3*I]' ) )

        print( q_lst )
        print( new_lst )
        print( chk_lst )
        assert chk_lst == new_lst


