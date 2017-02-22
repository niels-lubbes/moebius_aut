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


    def test__diff_mat__a00b( self ):
        mat = ring( 'matrix([(a,0),(0,1/a)])' )
        a = ring( 'a' )
        dmat = MARing.diff_mat( mat, a )
        print( mat )
        print( dmat )
        assert dmat == matrix( MARing.FF, [( 1, 0 ), ( 0, -1 / a ** 2 )] )


