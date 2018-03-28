'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Apr 4, 2017
@author: Niels Lubbes
'''

from moebius_aut.sage_interface import sage_matrix
from moebius_aut.sage_interface import sage_QQ

from moebius_aut.class_veronese import Veronese

from moebius_aut.class_ma_ring import ring
from moebius_aut.class_ma_ring import MARing


class TestClassVeronese( object ):


    def test__get_pmz_lst( self ):

        id_lst = Veronese.get_ideal_lst()
        p = Veronese.get_pmz_lst()
        x = ring( 'x0,x1,x2,x3,x4,x5' )
        dct = {x[i]:p[i] for i in range( 5 + 1 )}

        for id in id_lst:
            ids = id.subs( dct )
            print( id, '--->', ids )
            assert ids == 0


    def test__change_basis( self ):

        l_lst = l0, l1, l2, l3, l4, l5 = Veronese.change_basis( Veronese.get_ideal_lst() )
        print( 'l_lst =', l_lst )
        #
        # l_lst = [
        # x1^2 - x4^2 - x5^2,
        # x0*x1 - x2^2 - x3^2,
        # x2^2 + (2*I)*x2*x3 - x3^2 - x0*x4 + ((-I))*x0*x5,
        # x2^2 + ((-2*I))*x2*x3 - x3^2 - x0*x4 + (I)*x0*x5,
        # x1*x2 + (I)*x1*x3 - x2*x4 + (I)*x3*x4 + ((-I))*x2*x5 - x3*x5,
        # x1*x2 + ((-I))*x1*x3 - x2*x4 + ((-I))*x3*x4 + (I)*x2*x5 - x3*x5
        # ]
        #

        m0 = l0
        m1 = l1
        m2 = ( l2 + l3 ) * ring( '1/2' )
        m3 = ( l2 - l3 ) * ring( '1/2*I' )
        m4 = ( l4 + l5 ) * ring( '1/2' )
        m5 = ( l4 - l5 ) * ring( '1/2*I' )
        m_lst = [m0, m1, m2, m3, m4, m5]
        print( 'm_lst =', m_lst )

        x0, x1 = ring( 'x0,x1' )
        dct = {x1:-x0 - x1}
        n_lst = n0, n1, n2, n3, n4, n5 = [ m.subs( dct ) for m in m_lst]
        print( 'n_lst =', n_lst )

        print( n0 + 2 * n1 )
        chk = ring( 'x1^2-x0^2-2*x2^2-2*x3^2-x4^2-x5^2' )

        assert n0 + 2 * n1 == chk


    def test__get_aut_P5__abcdefghk( self ):
        chk_mat = ''
        chk_mat += '['
        chk_mat += '( k^2 , 2*g*h     , 2*g*k        , 2*h*k     , g^2 , h^2 ),'
        chk_mat += '( c*f , b*d + a*e , c*d + a*f    , c*e + b*f , a*d , b*e ),'
        chk_mat += '( c*k , b*g + a*h , c*g + a*k    , c*h + b*k , a*g , b*h ),'
        chk_mat += '( f*k , e*g + d*h , f*g + d*k    , f*h + e*k , d*g , e*h ),'
        chk_mat += '( c^2 , 2*a*b     , 2*a*c        , 2*b*c     , a^2 , b^2 ),'
        chk_mat += '( f^2 , 2*d*e     , 2*d*f        , 2*e*f     , d^2 , e^2 )'
        chk_mat += ']'

        chk_mat = 'matrix(' + chk_mat + ')'

        a, b, c, d, e, f, g, h, k = ring( 'a, b, c, d, e, f, g, h, k' )
        out = Veronese.get_aut_P5( [a, b, c, d, e, f, g, h, k] )
        assert out == ring( chk_mat )


    def test__get_c_lst_lst_dct( self ):

        k = ring( 'k' )
        for key in Veronese.get_c_lst_lst_dct():
            for c_lst in Veronese.get_c_lst_lst_dct()[key]:
                print( 'key   =', key )
                print( 'c_lst =', c_lst )
                print( 'det   =', sage_matrix( 3, 3, c_lst ).det() )

                n_lst = []
                for c in c_lst:
                    if c in sage_QQ:
                        n_lst += [c]
                    else:
                        n_lst += [c.subs( {k:0} )]
                print( 'n_lst =', n_lst )
                print

                if key != 'SO3(R)':
                    assert sage_matrix( 3, 3, c_lst ).det() == 1
                assert n_lst == [1, 0, 0, 0, 1, 0, 0, 0, 1]


    def test__get_invariant_q_lst__SO2( self ):
        h1, h2, a1, a2, a3, b1, b2, b3 = Veronese.get_c_lst_lst_dct()['SL3(C)']

        q_lst = Veronese.get_invariant_q_lst( h1 )
        out = sorted( list( set( q_lst ) ) )
        print( 'q_lst =', q_lst )
        print( 'out   =', out )

        q0, q1, q2, q3, q4, q5 = MARing.q()[:6]

        assert out == ring( '[-1/2*q5, 1/2*q5, -1/2*q4, 1/2*q4, -2*q3, q3, -q2, 2*q2]' )


    def test__get_invariant_qf__SL3( self ):
        sl3_lst = Veronese.get_c_lst_lst_dct()['SL3(C)']
        iqf_lst = Veronese.get_invariant_qf( sl3_lst )
        print( iqf_lst )
        assert iqf_lst == []


    def test__get_invariant_qf__SO3( self ):
        r1, r2, r3 = Veronese.get_c_lst_lst_dct()['SO3(R)']
        iqf_lst = Veronese.get_invariant_qf( [r1, r2, r3] )
        print( iqf_lst )
        assert iqf_lst != []


if __name__ == '__main__':
    # TestClassVeronese().test__get_c_lst_SL3()
    # TestClassVeronese().test__get_invariant_qf__SO3()
    pass



