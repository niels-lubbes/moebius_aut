'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''
from sage.all import *

from moebius_aut.class_dsegre import DSegre
from moebius_aut.class_ma_ring import ring
from moebius_aut.class_ma_ring import MARing


class TestClassDSegre:


    def test__get_ideal_lst__56_y( self ):
        chk_str = '['
        chk_str += 'y0^2 - y1*y2, y0^2 - y3*y4, y0^2 - y7*y8,'
        chk_str += 'y0*y1 - y3*y7, y0*y2 - y4*y8, y0*y3 - y1*y8,'
        chk_str += 'y0*y4 - y2*y7, -y1*y4 + y0*y7, -y2*y3 + y0*y8'
        chk_str += ']'

        out = DSegre.get_ideal_lst( [5, 6], 'y' )
        print( out )
        assert out == ring( chk_str )


    def test__change_basis_identity( self ):
        qf_lst = DSegre.get_ideal_lst()
        nqf_lst = DSegre.change_basis( qf_lst, 'identity' )

        assert qf_lst == nqf_lst


    def test__change_basis_leftright( self ):
        qf_lst = DSegre.get_ideal_lst()
        nqf_lst = DSegre.change_basis( qf_lst, 'leftright' )

        chk_lst = '['
        chk_lst += 'x0^2 - x1^2 - x2^2,'
        chk_lst += 'x0^2 - x3*x4,'
        chk_lst += 'x0^2 + I*x5*x6 - x5*x7 - x6*x8 - I*x7*x8,'
        chk_lst += 'x0^2 - I*x5*x6 - x5*x7 - x6*x8 + I*x7*x8,'
        chk_lst += 'x1^2 + 2*I*x1*x2 - x2^2 - I*x5*x6 - x5*x7 + x6*x8 - I*x7*x8,'
        chk_lst += 'x1^2 - 2*I*x1*x2 - x2^2 + I*x5*x6 - x5*x7 + x6*x8 + I*x7*x8,'
        chk_lst += 'x3^2 - x5^2 - x8^2,'
        chk_lst += 'x4^2 - x6^2 - x7^2,'
        chk_lst += 'x0*x1 + I*x0*x2 - x4*x5 - I*x4*x8,'
        chk_lst += 'x0*x1 - I*x0*x2 + I*x3*x6 - x3*x7,'
        chk_lst += 'x0*x3 - x1*x5 + I*x2*x5 - I*x1*x8 - x2*x8,'
        chk_lst += 'x0*x4 + I*x1*x6 - x2*x6 - x1*x7 - I*x2*x7,'
        chk_lst += 'x0*x1 + I*x0*x2 - I*x3*x6 - x3*x7,'
        chk_lst += 'x0*x1 - I*x0*x2 - x4*x5 + I*x4*x8,'
        chk_lst += 'x0*x3 - x1*x5 - I*x2*x5 + I*x1*x8 - x2*x8,'
        chk_lst += 'x0*x4 - I*x1*x6 - x2*x6 - x1*x7 + I*x2*x7,'
        chk_lst += '-x1*x3 - I*x2*x3 + x0*x5 + I*x0*x8,'
        chk_lst += '-x1*x4 + I*x2*x4 - I*x0*x6 + x0*x7,'
        chk_lst += '-x1*x4 - I*x2*x4 + I*x0*x6 + x0*x7,'
        chk_lst += '-x1*x3 + I*x2*x3 + x0*x5 - I*x0*x8'
        chk_lst += ']'

        for nqf in nqf_lst:
            print( nqf )
        assert nqf_lst == ring( chk_lst )


    def test__change_basis_rotate( self ):
        qf_lst = DSegre.get_ideal_lst()
        nqf_lst = DSegre.change_basis( qf_lst, 'rotate' )

        chk_lst = '['
        chk_lst += 'x0^2 - x1^2 - x2^2,'
        chk_lst += 'x0^2 - x3^2 - x4^2,'
        chk_lst += 'x0^2 - x5^2 - x6^2,'
        chk_lst += 'x0^2 - x7^2 - x8^2,'
        chk_lst += 'x1^2 + 2*I*x1*x2 - x2^2 - x5*x7 - I*x6*x7 - I*x5*x8 + x6*x8,'
        chk_lst += 'x1^2 - 2*I*x1*x2 - x2^2 - x5*x7 + I*x6*x7 + I*x5*x8 + x6*x8,'
        chk_lst += 'x3^2 + 2*I*x3*x4 - x4^2 - x5*x7 - I*x6*x7 + I*x5*x8 - x6*x8,'
        chk_lst += 'x3^2 - 2*I*x3*x4 - x4^2 - x5*x7 + I*x6*x7 - I*x5*x8 - x6*x8,'
        chk_lst += 'x0*x1 + I*x0*x2 - x3*x5 + I*x4*x5 - I*x3*x6 - x4*x6,'
        chk_lst += 'x0*x1 - I*x0*x2 - x3*x5 - I*x4*x5 + I*x3*x6 - x4*x6,'
        chk_lst += 'x0*x3 + I*x0*x4 - x1*x5 + I*x2*x5 - I*x1*x6 - x2*x6,'
        chk_lst += 'x0*x3 - I*x0*x4 - x1*x5 - I*x2*x5 + I*x1*x6 - x2*x6,'
        chk_lst += 'x0*x1 + I*x0*x2 - x3*x7 - I*x4*x7 - I*x3*x8 + x4*x8,'
        chk_lst += 'x0*x1 - I*x0*x2 - x3*x7 + I*x4*x7 + I*x3*x8 + x4*x8,'
        chk_lst += 'x0*x3 + I*x0*x4 - x1*x7 - I*x2*x7 + I*x1*x8 - x2*x8,'
        chk_lst += 'x0*x3 - I*x0*x4 - x1*x7 + I*x2*x7 - I*x1*x8 - x2*x8,'
        chk_lst += '-x1*x3 - I*x2*x3 - I*x1*x4 + x2*x4 + x0*x5 + I*x0*x6,'
        chk_lst += '-x1*x3 + I*x2*x3 + I*x1*x4 + x2*x4 + x0*x5 - I*x0*x6,'
        chk_lst += '-x1*x3 - I*x2*x3 + I*x1*x4 - x2*x4 + x0*x7 + I*x0*x8,'
        chk_lst += '-x1*x3 + I*x2*x3 - I*x1*x4 - x2*x4 + x0*x7 - I*x0*x8'
        chk_lst += ']'

        assert nqf_lst == ring( chk_lst )


    def test__get_pmz_lst( self ):
        out = DSegre.get_pmz_lst()
        print( out )
        assert len( out ) == 9


    def test__get_aut_P8__abcdefgh( self ):
        chk_mat = ''
        chk_mat += '['
        chk_mat += '(b*c*f*g + a*d*f*g + b*c*e*h + a*d*e*h, a*c*f*g + a*c*e*h, b*d*f*g + b*d*e*h, b*c*e*g + a*d*e*g, b*c*f*h + a*d*f*h, a*c*e*g, b*d*f*h, a*c*f*h, b*d*e*g),'
        chk_mat += '(2*a*b*f*g + 2*a*b*e*h, a^2*f*g + a^2*e*h, b^2*f*g + b^2*e*h, 2*a*b*e*g, 2*a*b*f*h, a^2*e*g, b^2*f*h, a^2*f*h, b^2*e*g),'
        chk_mat += '(2*c*d*f*g + 2*c*d*e*h, c^2*f*g + c^2*e*h, d^2*f*g + d^2*e*h, 2*c*d*e*g, 2*c*d*f*h, c^2*e*g, d^2*f*h, c^2*f*h, d^2*e*g),'
        chk_mat += '(2*b*c*e*f + 2*a*d*e*f, 2*a*c*e*f, 2*b*d*e*f, b*c*e^2 + a*d*e^2, b*c*f^2 + a*d*f^2, a*c*e^2, b*d*f^2, a*c*f^2, b*d*e^2),'
        chk_mat += '(2*b*c*g*h + 2*a*d*g*h, 2*a*c*g*h, 2*b*d*g*h, b*c*g^2 + a*d*g^2, b*c*h^2 + a*d*h^2, a*c*g^2, b*d*h^2, a*c*h^2, b*d*g^2),'
        chk_mat += '(4*a*b*e*f, 2*a^2*e*f, 2*b^2*e*f, 2*a*b*e^2, 2*a*b*f^2, a^2*e^2, b^2*f^2, a^2*f^2, b^2*e^2),'
        chk_mat += '(4*c*d*g*h, 2*c^2*g*h, 2*d^2*g*h, 2*c*d*g^2, 2*c*d*h^2, c^2*g^2, d^2*h^2, c^2*h^2, d^2*g^2),'
        chk_mat += '(4*a*b*g*h, 2*a^2*g*h, 2*b^2*g*h, 2*a*b*g^2, 2*a*b*h^2, a^2*g^2, b^2*h^2, a^2*h^2, b^2*g^2),'
        chk_mat += '(4*c*d*e*f, 2*c^2*e*f, 2*d^2*e*f, 2*c*d*e^2, 2*c*d*f^2, c^2*e^2, d^2*f^2, c^2*f^2, d^2*e^2)'
        chk_mat += ']'

        chk_mat = 'matrix(' + chk_mat + ')'

        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )
        out = DSegre.get_aut_P8( [a, b, c, d, e, f, g, h] )
        print( list( out ) )
        print( out )
        assert out == ring( chk_mat )


    def test__get_aut_P8__SO2xID( self ):
        chk_mat = ''
        chk_mat += '['
        chk_mat += '(1, 0, 0, 0, 0, 0, 0, 0, 0),'
        chk_mat += '(0, k^2 + 2*k + 1, 0, 0, 0, 0, 0, 0, 0),'
        chk_mat += '(0, 0, 1/(k^2 + 2*k + 1), 0, 0, 0, 0, 0, 0),'
        chk_mat += '(0, 0, 0, 1, 0, 0, 0, 0, 0),'
        chk_mat += '(0, 0, 0, 0, 1, 0, 0, 0, 0),'
        chk_mat += '(0, 0, 0, 0, 0, k^2 + 2*k + 1, 0, 0, 0),'
        chk_mat += '(0, 0, 0, 0, 0, 0, 1/(k^2 + 2*k + 1), 0, 0),'
        chk_mat += '(0, 0, 0, 0, 0, 0, 0, k^2 + 2*k + 1, 0),'
        chk_mat += '(0, 0, 0, 0, 0, 0, 0, 0, 1/(k^2 + 2*k + 1))'
        chk_mat += ']'
        chk_mat = 'matrix(' + chk_mat + ')'

        k = ring( 'k' )
        out = DSegre.get_aut_P8( [k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1] )
        print( list( out ) )
        print( out )
        t = ( k + 1 ) ** 2
        assert out == diagonal_matrix( [1, t, 1 / t, 1, 1, t, 1 / t, t, 1 / t ] )


    def test__qmat( self ):
        chk_mat = ''
        chk_mat += '['
        chk_mat += '(q0 + q1 + q2 + q3, 1/2*q8 + 1/2*q12, 1/2*q9 + 1/2*q13, 1/2*q10 + 1/2*q14, 1/2*q11 + 1/2*q15, 1/2*q16, 1/2*q17, 1/2*q18, 1/2*q19),'
        chk_mat += '(1/2*q8 + 1/2*q12, q4, (-1/2)*q0, (-1/2)*q16, (-1/2)*q18, 0, (-1/2)*q11, 0, (-1/2)*q14),'
        chk_mat += '(1/2*q9 + 1/2*q13, (-1/2)*q0, q5, (-1/2)*q19, (-1/2)*q17, (-1/2)*q10, 0, (-1/2)*q15, 0),'
        chk_mat += '(1/2*q10 + 1/2*q14, (-1/2)*q16, (-1/2)*q19, q6, (-1/2)*q1, 0, (-1/2)*q9, (-1/2)*q12, 0),'
        chk_mat += '(1/2*q11 + 1/2*q15, (-1/2)*q18, (-1/2)*q17, (-1/2)*q1, q7, (-1/2)*q8, 0, 0, (-1/2)*q13),'
        chk_mat += '(1/2*q16, 0, (-1/2)*q10, 0, (-1/2)*q8, 0, (-1/2)*q2, (-1/2)*q4, (-1/2)*q6),'
        chk_mat += '(1/2*q17, (-1/2)*q11, 0, (-1/2)*q9, 0, (-1/2)*q2, 0, (-1/2)*q7, (-1/2)*q5),'
        chk_mat += '(1/2*q18, 0, (-1/2)*q15, (-1/2)*q12, 0, (-1/2)*q4, (-1/2)*q7, 0, (-1/2)*q3),'
        chk_mat += '(1/2*q19, (-1/2)*q14, 0, 0, (-1/2)*q13, (-1/2)*q6, (-1/2)*q5, (-1/2)*q3, 0)'
        chk_mat += ']'
        chk_mat = 'matrix(' + chk_mat + ')'

        qmat = DSegre.get_qmat()
        for row in qmat:
            print( row )
        assert qmat == ring( chk_mat )


    def test__qmat__symmetric( self ):
        qmat = DSegre.get_qmat()
        print( qmat )
        assert qmat == qmat.T


    def test__qmat__qpol( self ):
        x = MARing.x()
        q = MARing.q()

        g_lst = DSegre.get_ideal_lst()
        chk_qpol = 0
        for i in range( len( g_lst ) ):
            chk_qpol += q[i] * g_lst[i]

        qmat = DSegre.get_qmat()
        qpol = list( vector( x ).row() * qmat * vector( x ).column() )[0][0]
        assert qpol == chk_qpol


    def test__get_invariant_q_lst( self ):
        k = ring( 'k' )

        ig_lst = DSegre.get_invariant_q_lst( [k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1] )
        print( ig_lst )
        assert ig_lst == ring( '[q8 + q12, -q9 - q13, q16, -q17, q18, -q19, q8 + q12, 4*q4, -q16, -q18, -q9 - q13, (-4)*q5, q19, q17, -q16, q19, q9, -q12, -q18, q17, -q8, q13, q16, -q8, (-2)*q4, -q17, q9, 2*q5, q18, -q12, (-2)*q4, -q19, q13, 2*q5]' )

        ig_lst = DSegre.get_invariant_q_lst( [1, k, 0, 1, 1, 0, 0, 1] )
        print( ig_lst )
        assert ig_lst == ring( '[2*q8 + 2*q12, 2*q4, q1 + q2 + q3, (-1/2)*q11 + 1/2*q15, 1/2*q10 + (-1/2)*q14, 2*q4, 1/2*q8 + 1/2*q12, (-1/2)*q18, (-1/2)*q16, q1 + q2 + q3, 1/2*q8 + 1/2*q12, q9 + q13, (-1/2)*q10 + 1/2*q14, 1/2*q11 + (-1/2)*q15, 1/2*q16, 1/2*q18, (-1/2)*q10 + 1/2*q14, -q8 - q12, (-1/2)*q1 - q2, -q4, 1/2*q11 + (-1/2)*q15, -q8 - q12, -q4, (-1/2)*q1 - q3, 1/2*q16, -q4, (-1/2)*q8, (-1/2)*q11 + 1/2*q15, (-1/2)*q18, (-1/2)*q1 - q2, (-1/2)*q8, (-1/2)*q9 + (-1/2)*q13, 1/2*q18, -q4, (-1/2)*q12, 1/2*q10 + (-1/2)*q14, (-1/2)*q16, (-1/2)*q1 - q3, (-1/2)*q9 + (-1/2)*q13, (-1/2)*q12]' )


    def test__get_invariant_ideal__SO2xSO2( self ):
        k = ring( 'k' )
        c_lst_lst = []
        c_lst_lst += [[k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1]]
        c_lst_lst += [[1, 0, 0, 1, k + 1, 0, 0, 1 / ( k + 1 )]]

        iqf_lst = DSegre.get_invariant_ideal( c_lst_lst )
        print( iqf_lst )
        assert iqf_lst == ring( '[x0^2 - x7*x8, x0^2 - x5*x6, x0^2 - x3*x4, x0^2 - x1*x2]' )


    def test__get_c_lst_lst_dct( self ):
        k = ring( 'k' )
        c_lst_lst_dct = DSegre.get_c_lst_lst_dct()
        for key in c_lst_lst_dct:
            for c_lst_lst in c_lst_lst_dct[key]:
                assert len( c_lst_lst ) == key
                for c_lst in c_lst_lst:
                    assert len( c_lst ) == 8
                    # evaluate at k=0
                    n_lst = []
                    for c in c_lst:
                        if c in ZZ:
                            n_lst += [ c ]
                        else:
                            n_lst += [c.subs( {k:0} ) ]
                    print( n_lst )
                    assert n_lst == [1, 0, 0, 1, 1, 0, 0, 1]
                    for c in c_lst:
                        assert c in MARing.FF


    def test__to_str( self ):

        c_lst_lst_dct = DSegre.get_c_lst_lst_dct()
        s = '[ '
        for c_lst_lst in c_lst_lst_dct[2]:
             s += DSegre.to_str( c_lst_lst ) + ', '
        s = s[:-2] + ' ]'
        print( s )
        assert s == '[ < g1, g3 >, < g1, t1 >, < g1, t3 >, < g3, t3 >, < g1xt3, t1 >, < g3xt3, t1 >, < g1xt1, g3xt3 > ]'


    def test__get_aut_P8__action_of_involution( self ):

        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )

        # left-right involution
        #
        M = DSegre.get_aut_P8( [a, b, c, d] + [e, f, g, h] )
        L = matrix( [
            ( 1, 0, 0, 0, 0, 0, 0, 0, 0 ),
            ( 0, 0, 1, 0, 0, 0, 0, 0, 0 ),
            ( 0, 1, 0, 0, 0, 0, 0, 0, 0 ),
            ( 0, 0, 0, 1, 0, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 1, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 0, 0, 1 ),
            ( 0, 0, 0, 0, 0, 0, 0, 1, 0 ),
            ( 0, 0, 0, 0, 0, 0, 1, 0, 0 ),
            ( 0, 0, 0, 0, 0, 1, 0, 0, 0 )
            ] )
        assert L == ~L
        LML = ~L * M * L
        LML_chk = DSegre.get_aut_P8( [d, c, b, a] + [e, f, g, h] )
        assert LML_chk == LML


        # rotate
        #
        M = DSegre.get_aut_P8( [a, b, c, d] + [e, f, g, h] )
        R = matrix( [
            ( 1, 0, 0, 0, 0, 0, 0, 0, 0 ),
            ( 0, 0, 1, 0, 0, 0, 0, 0, 0 ),
            ( 0, 1, 0, 0, 0, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 1, 0, 0, 0, 0 ),
            ( 0, 0, 0, 1, 0, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 1, 0, 0 ),
            ( 0, 0, 0, 0, 0, 1, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 0, 0, 1 ),
            ( 0, 0, 0, 0, 0, 0, 0, 1, 0 )
            ] )
        assert R == ~R
        RMR = ~R * M * R
        RMR_chk = DSegre.get_aut_P8( [d, c, b, a] + [h, g, f, e] )
        assert RMR_chk == RMR



