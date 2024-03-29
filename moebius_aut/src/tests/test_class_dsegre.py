'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''

from moebius_aut.sage_interface import sage_diagonal_matrix
from moebius_aut.sage_interface import sage_matrix
from moebius_aut.sage_interface import sage_vector
from moebius_aut.sage_interface import sage_var
from moebius_aut.sage_interface import sage_sin
from moebius_aut.sage_interface import sage_cos
from moebius_aut.sage_interface import sage_ZZ
from moebius_aut.sage_interface import sage__eval
from moebius_aut.sage_interface import sage_diff

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


    def test__get_aut_P8__moduli( self ):

        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )
        M = DSegre.get_aut_P8( [a, b, c, d, e, f, g, h] )

        print( M )

        x = ring( '[x0, x1, x2, x3, x4, x5, x6, x7, x8]' )
        y = ring( '[y0, y1, y2, y3, y4, y5, y6, y7, y8]' )
        My = M * sage_vector( y )
        dct = {x[i]:My[i] for i in range( 9 )}

        pol = ring( 'x0^2-x1*x2' )
        print( pol.subs( dct ) )

        # We may use this to check the dimension of the moduli space
        # of invariant quadratic forms, since coefficients of some of the
        # terms are necessarily 0.
        #
        chk = ''
        chk += '(b^2*c^2*f^2*g^2 - 2*a*b*c*d*f^2*g^2 + a^2*d^2*f^2*g^2 + 2*b^2*c^2*e*f*g*h - 4*a*b*c*d*e*f*g*h + 2*a^2*d^2*e*f*g*h + b^2*c^2*e^2*h^2 - 2*a*b*c*d*e^2*h^2 + a^2*d^2*e^2*h^2 )*y0^2  + '
        chk += '(-b^2*c^2*f^2*g^2 + 2*a*b*c*d*f^2*g^2 - a^2*d^2*f^2*g^2 - 2*b^2*c^2*e*f*g*h + 4*a*b*c*d*e*f*g*h - 2*a^2*d^2*e*f*g*h - b^2*c^2*e^2*h^2 + 2*a*b*c*d*e^2*h^2 - a^2*d^2*e^2*h^2)*y1*y2 + '
        chk += '(2*b^2*c^2*e*f*g*h - 4*a*b*c*d*e*f*g*h + 2*a^2*d^2*e*f*g*h                                                                                                                 )*y3*y4 + '
        chk += '(-b^2*c^2*e*f*g*h + 2*a*b*c*d*e*f*g*h - a^2*d^2*e*f*g*h                                                                                                                    )*y5*y6 + '
        chk += '(-b^2*c^2*e*f*g*h + 2*a*b*c*d*e*f*g*h - a^2*d^2*e*f*g*h                                                                                                                    )*y7*y8 + '
        chk += '(b^2*c^2*e^2*g^2 - 2*a*b*c*d*e^2*g^2 + a^2*d^2*e^2*g^2                                                                                                                     )*y3^2  + '
        chk += '(b^2*c^2*f^2*h^2 - 2*a*b*c*d*f^2*h^2 + a^2*d^2*f^2*h^2                                                                                                                     )*y4^2  + '
        chk += '(2*b^2*c^2*e*f*g^2 - 4*a*b*c*d*e*f*g^2 + 2*a^2*d^2*e*f*g^2 + 2*b^2*c^2*e^2*g*h - 4*a*b*c*d*e^2*g*h + 2*a^2*d^2*e^2*g*h                                                     )*y0*y3 + '
        chk += '(2*b^2*c^2*f^2*g*h - 4*a*b*c*d*f^2*g*h + 2*a^2*d^2*f^2*g*h + 2*b^2*c^2*e*f*h^2 - 4*a*b*c*d*e*f*h^2 + 2*a^2*d^2*e*f*h^2                                                     )*y0*y4 + '
        chk += '(-b^2*c^2*e*f*g^2 + 2*a*b*c*d*e*f*g^2 - a^2*d^2*e*f*g^2 - b^2*c^2*e^2*g*h + 2*a*b*c*d*e^2*g*h - a^2*d^2*e^2*g*h                                                            )*y2*y5 + '
        chk += '(-b^2*c^2*f^2*g*h + 2*a*b*c*d*f^2*g*h - a^2*d^2*f^2*g*h - b^2*c^2*e*f*h^2 + 2*a*b*c*d*e*f*h^2 - a^2*d^2*e*f*h^2                                                            )*y1*y6 + '
        chk += '(-b^2*c^2*f^2*g*h + 2*a*b*c*d*f^2*g*h - a^2*d^2*f^2*g*h - b^2*c^2*e*f*h^2 + 2*a*b*c*d*e*f*h^2 - a^2*d^2*e*f*h^2                                                            )*y2*y7 + '
        chk += '(-b^2*c^2*f^2*h^2 + 2*a*b*c*d*f^2*h^2 - a^2*d^2*f^2*h^2                                                                                                                    )*y6*y7 + '
        chk += '(-b^2*c^2*e*f*g^2 + 2*a*b*c*d*e*f*g^2 - a^2*d^2*e*f*g^2 - b^2*c^2*e^2*g*h + 2*a*b*c*d*e^2*g*h - a^2*d^2*e^2*g*h                                                            )*y1*y8 + '
        chk += '(-b^2*c^2*e^2*g^2 + 2*a*b*c*d*e^2*g^2 - a^2*d^2*e^2*g^2                                                                                                                    )*y5*y8   '

        assert pol.subs( dct ) == ring( chk )


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
        assert out == sage_diagonal_matrix( [1, t, 1 / t, 1, 1, t, 1 / t, t, 1 / t ] )


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
        qpol = list( sage_vector( x ).row() * qmat * sage_vector( x ).column() )[0][0]
        assert qpol == chk_qpol


    def test__get_invariant_q_lst( self ):
        k = ring( 'k' )

        ig_lst = DSegre.get_invariant_q_lst( [k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1] )
        print( ig_lst )
        assert ig_lst == ring( '[q8 + q12, -q9 - q13, q16, -q17, q18, -q19, q8 + q12, 4*q4, -q16, -q18, -q9 - q13, (-4)*q5, q19, q17, -q16, q19, q9, -q12, -q18, q17, -q8, q13, q16, -q8, (-2)*q4, -q17, q9, 2*q5, q18, -q12, (-2)*q4, -q19, q13, 2*q5]' )

        ig_lst = DSegre.get_invariant_q_lst( [1, k, 0, 1, 1, 0, 0, 1] )
        print( ig_lst )
        assert ig_lst == ring( '[2*q8 + 2*q12, 2*q4, q1 + q2 + q3, (-1/2)*q11 + 1/2*q15, 1/2*q10 + (-1/2)*q14, 2*q4, 1/2*q8 + 1/2*q12, (-1/2)*q18, (-1/2)*q16, q1 + q2 + q3, 1/2*q8 + 1/2*q12, q9 + q13, (-1/2)*q10 + 1/2*q14, 1/2*q11 + (-1/2)*q15, 1/2*q16, 1/2*q18, (-1/2)*q10 + 1/2*q14, -q8 - q12, (-1/2)*q1 - q2, -q4, 1/2*q11 + (-1/2)*q15, -q8 - q12, -q4, (-1/2)*q1 - q3, 1/2*q16, -q4, (-1/2)*q8, (-1/2)*q11 + 1/2*q15, (-1/2)*q18, (-1/2)*q1 - q2, (-1/2)*q8, (-1/2)*q9 + (-1/2)*q13, 1/2*q18, -q4, (-1/2)*q12, 1/2*q10 + (-1/2)*q14, (-1/2)*q16, (-1/2)*q1 - q3, (-1/2)*q9 + (-1/2)*q13, (-1/2)*q12]' )


    def test__get_invariant_qf__SO2xSO2( self ):
        k = ring( 'k' )
        c_lst_lst = []
        c_lst_lst += [[k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1]]
        c_lst_lst += [[1, 0, 0, 1, k + 1, 0, 0, 1 / ( k + 1 )]]

        iqf_lst = DSegre.get_invariant_qf( c_lst_lst )
        print( iqf_lst )
        assert iqf_lst == ring( '[x0^2 - x7*x8, x0^2 - x5*x6, x0^2 - x3*x4, x0^2 - x1*x2]' )


    def test__get_invariant_qf__5678_SO2xSO2( self ):
        k = ring( 'k' )
        c_lst_lst = []
        c_lst_lst += [[k + 1, 0, 0, 1 / ( k + 1 ), 1, 0, 0, 1]]
        c_lst_lst += [[1, 0, 0, 1, k + 1, 0, 0, 1 / ( k + 1 )]]

        iqf_lst = DSegre.get_invariant_qf( c_lst_lst, [5, 6, 7, 8] )
        print( iqf_lst )
        assert iqf_lst == ring( '[x0^2 - x1*x2, x0^2 - x3*x4]' )


    def test__get_c_lst_lst_lst( self ):
        k = ring( 'k' )
        for c_lst_lst in DSegre.get_c_lst_lst_lst():
            for c_lst in c_lst_lst:
                # evaluate at k=0
                n_lst = []
                for c in c_lst:
                    if c in sage_ZZ:
                        n_lst += [ c ]
                    else:
                        n_lst += [c.subs( {k:0} ) ]
                print( n_lst )
                assert n_lst == [1, 0, 0, 1, 1, 0, 0, 1]
                for c in c_lst:
                    assert c in MARing.FF


    def test__to_str( self ):

        s = '[ '
        for c_lst_lst in DSegre.get_c_lst_lst_lst():
             s += DSegre.to_str( c_lst_lst ) + ', '
        s = s[:-2] + ' ]'
        print( s )

        chk = '[ '
        chk += '< t1, q1, s1 >, < t1, s1 >, < t1 >, < s1 >, < r1 >, '
        chk += '< t1, q1, s1, t2, q2, s2 >, < t1, q1, s1, t2, s2 >, '
        chk += '< t1, q1, s1, t2 >, < t1, q1, s1, s2 >, < t1, q1, s1, r2 >, '

        chk += '< t1, s1, t2, s2 >, < t1, s1, t2 >, < t1, s1, s2 >, < t1, s1, r2 >, '
        chk += '< t1, t2 >, < t1, s2 >, < t1, r2 >, < s1, s2 >, < s1, r2 >, < r1, r2 >, '

        chk += '< t1+t2, g1+q2, s1+s2 >, < t1+t2, s1+s2 >, < t1+t2 >, < s1+s2 >, < r1+r2 >, '
        chk += '< t1+s2 >, < t1+r2 >, < s1+r2 >, '

        chk += '< s1+s2, t1, t2 >, < s1+t2, t1 >, < s1+s2, t1 >, < s1+r2, t1 > '
        chk += ']'

        assert s == chk


    def test__get_aut_P8__action_of_involution( self ):
        '''
        OUPUT:
            -   We check how the action of the involution J: P^8---->P^8 acts on 
                elements "c_lst" in Aut(P^1xP^1) by conjugation. See documentation 
                "DSegre.change_basis()" and "DSegre.get_aut_P8()" for our internal 
                implementation of J and "c_lst" respectively.              
        '''
        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )

        # left-right involution
        #
        M = DSegre.get_aut_P8( [a, b, c, d] + [e, f, g, h] )
        L = sage_matrix( [
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
        R = sage_matrix( [
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


    def test__get_aut_P8__rotation_matrix( self ):
        '''
        OUTPUT:
            - We verify that for the 1-parameter subgroup of rotations 
              with matrix
              
                  [ cos(k) -sin(k) ]
                  [ sin(k)  cos(k) ]
            
              it is for our Lie algebra methods sufficient to consider 
              1-parameter subgroups that have the same tangent vector
              at the identity. Note that for k=0 we need to get the identity
              and the determinant should be 1.   
        '''
        k = ring( 'k' )
        I = ring( 'I' )
        a, b, c, d, e, f, g, h = ring( 'a, b, c, d, e, f, g, h' )
        x = MARing.x()
        q = MARing.q()
        r = MARing.r()

        #
        # We consider the representation of the
        # following matrix into P^8
        #
        #   (  [ 1 -k ] ,  [ 1 0 ] )
        #   (  [ k  1 ]    [ 0 1 ] )
        #
        # We define A to be the tangent vector at the
        # identity of this representation
        #
        N = DSegre.get_aut_P8( [1, -k, k , 1] + [1, 0, 0, 1] )
        A = MARing.diff_mat( N, k ).subs( {k:0} )


        #
        # We consider the representation of the
        # following matrix into P^8
        #
        #   (  [ cos(k) -sin(k) ] ,  [ 1 0 ] )
        #   (  [ sin(k)  cos(k) ]    [ 0 1 ] )
        #
        # We define B to be the tangent vector at the
        # identity of this representation
        #
        M = DSegre.get_aut_P8( [a, b, c, d] + [e, f, g, h] )

        a, b, c, d, e, f, g, h = sage_var( 'a, b, c, d, e, f, g, h' )
        k = sage_var( 'k' )
        M = sage__eval( str( list( M ) ), {'a':a, 'b':b, 'c':c, 'd':d, 'e':e, 'f':f, 'g':g, 'h':h, 'k':k} )
        M = sage_matrix( M )
        M = M.subs( {a:sage_cos( k ), b:-sage_sin( k ), c:sage_sin( k ), d:sage_cos( k ), e:1, f:0, g:0, h:1} )

        # differentiate the entries of M wrt. k
        dmat = []
        for row in M:
            drow = []
            for col in row:
                drow += [ sage_diff( col, k ) ]
            dmat += [drow]
        M = sage_matrix( dmat )
        B = M.subs( {k:0} )

        assert str( A ) == str( B )



if __name__ == '__main__':
    # TestClassDSegre().test__get_aut_P8__action_of_involution()
    # TestClassDSegre().test__get_invariant_qf__5678_SO2xSO2()
    # TestClassDSegre().test__get_c_lst_lst_lst()
    # TestClassDSegre().test__to_str()
    TestClassDSegre().test__get_aut_P8__moduli()
