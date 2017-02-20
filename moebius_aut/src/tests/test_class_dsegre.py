'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Feb 20, 2017
@author: Niels Lubbes
'''
from sage.all import *

from moebius_aut.class_dsegre import DSegre
from moebius_aut.class_ma_ring import ring
# from moebius_aut.class_ma_ring import MARing

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


    def test__get_pmz_lst( self ):
        out = DSegre.get_pmz_lst()
        print( out )
        assert len( out ) == 9


    def test__get_aut_P8( self ):
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

        out = DSegre.get_aut_P8()
        print( list( out ) )
        print( out )
        assert False
        assert out == ring( chk_mat )





