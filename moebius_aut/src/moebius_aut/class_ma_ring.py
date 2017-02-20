'''
Created on Nov 21, 2016

@author: niels
'''
from sage.all import *

def ring( s ):
    return MARing.ring( s )

class MARing:

    a = PolynomialRing( QQ, 'a' ).gens()[0]
    NF = NumberField( [a ** 2 + 1], 'I' )
    I = NF.gens()[0]

    va_lst = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'k']
    FF = FractionField( PolynomialRing( NF, va_lst ) )

    vx_lst = [ 'x' + str( i ) for i in range( 9 )]
    vy_lst = [ 'y' + str( i ) for i in range( 9 )]
    vz_lst = [ 'z' + str( i ) for i in range( 9 )]
    vs_lst = ['s', 't', 'u', 'w']
    vq_lst = [ 'q' + str( i ) for i in range( 20 ) ]
    vr_lst = [ 'r' + str( i ) for i in range( 20 ) ]

    R = PolynomialRing( FF, vx_lst + vy_lst + vz_lst + vs_lst + vq_lst + vr_lst )

    gens_dct = NF.gens_dict()
    gens_dct.update( FF.gens_dict() )
    gens_dct.update( R.gens_dict() )

    # R.inject_variables();NF.inject_variables();FF.inject_variables()

    @staticmethod
    def x():
        return MARing.ring( 'x0, x1, x2, x3, x4, x5, x6, x7, x8' )

    @staticmethod
    def y():
        return MARing.ring( 'y0, y1, y2, y3, y4, y5, y6, y7, y8' )

    @staticmethod
    def z():
        return MARing.ring( 'z0, z1, z2, z3, z4, z5, z6, z7, z8' )

    @staticmethod
    def q():
        return MARing.ring( 'q0,q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,q13,q14,q15,q16,q17,q18,q19' )

    @staticmethod
    def r():
        return MARing.ring( 'r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19' )

    @staticmethod
    def ring( expr ):
        return sage_eval( str( expr ), MARing.gens_dct )


    @staticmethod
    def random_int( val ):
        '''
        INPUT:
            - "val" -- An integer.
        OUTPUT:
            - A random element in the interval [-val,val]
        '''
        return int( ZZ.random_element( -val, val + 1 ) )

    @staticmethod
    def random_elt( lst ):
        '''
        INPUT:
            - "lst" -- A list.
        OUTPUT:
            - A random element in "lst".
        '''
        idx = int( ZZ.random_element( 0, len( lst ) ) )
        return lst[idx]
