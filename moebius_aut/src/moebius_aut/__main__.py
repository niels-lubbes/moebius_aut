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
                J = DSegre.change_basis( J )

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
    usecase__double_segre()

    mt.stop_timer()
    print
    print( 'The End' )

