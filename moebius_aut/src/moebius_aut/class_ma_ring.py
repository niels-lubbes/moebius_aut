'''
Created on Nov 21, 2016

@author: niels
'''
from sage.all import *

from class_ma_tools import MATools

def ring( s ):
    return MARing.ring( s )

class MARing:
    '''
    This class represents a polynomial ring over a fraction field:
    
        FF[s,t,u,w,x#,y#,z#,q#,r#]
    
    where # is in [0,19]
    
        FF := NF(a,b,c,d,e,f,g,h,k) 
    and
     
        NF := QQ(I) 
    
    with I^2==-1
    
    The r% variables are in order to cast solution dictionaries as 
    returned by the Sage "solve" method. 
    '''

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
        return MARing.ring( [ring( 'r' + str( i ) ) for i in range( 20 )] )


    @staticmethod
    def solve( pol_lst, var_lst ):
        '''
        INPUT:
            - "pol_lst" -- List of polynomials in "MARing.R" 
                           but not in the variables r0,r1,....
            - "var_lst" -- List of variables at most 20 of "MARing.R".
        
        OUTPUT:
            - A dictionary of parametrized solutions of the  
              zeroset of the ideal generated by the polynomials in 
              "pol_lst". The solutions are parametrized 
              in terms of r0,r1,... 
        
        METHOD:
            - We call Sage "solve" method, but reset the indeterminates
              to be r0,...,r19, otherwise for each new call new 
              r# variables are introduced. The "reset" and "restore"
              functions do not work. See also:
              <https://ask.sagemath.org/question/23719/how-to-reset-variables/>              
              The indeterminates of the "solve" method are either c or r. 
        '''
        # cast to Sage symbolic ring SR and
        # call solve.
        #
        spol_lst = [ SR( str( pol ) ) for pol in pol_lst ]
        svar_lst = [ SR( str( var ) ) for var in var_lst ]
        dct = solve( spol_lst, svar_lst, solution_dict = True )[0]
        MATools.p( 'dct =', dct )

        # construct a list with integers i such
        # ri is occuring in the "dct" dictionary.
        #
        end_char_lst = []
        s = str( dct )
        rn_lst = []
        i = 0
        while i < len( s ):

            if s[i] == 'r':
                i += 1
                n = ''
                while s[i] in [str( j ) for j in range( 10 )]:
                    n += s[i]
                    i += 1

                # store end character: eg. 'r10}' OR 'r10*' OR 'r10,' OR 'r10 '
                if s[i] not in end_char_lst:
                    end_char_lst += [ s[i] ]

                if int( n ) not in rn_lst:
                    rn_lst += [int( n )]
            i += 1

        MATools.p( 'end_char_lst =', end_char_lst )
        MATools.p( 'rn_lst =', rn_lst )

        # replace each occurence of r<i> with r<t>
        # where t=0,...,t19
        #
        t = 0
        for n in sorted( rn_lst ):

            #
            # note that we want to prevent for example:
            #     r1-->r0 but also r10-->r00, for this
            # reason we consider the character which indicates
            # the end of the r-variable string so that:
            #    'r1}'-->'r0}', 'r1,'-->'r0,', 'r10'--->'r10'
            #
            for end_char in end_char_lst:
                s = s.replace( 'r' + str( n ) + end_char, 'r' + str( t ) + end_char )
            t += 1


        # cast to solution dictionary with
        # keys and values in "MARing.R"
        #
        sol_dct = ring( s )
        MATools.p( 'sol_dct =', sol_dct )

        return sol_dct


    @staticmethod
    def ring( expr ):
        return sage_eval( str( expr ), MARing.gens_dct )


    @staticmethod
    def diff_mat( mat , var ):
        '''
        INPUT:
            - "mat" -- A matrix defined over "MARing".
            - "var" -- An indeterminate in "MARing". 
        OUTPUT:
            - Returns a matrix whose entries are
              differentiated wrt. "var".              
        '''

        dmat = []
        for row in mat:
            drow = []
            for col in row:
                drow += [ diff( col, var ) ]
            dmat += [drow]

        return matrix( dmat )


    @staticmethod
    def get_sig( pol ):
        '''
        INPUT:
            - "pol" -- A quadratic form in the subring QQ[x0,...,x9] of the "MARing".
        
        OUTPUT:
            - An ordered pair of integers which denotes the signature 
              of the 9x9 Gramm matrix of the quadratic form.
                            
        EXAMPLE:
            - "get_sig(MARing.ring('-x0^2+x1^2+x2^2+x3^2+x4^2'))==[1,4]"
              The return value [1,4] means that the matrix is conjugate 
              to matrix with diagonal either (1,-1,-1,-1,-1,0,0,0,0) 
              or (-1,1,1,1,1,0,0,0,0).
        '''

        x = MARing.x()

        M = invariant_theory.quadratic_form( pol, x ).as_QuadraticForm().matrix()
        M = matrix( QQ, M )
        D, V = M.eigenmatrix_right()

        # determine signature of quadric
        #
        num_neg = len( [ d for d in D.diagonal() if d < 0 ] )
        num_pos = len( [ d for d in D.diagonal() if d > 0 ] )

        return sorted( [ num_neg, num_pos ] )


    @staticmethod
    def get_rand_sigs( pol_lst, num_tries = 10 ):
        '''
        INPUT:
            - "pol_lst"   -- A list of polynomials in the ring "MARing.R".
            - "num_tries" -- A positive integer. 
        OUTPUT:
            - Returns a sorted list of pairs of integers that represent 
              signatures of random quadratic forms in the ideal generated 
              by "pol_lst". For each subset of the list of quadrics, this 
              method computes "num_tries" of random linear combination of 
              this subset of quadratic forms. 
        '''
        q_lst = [ pol for pol in pol_lst if pol.total_degree() == 2 ]

        sig_lst = []
        for s_lst in Subsets( len( q_lst ) ):
            for i in range( num_tries ):
                quad = sum( [  QQ.random_element() * q_lst[s - 1] for s in s_lst ] )
                if quad in QQ:
                    continue
                sig = MARing.get_sig( quad )
                if sig not in sig_lst:
                    sig_lst += [sig ]
                    MATools.p( 'sig =', sig, '\t\t quad =', quad )
                    MATools.p( '\t\t sub_lst', [q_lst[s - 1] for s in s_lst ] )

        return sorted( sig_lst )



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


    @staticmethod
    def replace_conj_pairs( q_lst ):
        '''
        INPUT:
            - "q_lst" -- A list of polynomials in ".R".                               
        OUTPUT:
            - If there is a pair of complex conjugate polynomials 
              (q1,q2) with q1 and q2 in "q_lst" then this pair is 
              replaced by ((q1+q2)/2,((q1-q2)*I)/(-2)).  
        '''
        I = ring( 'I' )

        c_lst = []
        r_lst = []
        for q in q_lst:
            if 'I' in str( q ):
                c_lst += [q]
            else:
                r_lst += [q]

        pair_lst = []
        for c1 in c_lst:
            for c2 in c_lst:
                A = ( c1 + c2 ) / 2
                B = ( c1 - c2 ) * I / 2
                if c1 in pair_lst and c2 in pair_lst:
                    continue
                elif 'I' in str( A ) or 'I' in str( B ):
                    continue
                else:
                    r_lst += [A, B]
                    pair_lst += [c1, c2]

        nc_lst = []
        for c in c_lst:
            if c not in pair_lst:
                nc_lst += [c]

        return sorted( r_lst + nc_lst )
