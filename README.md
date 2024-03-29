# Moebius-aut 


## Introduction

Moebius-aut is a python library that computes G-invariant quadratic forms in the ideal of surfaces.
In this implementation we consider surfaces that are projections of the double Segre surface 
or Veronese surface. The group G is assumed to be a subgroup of SL(2)xSL(2) or SL(3).  

This library depends on [SageMath](https://SageMath.org) libraries.

## Installation

* Install Sage from [SageMath](https://SageMath.org).
We assume that `sage` is accessible from your commandline interface.

* Install the `moebius_aut` package: 
```    
sage -pip install moebius_aut
```    
If you do not have root access use the following command instead:
```    
sage -pip install --user moebius_aut
```    

* We advice to upgrade the `moebius_aut` package regularly:
```
sage -pip install --upgrade moebius_aut
```
 
* To execute some [usecases](https://github.com/niels-lubbes/moebius_aut/blob/master/moebius_aut/src/moebius_aut/__main__.py) type:
```    
sage -python -m moebius_aut
```

* For showing which files were installed 
or for uninstalling the `moebius_aut` package, 
use one of the following commands:
```
sage -pip show --files moebius_aut
sage -pip uninstall moebius_aut
```


## Examples

For running the examples below, either copy paste the code into the Sage interface or run them as a Python module:

    sage -python -m my_module_name.py

See [this file](https://github.com/niels-lubbes/moebius_aut/blob/master/moebius_aut/src/moebius_aut/__main__.py) 
for more example usecases. 
See the [source code](https://github.com/niels-lubbes/moebius_aut/blob/master/moebius_aut/src/moebius_aut)
the io-specification of each function.
The [test functions](https://github.com/niels-lubbes/moebius_aut/blob/master/moebius_aut/src/tests)
might be informative for how to call each function.


### Example 1: Computing G-invariant quadratic forms in the ideal of the double Segre surface

We show in this example how to compute quadratic forms in the ideal of the 
double Segre surface that are invariant under the group SL(2)xSL(2).
The double Segre surface (also known as the Veronese-Segre surface) is an embedding of 
P^1xP^1 into projective 8-space P^8. 

```python
from moebius_aut.class_dsegre import DSegre
from moebius_aut.class_ma_ring import MARing

t, q, s, r, e, T = DSegre.get_gens_sl2()
c_lst_lst = [t + e, q + e, s + e, e + t, e + q, e + s ] # generators of Lie algebra sl2+sl2
exc_idx_lst = []  # double Segre surface
involution = 'identity'

iq_lst = DSegre.get_invariant_qf( c_lst_lst, exc_idx_lst )
iq_lst = DSegre.change_basis( iq_lst, involution )
print( iq_lst )
print( MARing.get_rand_sigs( iq_lst, 1 ) )
```
Output:
    
    [(-2)*x0^2 + 2*x1*x2 + 2*x3*x4 - x5*x6 - x7*x8]
    [[4, 5]]

We may conclude from the output that  there is a hyperquadric Q in P^8 
of signature (4,5) such that Aut(P^1xP^1)=SL(2)xSL(2) is contained in Aut(Q). 
    
    
### Internal representation of projections of the double Segre surface and their real structures

In this section we explain the `exc_idx_lst` and `involution` parameters.

We fix the following monomial parametrization of the double Segre surface

    (s,u) |--> (1:s:s^{-1}:u:u^{-1}:s*u:s^{-1}*u^{-1}:s*u^{-1}:s^{-1}*u)
               =
               (x0:x1:x2:x3:x4:x5:x6:x7:x8)

We consider the exponents of the monomials as coordinates in a lattice.
Thus, for example, `1=s^0*t^0` corresponds to the coordinate (0,0) which corresponds to `x0`
and the monomial `s*u^{-1}` corresponds to the coordinate (1,-1) which corresponds to `x7`, etc.
We schematically denote this correspondence by the following square:

    [8 3 5]    
    [2 0 1]    
    [6 4 7]

The antiholomorphic involution defined by the real structure acts as unimodular involution
on the square:     

    identity:    ( a, b ) |--> ( a, b)
    leftright:   ( a, b ) |--> (-a, b)
    rotate:      ( a, b ) |--> (-a,-b)  
    diagonal:    ( a, b ) |--> ( b, a)
    
Moreover, the real structure acts as an antiholomorphic involution on P^8. 
We consider antiholomorphic involutions that are defined by one of the following maps 
composed with complex conjugation: 
    
    identity: (x0:...:x8) |--> (x0:...:x8) 
    
    leftright:  x3 |--> x3, x0 |--> x0, x4 |--> x4,                     
                x1 |--> x1 + I*x2, x2 |--> x1 - I*x2, 
                x5 |--> x5 + I*x8, x8 |--> x5 - I*x8,                                
                x7 |--> x7 + I*x6, x6 |--> x7 - I*x6 
                         
    rotate: x0 |--> x0,
            x1 |--> x1 + I*x2, x2 |--> x1 - I*x2, 
            x3 |--> x3 + I*x4, x4 |--> x3 - I*x4, 
            x5 |--> x5 + I*x6, x6 |--> x5 - I*x6, 
            x7 |--> x7 + I*x8, x8 |--> x7 - I*x8      
                    
    diagonal:   x0 |--> x0, x5 |--> x5, x6 |--> x6,
                x3 |--> x3 + I*x1,  x1 |--> x3 - I*x1,
                x8 |--> x8 + I*x7,  x7 |--> x8 - I*x7,
                x2 |--> x2 + I*x4,  x4 |--> x2 - I*x4  

For example we obtain a parametrization of a smooth Del Pezzo surface of degree 6, 
if we set `involution='rotate'` and compose the parametrization of the double Segre surface 
with the linear projection 

    (x0:x1:x2:x3:x4:x5:x6:x7:x8) |--> (x0:x1:x2:x3:x4:x7:x8) 

We schematically denote this as    
    
    -= del Pezzo surface of degree 6 =-
        [* *  ]   
        [* * *]   
        [  * *]
    exc_idx_lst = [5,6]   
    involution='rotate'

Below we give an overview of choices for normal toric surfaces that contain
at least circles through each point

    -= double Segre surface =-                 -= del Pezzo surface of degree 6 =-   
        [* * *]                                    [* *  ]                           
        [* * *]                                    [* * *]                           
        [* * *]                                    [  * *]                           
    exc_idx_lst = []                           exc_ idx_lst = [5,6]                   
    involution = 'identity'                    involution = 'rotate'                 

    -= weak del Pezzo surface of degree 6 =-   -= Veronese surface =-                       
        [  *  ]                                    [*    ]                                  
        [* * *]                                    [* *  ]                                  
        [* * *]                                    [* * *]                                  
    exc_idx_lst = [5,8]                        exc_idx_lst = [1,3,5]                        
    involution='leftright'                     involution in ['identity','diagonal']        

    -= ring cyclide =-                         -= spindle cyclide =-        
        [  *  ]                                    [  *  ]               
        [* * *]                                    [* * *]               
        [  *  ]                                    [  *  ]               
    exc_idx_lst = [5,6,7,8]                    exc_idx_lst = [5,6,7,8]   
    involution='rotate'                        involution='leftright'    

    -= horn cyclide =-                         -= 2-sphere =-               
        [  *  ]                                    [     ]              
        [  *  ]                                    [* *  ]              
        [* * *]                                    [* *  ]              
    exc_idx_lst = [1,2,5,8]                    exc_idx_lst = [1,3,5,7,8]
    involution='leftright'                     involution='diagonal'   


### Internal representation of Lie subalgebra's of the Lie algebra sl2+sl2.

In this section we explain the `c_lst_lst` parameter.

If the real structure acts as complex conjugation (thus we set `involution='identity'`), 
then the Lie algebra sl2 is generated by the following 2x2 matrices    

    t = [ 0 1 ]   q = [ 0 0 ]   s = [ 1  0 ]  
        [ 0 0 ]       [ 1 0 ]       [ 0 -1 ]    

We also consider elements in sl2 that generate subalgebra's
    
    e = [ 0 0 ]   r = [ 0 -1 ]
        [ 0 0 ]       [ 1  0 ]

We obtain elements in the direct sum sl2+sl2 via of generators in sl2
    
    t+e = ( [ 0 1 ] , [ 0 0 ] )   t+t = ( [ 0 1 ] , [ 0 1 ] )    
          ( [ 0 0 ]   [ 0 0 ] )         ( [ 0 0 ]   [ 0 0 ] )
    
The direct sum sl2+sl2 is generated by

    c_lst_lst = [t + e, q + e, s + e, e + t, e + q, e + s ]

The subalgebra so2+so2 is generated by

    c_lst_lst = [ r+e, e+r ]

The real structure may act different as complex conjugation. 
For example if `involution='rotate'`, then the action on sl2+sl2 is    

    ( [ a b ] , [ e f ] )   |-->    ( [ d c ] , [ h g ] )
    ( [ c d ]   [ g h ] )           ( [ b a ]   [ f e ] )

followed by complex conjugation.
In this case the subalgebra so2+so2 is generated by

    [ S+e, e+S ]    where   S = [ I  0 ]
                                [ 0 -I ]

However, since --- for our purposes --- we consider the matrices up to scalar multiplication, 
we have that the Lie algebra of so2+so2 is generated by

    c_lst_lst = [ s+e, e+s ]

If we set `involution='leftright'`, then the Lie subalgebra so2+se1 is generated by

    c_lst_lst = [ s+e, e+t ]

If we set `involution='identity'`, then the Lie subalgebra se1+se1 is generated by

    c_lst_lst = [ t+e, e+t ]

If we set `involution='leftright'`, then the Lie subalgebra so2+sx1 is generated by

    c_lst_lst = [ s+e, e+s ]

Notice that so2, se1 and sx1 are the Lie algebra's of the groups 
SO(2), SE(1) and SX(1) respectively.
Here SX(1) is the group defined by the automorphisms of the projective line that preserve
two real points.  


### Example 2: Computing G-invariant quadratic forms in the ideal of a projection of the double Segre surface

We show in this example how to compute quadratic forms in the ideal of a linear projection of the 
double Segre surface that are invariant under the group SO(2)xSO(2). Moreover, we assume that 
the real structure is non-trivial.

```python
from moebius_aut.class_dsegre import DSegre
from moebius_aut.class_ma_ring import MARing
 
t, q, s, r, e, T = DSegre.get_gens_sl2()
exc_idx_lst = [5, 6]  # del Pezzo surface of degree 6
c_lst_lst = [ s + e, e + s ] # so2+so2
involution = 'rotate'

iq_lst = DSegre.get_invariant_qf( c_lst_lst, exc_idx_lst )
print( iq_lst )
iq_lst = DSegre.change_basis( iq_lst, involution )
print( iq_lst )
sig_lst = MARing.get_rand_sigs( iq_lst, 10 )
print( sig_lst )
```
Output:
    
    [x0^2 - x7*x8, x0^2 - x3*x4, x0^2 - x1*x2] 
    [x0^2 - x7^2 - x8^2, x0^2 - x3^2 - x4^2, x0^2 - x1^2 - x2^2]
    [[1, 2], [1, 4], [1, 6], [2, 3], [2, 5], [3, 4]]

The output shows that there exists a del Pezzo surface of degree 6 in the projective 5-sphere S^5 
that is invariant under a 2-dimensional subgroup SO(2)xSO(2) of Aut(S^5).
Notice that S^5 is a quadric of signature (1,5+1). 


### Example 3: Computing G-invariant quadratic forms in the ideal of the Veronese surface

We compute quadratic forms in the ideal of the Veronese surface that are
invariant under SO(3) and under SL(3). We assume that the real structure is
defined by complex conjugation. 

```python
from moebius_aut.class_ma_ring import MARing
from moebius_aut.class_veronese import Veronese

so3_lst = Veronese.get_c_lst_lst_dct()['SO3(R)']
iq_lst = Veronese.get_invariant_qf( so3_lst )
print( iq_lst )
print( MARing.get_rand_sigs( iq_lst, 1 ) )

sl3_lst = Veronese.get_c_lst_lst_dct()['SL3(C)']
print( Veronese.get_invariant_qf( sl3_lst ) )
```
Output: 

    [x1^2 + x2^2 + x3^2 - x0*x4 - x0*x5 - x4*x5]
    [[1, 5]]    
    []
    
The output shows that there exists a quadric of signature (1,5) in the ideal of the Veronese surface, 
such that this quadric is invariant under the group SO(3). 
Moreover, we have shown that there exists no quadric that is invariant under the full automorphism group SL(3) 
of the Veronese surface.
      

