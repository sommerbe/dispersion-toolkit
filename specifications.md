## File format of a point set sequence

The bounded domain is spanned by a tuple of smallest possible coordinates,
````
(low_0 low_1 ... low_n),
````
and a tuple of greatest possible coordinates
````
(up_0 up_1 ... up_n).
````

An n-dimensional point set P_i with cardinality k, having points p_j with coordinates c_{i,j,s} is serialised with ASCII encoding according to the format
````
#d low_0 low_1 ... low_n up_0 up_1 ... up_n
c_{i,0,0} c_{i,0,1} ... c_{i,0,n}
c_{i,1,0} c_{i,1,1} ... c_{i,1,n}
                     .
                     .
                     .
c_{i,k,0} c_{i,k,1} ... c_{i,k,n}
#eos
````
where columns are separated by an ASCII character named *delimiter*, here
````
' '.
````

Each line of the serialised result may additionally contain a commenting line,
````
# this is a comment.
````

A point set sequence of length m is the concatenation of serialised point sets, i.e.
````
#d low_0 low_1 ... low_n0 up_0 up_1 ... up_n0
c_{0,0,0}  c_{0,0,1}  ... c_{0,0,n0}
                       .
c_{0,k0,0} c_{0,k0,1} ... c_{0,k0,n0}
#eos
#d low_0 low_1 ... low_n1 up_0 up_1 ... up_n1
c_{1,0,0}  c_{1,0,1}  ... c_{1,0,n1}
                       .
c_{1,k1,0} c_{1,k1,1} ... c_{1,k1,n1}
#eos
                       .
                       .
#d low_0 low_1 ... low_nm up_0 up_1 ... up_nm
c_{m,0,0}  c_{m,0,1}  ... c_{m,0,nm}
                       .
c_{m,km,0} c_{m,km,1} ... c_{m,km,nm}
#eos
````