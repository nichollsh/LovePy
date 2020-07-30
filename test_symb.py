from sympy.interactive.printing import init_printing
init_printing(use_unicode=False, wrap_line=False, no_global=True)

import numpy as np 
from numpy import pi
from sympy import *
from sympy import Matrix




G, n, r, R, mu, g, den, pi = symbols('G n r R mu g den pi')
g = 4*pi*G*den*r/3

D = zeros(6,6)

D[0, 0] = ( n+1 ) / r**( n+1 )
D[1, 1] = n * ( n+1 ) * r**(-n + 1) / (2 * (2*n - 1) )
D[2, 2] = r**( -(n-1) )
D[3, 3] = n * r**n
D[4, 4] = n * ( n+1 ) * r**(n+2) / ( 2 * (2*n + 3) )
D[5, 5] = -r**(n+1)

D /= 2*n + 1

Ybar = zeros(6,6)

# First column:

Ybar[0, 0] = den*g*r/mu - 2*( n + 2 )
Ybar[1, 0] = -den*g*r/mu + 2*( n**2 + 3*n - 1) / ( n+1 )
Ybar[2, 0] = 4*pi*G*den
Ybar[3, 0] = den*g*r/mu + 2*( n - 1 )
Ybar[4, 0] = -den*g*r/mu - 2*(n**2 - n - 3) / n
Ybar[5, 0] = 4*pi*G*den*r

# Second column

Ybar[0, 1] = 2*n*( n+2 )
Ybar[1, 1] = -2*( n**2 - 1 )
Ybar[3, 1] =  2*( n**2 - 1 )
Ybar[4, 1] = -2*n*( n + 2 )

# Third column

Ybar[0, 2] = -r/mu
Ybar[1, 2] =  r/mu
Ybar[3, 2] = -r/mu
Ybar[4, 2] =  r/mu

# Fourth column

Ybar[0, 3] = n*r / mu
Ybar[1, 3] = ( 2-n ) * r / mu
Ybar[3, 3] =  -( n+1 ) * r / mu
Ybar[4, 3] = ( n+3 ) * r /mu

# Fifth column

Ybar[0, 4] =  den * r / mu
Ybar[1, 4] = -den * r / mu
Ybar[3, 4] =  den * r / mu
Ybar[4, 4] = -den * r / mu
Ybar[5, 4] = 2*n + 1 

# Sixth column 

Ybar[2, 5] = -1.0
Ybar[5, 5] = -r




Y = zeros(6,6)

# First column:

Y[0, 0] = n*r**( n+1 ) / ( 2*( 2*n + 3) )
Y[1, 0] = ( n+3 )*r**( n+1 ) / ( 2*( 2*n+3 ) * ( n+1 ) )
Y[2, 0] = ( n*den*g*r + 2*( n**2 - n - 3)*mu ) * r**n / ( 2*( 2*n + 3) )
Y[3, 0] = n *( n+2 ) * mu * r**n / ( ( 2*n + 3 )*( n+1 ) )
Y[5, 0] = 2*pi*G*den*n*r**( n+1 ) / ( 2*n + 3 )

# Second column

Y[0, 1] = r**( n-1 )
Y[1, 1] = r**( n-1 ) / n
Y[2, 1] = ( den*g*r + 2*( n-1 )*mu ) * r**( n-2 )
Y[3, 1] = 2*( n-1 ) * mu * r**( n-2 ) / n
Y[5, 1] = 4*pi*G*den*r**( n-1 )

# Third column

Y[2, 2] = -den * r**n
Y[4, 2] = -r**n
Y[5, 2] = -( 2*n + 1) * r**( n-1 )

# Fourth column

Y[0, 3] = ( n+1 ) * r**( -n ) / ( 2*( 2*n - 1 ) )
Y[1, 3] = ( 2 - n ) * r**( -n ) / ( 2*n * ( 2*n - 1 ) )
Y[2, 3] = ( ( n+1 ) * den*g*r - 2 * ( n**2 + 3*n - 1 ) * mu ) / ( 2*( 2*n -1 ) * r**(n+1) )
Y[3, 3] = ( n**2 - 1 ) * mu / ( n*( 2*n - 1 ) * r**( n+1) )
Y[5, 3] = 2*pi*G*den*( n+1 ) / ( ( 2*n - 1 ) * r**n )

# Fifth column

Y[0, 4] =  r**( -n - 2 )
Y[1, 4] = -r**( -n - 2 ) / ( n+1 )
Y[2, 4] = (den*g*r - 2*( n+2 )*mu ) / r**(n+3)
Y[3, 4] = 2 * ( n+2 ) * mu / ( (n+1) * r**(n+3) )
Y[5, 4] = 4*pi*G*den / r**( n+2 )

# Sixth column 

Y[2, 5] = - den / r**( n+1 )
Y[4, 5] = - 1.0 / r**( n+1 )


# Yinv = Y.inv()

Yinv2 = (D*Ybar).applyfunc(simplify)

# I = Y*Yinv

print("Evaluating")
# G, n, r, R, mu, g, den, pi = symbols('G n r R mu g den pi')
print(Yinv2.evalf(subs={G:6.67e-11, n:2, r:1000e3, mu:1e9, den:1e3, pi:np.pi}))

# print(Y*Yinv)



