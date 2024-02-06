# import libraries
import numpy as np
    
# define routines
def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance'''

    count = 0
    while (count <Nmax):
       count = count +1
       x1 = f(x0)
       if (abs(x1-x0) <tol):
          xstar = x1
          ier = 0
          return [xstar,ier]
       x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier]
    
"""

# use routines 
f1 = lambda x: 1+0.5*np.sin(x)
''' 
fixed point is alpha1 = 1.4987....
'''

f2 = lambda x: 3+2*np.sin(x)
''' 
fixed point is alpha2 = 3.09... 
'''

Nmax = 100
tol = 1e-6

''' test f1 '''
x0 = 0.0
[xstar,ier] = fixedpt(f1,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f1(xstar):',f1(xstar))
print('Error message reads:',ier)
    
''' test f2 '''
x0 = 0.0
[xstar,ier] = fixedpt(f2,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f2(xstar):',f2(xstar))
print('Error message reads:',ier)

"""


f_a = lambda x: x * (1 + ((7 - x**5) / x**2))**3
f_b = lambda x: x - ((x**5 - 7) / x**2)
f_c = lambda x: x - ((x**5 - 7) / (5 * x**4))
f_d = lambda x: x - ((x**5 - 7) / 12)

fixed_point = 7**(1/5)

tol = 1e-10
Nmax = 100

x0 = 1
[xstar,ier] = fixedpt(f_d,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f(xstar):',f_d(xstar))
print('Error message reads:',ier)