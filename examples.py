from __future__ import division
from polynomial import *

p1=Polynomial([0,2,4]) # polynomial of degree 2: 4*x^1 + 2*x^2 + 0
print repr(p1) # interpreter representation

# auto trimming of highest degrees with zero coeficients
p1_trim=Polynomial([0,2,4,0]) # polynomial of degree 2: 4*x^1 + 2*x^2 + 0
print 'p1 and p1_trim have the same representation: '
print 'represent in console --> p1: ', p1, 'p1_trim: ', p1_trim

# evaluate polynomial at value:
print 'evalute at 2: ', p1(2)
# polynomial stores dictionary of all terms indexed by term order (exponent)
print 'terms: ', p1.terms

p2 = Polynomial([3,1,2,8])

print '--- Addition ---'
print p1, ' + ', p2, ' : '
p3 = p1+p2
print 'result: ', p3

print '--- Substraction ---'
print p1, ' - ', p2, ' : '
p4= p1-p2
print 'result: ', p4

print '--- Polynomial Multiplication ---'
print p1, ' * ', p2, ' : '
p5 = p1*p2
print 'result: ', p5
print '--- Polynomial Division ---'
print p2, ' / ', p1, ' : '
p6 = p2/p1
print 'result: ', p6

print '--- Add quotient and reminder back ---'
print 'quotient:: ', p6[0]
print 'reminder:: ', p6[1]
print 'original:: ', p2
print 'reconstruct:: ', p6[0]*p1 + p6[1]

