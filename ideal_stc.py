'''
This file contains an implementation of the string-to-constant-weight-vector conversion function described in
Alessandro Barenghi and Gerardo Pelosi. “Constant weight strings in constant time: a building block for code-based post-quantum cryptosystems”. In: Proceedings of the
17th ACM International Conference on Computing Frontiers (2020).
Note that there are slight differences from the algorithm described in the figures in the paper; this follows more closely the prose description in the paper than the pseudocode.
'''

from sage.all_cmdline import *
from math import ceil, log2
from random import randrange
import auxiliary

#Takes as input a boolean cond, and two values t and f 
#Returns t if cond is true, f otherwise
def ctcond(cond, t, f):
	if cond:
		return t
	else:
		return f

#Stores an element in an array at the given position        
def ctstore(arr, pos, v):
    #print("ctstore pos", pos)
    arr[pos] = v

#Returns the element stored at a given position in the array    
def ctload(arr, pos):
    #print("ctload pos", pos)
    return arr[pos]

#Returns a random value in the given range    
def ctrand(v):
    #print(v)
    if v == 0:
        return 0
    return randrange(int(v) + 1)

#Chooses appropriate values for l, d based on n, t    
def fix_l_d(n, t):
    u = int((log2(n - t) - 1) / 2)
    d = int(2 ** u)
    lim1 = t * (1 + int(log2(d)))
    lim2 = int((n - t) / d)
    low = min(lim1, lim2)
    l = int((low - 1) / 8) * 8 
    assert (l * d) <= (n - t)
    return l, d

#Converts a binary string B to a vector of weight t (represented as run-length encodings)   
def StC(B, d, n, t):
    l = len(B)
    lambdaVec = [0] * t
    
    qdone = 0
    rdone = 0
    lambdadone = 0
    
    q = 0
    r = 0
    lam = 0
    
    idx = 0
    rbitctr = 0
    remainingpos = n - t
    
    for ind in range(len(B)):
        ch = B[ind]
        b = int(ch)
        qdone = qdone | (1 - b)
        q = q + (b & (1 - qdone))
        rbitctr = rbitctr + qdone 
        rdone = (rbitctr == (int(log2(d)) + 1))
        r = 2 * r + (b & qdone)
        lam = q * d + r 
        ctstore(lambdaVec, idx, lam)
        lambdadone = qdone & rdone 
        idx = idx + ctcond(lambdadone, 1, 0)
        remainingpos = remainingpos - ctcond(lambdadone, lam, 0)
        q = ctcond(lambdadone, 0, q)
        qdone = ctcond(lambdadone, 0, qdone)
        r = ctcond(lambdadone, 0, r)
        rbitctr = ctcond(lambdadone, 0, rbitctr)
    
    #Now bar(lambda)
    casepq = 1 - qdone #Part of the quotient is in s
    rpq = ctrand(remainingpos - lam)
    lam = lam + ctcond(casepq, rpq, 0)
    casecq = qdone & (rbitctr == 1) #Whole quotient is in s
    rcq = min(d - 1, remainingpos - q * d)
    rcq = ctrand(rcq)
    lam = lam + ctcond(casecq, rcq, 0)
    casepr = qdone & (1 - rdone) & (rbitctr > 1) #Whole quotient + part of remainder is in s
    pow2r = int(2 ** (log2(d) - rbitctr + 1))
    rpr = q * d + r * pow2r + ctrand(pow2r - 1)
    lam = ctcond(casepr, rpr, lam)
    lambdaVec[idx] =  lam
    idx = idx + ctcond(lambdadone, 0, 1)
    remainingpos = remainingpos - lam
    
    for i in range(t):
        r = ctrand(remainingpos)
        lam = ctload(lambdaVec, i)
        lam = ctcond((i > idx), r, lam)
        ctstore(lambdaVec, i, lam)
        remainingpos = remainingpos - ctcond((i > idx), lam, 0)
    
    return lambdaVec

#Converts a vector of weight t (represented as run-length encodings) to a binary string of fixed length   
def CtS(lambdaVec, d, n, t, l):
    B = ''
    for lam in lambdaVec:
        sub = ''
        q = int(lam / d)
        r = bin(lam % d)[2:]
        while len(r) < int(log2(d)):
            r = '0' + r
        unary_q = ''
        for j in range(q):
            unary_q = unary_q + '1'
        sub = sub + unary_q + '0' + r 
        B = B + sub 
    relevant_B = B[:l]
    return relevant_B

#Test of the correctness and invertibility of the conversion function    
def test(n, t):
    num_iter = 10000
    for i in range(num_iter):
        if (i % (num_iter / 10)) == 0:
            print(i, "iterations...")
        l, d = fix_l_d(n, t)
        B = ''
        for j in range(l):
            B = B + str(randrange(2))
        lv = StC(B, d, n, t)
        assert sum(lv) <= (n - t)
        z = auxiliary.positional_to_vector(lv, n)
        assert vector(z).hamming_weight() == t
        b = CtS(lv, d, t, n, l)
        assert b == B
                    
#test(1024, 38)
#test(2048, 69)
#test(4096, 128)
#lv = StC('100110111000011111001110111110', 32, 4096, 128)        #
#print(lv)
#print(CtS(lv, 32, 128, 4096, 30))