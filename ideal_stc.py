from sage.all_cmdline import *
from math import ceil, log2
from random import randrange
import classic

def ctcond(cond, t, f):
	if cond:
		return t
	else:
		return f
        
def ctstore(arr, pos, v):
    arr[pos] = v
    
def ctload(arr, pos):
    return arr[pos]
    
def ctrand(v):
    if v == 0:
        return 0
    return randrange(int(v) + 1)
    
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
        #print(ind, b, idx)
        qdone = qdone | (1 - b)
        #print(qdone)
        q = q + (b & (1 - qdone))
        #print(q)
        rbitctr = rbitctr + qdone 
        rdone = (rbitctr == (int(log2(d)) + 1))
        r = 2 * r + (b & qdone)
        lam = q * d + r 
        ctstore(lambdaVec, idx, lam)
        #print(lambdaVec)
        lambdadone = qdone & rdone 
        idx = idx + ctcond(lambdadone, 1, 0)
        #print(qdone, rdone, lambdadone, q, r)
        remainingpos = remainingpos - ctcond(lambdadone, lam, 0)
        q = ctcond(lambdadone, 0, q)
        qdone = ctcond(lambdadone, 0, qdone)
        r = ctcond(lambdadone, 0, r)
        rbitctr = ctcond(lambdadone, 0, rbitctr)
        #print(qdone, rdone, lambdadone, q, r)
    
    #Now bar(lambda)
    casepq = 1 - qdone #Part of the quotient is in s
    rpq = ctrand(remainingpos - lam)
    lam = lam + ctcond(casepq, rpq, 0)
    #lambdaVec[idx] =  lam
    #print(lam, casepq, rpq)
    casecq = qdone & (rbitctr == 1) #Whole quotient is in s
    #rcq = ((rpq > (d - 1)) & (d - 1)) | ((1 - (rpq > (d - 1))) & rpq)
    rcq = min(d - 1, remainingpos - q * d)
    rcq = ctrand(rcq)
    lam = lam + ctcond(casecq, rcq, 0)
    #lambdaVec[idx] =  lam
    #print(lam, casecq, rcq)
    casepr = qdone & (1 - rdone) & (rbitctr > 1) #Whole quotient + part of remainder is in s
    pow2r = int(2 ** (log2(d) - rbitctr + 1))
    #print(q, d, r, pow2r)
    rpr = q * d + r * pow2r + ctrand(pow2r - 1)
    lam = ctcond(casepr, rpr, lam)
    lambdaVec[idx] =  lam
    #print("Case:", casepq, casecq, casepr, lam, rpq, rcq, rpr)
    #print(lambdaVec)
    idx = idx + ctcond(lambdadone, 0, 1)
    remainingpos = remainingpos - ctcond(lambdadone, lam, 0)
    
    for i in range(t):
        r = ctrand(remainingpos)
        lam = ctload(lambdaVec, i)
        lam = ctcond((i > idx), r, lam)
        ctstore(lambdaVec, i, lam)
        remainingpos = remainingpos - ctcond((i > idx), lam, 0)
    
    return lambdaVec
    
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
    
def test():
    n = 4096
    t = 128
    for i in range(100):
        u = randrange(int(log2(n - t) - 1)) + 1
        d = int(2 ** u)
        #print(d)
        lim1 = t * (1 + int(log2(d)))
        lim2 = int((n - t) / d)
        low = min(lim1, lim2)
        l = randrange(low)
        assert (l * d) <= (n - t)
        B = ''
        for j in range(l):
            B = B + str(randrange(2))
        lv = StC(B, d, n, t)
        b = CtS(lv, d, t, n, l)
        if not (b == B):
            print("Error:", i, n, t, d, l, B, lv, b)
            
test()
#lv = StC('100110111000011111001110111110', 32, 4096, 128)        #
#print(lv)
#print(CtS(lv, 32, 128, 4096, 30))