'''
An implementation of the string-to-constant-weight-vector conversion protocol found in 
Nicolas Sendrier. “Encoding information into constant weight words”. In: Proceed-
ings. International Symposium on Information Theory, 2005. ISIT 2005. 2005,
pp. 435–438. doi: 10.1109/ISIT.2005.1523371
'''

from sage.all_cmdline import *   # import sage library
from math import ceil, log2
from random import randrange
import classic

#Takes as input an integer x and a bitlength u 
#Returns the u least significant bits of the binary conversion of the integer
def base2(x, u):
    binx = bin(x)[2:]
    #print(binx, len(binx), u)
    if u > 0:
        if len(binx) >= u:
            return binx[-u:]
        else:
            while len(binx) < u:
                binx = '0' + binx
            return binx
    else:
        return None

#Takes as input a bitstring B, a number of bits to read u, and a starting index 
#Returns a substring of u bits read from index start, converted to a decimal integer
def read_bits(B, u, start):
    #print('read_bits', B, u, start)
    if start < 0:
        return 0
    substr = B[start:start+u]
    if len(substr) > 0:
        return int(substr, 2)
    else:
        return 0
 
#Takes as input n, t 
#Returns an optimal value of d  
def best_d(n, t):
    d = ceil((n - ((t - 1) / 2)) * (1 - (1 / (2 ** (1 / t)))))
    #ex = int(ceil(log2(d))) - 1
    #d = 2 ** ex
    assert (1 <= d) and (d <= (n - t))
    return d

#An implementation of f_d() function in the paper    
def encode_fd(delta, d):
    u = ceil(log2(d))
    #print('u', u)
    limit = (2 ** u) - d
    #print('delta', delta, 'limit', limit)
    if delta < limit:
        u = u - 1
    else:
        delta = delta + limit
    #print('delta', delta)
    #print('base2', base2(delta, u))
    res = base2(delta, u)
    if res is None:
        return ''
    else:
        return res

#An implementation of the inverse of the f_d function in the paper    
def decode_fd(d, B, start):
    u = ceil(log2(d))
    #print('d', d, 'u', u)
    delta = read_bits(B, u-1, start)
    start = start + u - 1
    limit = (2 ** u) - d
    if delta >= limit:
        delta = (2 * delta) + read_bits(B, 1, start) - (2 ** u) + d
        start = start + 1
    return delta, start

#A recursive function to convert a constant-weight vector (expressed as run-length encodings) to a binary string    
def CWtoB(n, t, delta_tuple):
    if (t == 0) or (n <= t):
        return ''
    d = best_d(n, t)
    delta_1 = delta_tuple[0]
    #print('n', n, 't', t, 'd', d, 'delta_1', delta_1)
    if delta_1 >= d:
        new_delta_lst = list(delta_tuple)
        new_delta_lst[0] = delta_1 - d
        new_delta_tuple = tuple(new_delta_lst)
        #print('new_delta_tuple', new_delta_tuple)
        res = '1' + CWtoB(n - d, t, new_delta_tuple)
        #print('n', n, 't', t, 'd', d, 'res', res)
        return res
    else:
        enc = encode_fd(delta_1, d)
        s = '0'
        if not (enc is None):
            s = '0' + enc
        #print('s', s)
        new_delta_lst = list(delta_tuple)
        new_delta_lst = new_delta_lst[1:]
        new_delta_tuple = tuple(new_delta_lst)
        #print('new_delta_tuple', new_delta_tuple)
        res = s + CWtoB(n - delta_1 - 1, t - 1, new_delta_tuple)
        #print('n', n, 't', t, 'd', d, 'res', res, 'delta_1', delta_1, n-delta_1 - 1)
        return res

#A recursive function to convert an arbitrary binary string to a list of run-length encodings representing a vector of weight t        
def BtoCW(n, t, delta, B, start):
    #print('n', n, 't', t, 'delta', delta, 'B', B)
    if t == 0:
        return []
    elif n <= t:
        res = [delta] + BtoCW(n - 1, t - 1, 0, B, start)
        #print('n', n, 't', t, 'delta', delta, 'B', B, 'res', res)
        return res
    else:
        d = best_d(n, t)
        next_bit = read_bits(B, 1, start)
        #print('next_bit', next_bit, 'd', d)
        start = start + 1
        if next_bit == 1:
            res = BtoCW(n - d, t, delta + d, B, start)
            #print('n', n, 't', t, 'd', d, 'delta', delta, 'B', B, 'res', res)
            return res
        else:
            i, start = decode_fd(d, B, start)
            #print('i', i, 'd', d)
            res = [delta + i] + BtoCW(n - i - 1, t - 1, 0, B, start)
            #print('n', n, 't', t, 'd', d, 'i', i, 'start', start, 'delta', delta, 'B', B, 'res', res)
            return res

#Test that the functions for f_d and its inverse work correctly for random inputs
def test_decode_encode_fd():
    for i in range(100):
        d = randrange(2048)
        delta = randrange(d)
        B = encode_fd(delta, d)
        res, index = decode_fd(d, B, 0)
        assert res == delta

#Test that the functions f_d and its inverse work with bijection in the other direction        
def test_encode_decode_fd():
    for i in range(100):
        exp = randrange(12)
        d = 2 ** exp
        u = ceil(log2(d))
        r = random_matrix(Integers(2), 1, u)
        B = ''
        for ele in r[0]:
            B = B + str(ele)
        delta, index = decode_fd(d, B, 0)
        res = encode_fd(delta, d)
        if not (res == B):
            print(d, u, delta, res, B)

#Test that the string-to-constant-weight-vector function always returns a vector of hamming weight t
def test_BtoCW_weight():
    n = 2048
    t = 29
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    bitlength = int(ceil(log2(nct)))
    #print(bitlength)
    for i in range(100):
        r = random_matrix(Integers(2), 1, bitlength)
        B = ''
        for ele in r[0]:
            B = B + str(ele)
        #print(len(B))
        B = int(B, 2) % nct
        B = bin(B)[2:]
        while len(B) < bitlength:
            B = '0' + B
        #print(len(B))
        delta_lst = BtoCW(n, t, 0, B, 0)
        vec = matrix(1, n)
        delta_agg = []
        for i in range(len(delta_lst)):
            s = delta_lst[i]
            for j in range(i):
                s = s + delta_lst[j] + 1
            delta_agg.append(s)
        for ele in delta_agg:
            vec[0, ele] = 1
        #print(delta_lst, delta_agg, vec)
        assert vector(vec).hamming_weight() == t

#Test that the function for string-to-constant-weight-vector conversion always returns the same output on the same input string 
def test_BtoCW_unique():
    n = 2048
    t = 29
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    bitlength = int(ceil(log2(nct)))
    r = random_matrix(Integers(2), 1, bitlength)
    B = ''
    for ele in r[0]:
        B = B + str(ele)
    #print(len(B))
    B = int(B, 2) % nct
    B = bin(B)[2:]
    while len(B) < bitlength:
        B = '0' + B
    #print(bitlength)
    delta_lst = BtoCW(n, t, 0, B, 0)
    for i in range(99):
        #print(len(B))
        delta_lst_2 = BtoCW(n, t, 0, B, 0)
        assert delta_lst == delta_lst_2

#Test that the function for constant-weight-vector-to-string conversion always returns the same output on the same input vector
def test_CWtoB_unique():
    n = 2048
    t = 29
    r = random_matrix(Integers(2), 1, n)
    classic.select_error(r, t, n)
    delta_lst = []
    ctr = 0
    for ele in r[0]:
        if ele == 1:
            delta_lst.append(ctr)
            ctr = 0
        else:
            ctr = ctr + 1
    delta_tuple = tuple(delta_lst)
    B = CWtoB(n, t, delta_tuple)
    for i in range(99):
        B_1 = CWtoB(n, t, delta_tuple)
        assert B == B_1

#Test that BtoCW(CWtoB()) returns the initial error vector  
def test_conversion_bijective():
    n = 2048
    t = 29
    #n = 30
    #t = 5
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    bitlength = int(ceil(log2(nct)))
    print(bitlength)
    for i in range(10):
        print("Test", i)
        #r = random_matrix(Integers(2), 1, bitlength)
        r = random_matrix(Integers(2), 1, n)
        classic.select_error(r, t, n)
        delta_lst = []
        ctr = 0
        for ele in r[0]:
            if ele == 1:
                delta_lst.append(ctr)
                ctr = 0
            else:
                ctr = ctr + 1
        delta_tuple = tuple(delta_lst)
        B = CWtoB(n, t, delta_tuple)
        res = BtoCW(n, t, 0, B, 0)
        assert (res == delta_lst)

#Test that CWtoB(BtoCW()) returns the initial binary string 
#Note that this does not work as Sendrier's method has hard-to-predict failures on binary strings
def test_reverse_bijection():
    n = 30
    t = 5
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    bitlength = int(ceil(log2(nct))) - 1
    for i in range(1):
        r = random_matrix(Integers(2), 1, bitlength)
        B = ''
        B = '00010011111100001'
        delta_lst = BtoCW(n, t, 0, B, 0)
        res = CWtoB(n, t, tuple(delta_lst))
        if not (res == B):
            print(delta_lst)
            print(res, B)
        
'''
test_decode_encode_fd()
print("decode_fd(encode_fd()) works")
test_encode_decode_fd()
print("encode_fd(decode_fd()) works")


test_BtoCW_weight()
print("BtoCW returns values of correct weight")
test_BtoCW_unique()
print("BtoCW always returns the same output on an input")
test_CWtoB_unique()
print("CWtoB always returns the same output on an input")

test_conversion_bijective()
print("BtoCW(CWtoB()) works")

test_reverse_bijection()
'''