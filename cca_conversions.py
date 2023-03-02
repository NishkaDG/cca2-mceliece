from sage.all_cmdline import *   # import sage library
from Crypto.Hash import SHAKE256
from math import ceil
import numpy as np
import timeit

import classic
import auxiliary

error_vec_list = []

def pad(num, n):
    res = []
    bin_val = bin(num)[2:]
    if len(bin_val) < n:
        for i in range(n - len(bin_val)):
            res.append(0)
    else:
        bin_val = bin_val[:n]
    for i in bin_val:
        res.append(int(i))
    #print(res)
    return res

def concat(vec1, vec2):
    cc = bytearray(b'\x00')
    for i in vec1[0]:
        cc.append(i)
    for i in vec2[0]:
        cc.append(i)
    return cc

def H(bit_string, n, t):
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    bitlength = int(ceil(log(nct, 2).n()))
    bytelength = int(ceil(bitlength / 8))
    
    shake = SHAKE256.new()
    shake.update(bit_string)
    h = shake.read(bytelength).hex()
    h = int(h, 16) % nct
    #print(nct)
    return h

def H1(bit_string, k):
    bytelength = int(ceil(k / 8))
    #print(k, 2**k)
    
    shake = SHAKE256.new()
    shake.update(bit_string)
    h = shake.read(bytelength).hex()
    h = int(h, 16) % (2 ** k)
    return h

def R(r, k):
    bytelength = int(ceil(k / 8))
    shake = SHAKE256.new()
    shake.update(r)
    h_hex = shake.read(bytelength).hex()
    h_list = pad(int(h_hex, 16), k)
    h_vec = matrix(Integers(2), h_list)
    return h_vec

def fujisaki_okamoto_encrypt(m, n, k, pk):
    t = pk[1]
    #Generate r
    r = random_matrix(Integers(2), 1, k)
    in1 = concat(r, m)
    z1 = bin(H(in1, n, t))[2:] #BtoCW takes binary strings as input 
    z2 = auxiliary.BtoCW(n, t, 0, z1, 0)
    #z = auxiliary.positional_to_vector(z2, n)
    #z = Conv(z, n)
    #print('r', r, 'z', z)
    c1 = classic.encrypt(r, z, pk)
    in2 = concat(r, matrix(Integers(2), [0]))
    c2 = R(in2, k) + m
    return c1, c2
	
def alt_fujisaki_okamoto_encrypt(m, n, k, pk):
    t = pk[1]
    r = random_matrix(Integers(2), 1, n)
    classic.select_error(r, t, n)
    in1 = concat(r, m)
    z = matrix(Integers(2), [pad(H1(in1, k), k)])
    c1 = classic.encrypt(z, r, pk)
    in2 = concat(r, matrix(Integers(2), [0]))
    c2 = R(in2, k) + m
    return c1, c2
    
def kobara_imai_encrypt(m, const, pk):
	return True
    
def generate_all_error_vecs(n, t):
    limit = 1 << n
    val = (1 << t) - 1
    list_all = []
    while val < limit:
        yield int("{0:0{1}b}".format(val,n), 2)
        minbit=val&-val #rightmost 1 bit
        fillbit = (val+minbit)&~val  #rightmost 0 to the left of that bit
        val = val+minbit | (fillbit//(minbit<<1))-1
    
    
def test_original_f_o():
    global error_vec_list
    n = 30
    t = 20
    k = 15
    #all_error_gen = generate_all_error_vecs(n, t)
    #error_vec_list = np.fromiter(all_error_gen, int)
    #print("Error vectors generated...")

    #print(error_vec_list)

    #print(H(b'hellp', 50, 30))
        
    #print(R(b'123', 4))

    m = random_matrix(Integers(2), 1, k)
    pk, sk = classic.keygen(n, t, k)
    print("Classic McEliece key generated...")

    start = timeit.default_timer()
    c1, c2 = fujisaki_okamoto_encrypt(m, n, k, pk)
    stop = timeit.default_timer()
    #print(c1)
    #print(c2)
    print(stop - start)
    
def test_alt_f_o():
    #n = 1024
    #t = 50
    #k = 524
    n = 2048
    t = 29
    k = 2000
    m = random_matrix(Integers(2), 1, k)
    pk, sk = classic.keygen(n, t, k)
    print("Classic McEliece key generated...")

    start = timeit.default_timer()
    c1, c2 = alt_fujisaki_okamoto_encrypt(m, n, k, pk)
    stop = timeit.default_timer()
    #print(c1)
    #print(c2)
    print(stop - start)

#test_original_f_o()
test_alt_f_o()