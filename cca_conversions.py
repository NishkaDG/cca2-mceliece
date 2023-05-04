'''
This file contains various IND-CCA2 secure transformations of the Classic McEliece protocol, using a variety of conversion functions. 
References:
 - D. Engelbert, R. Overbeck, and A. Schmidt. A Summary of McEliece-Type Cryptosystems and their Security. Cryptology ePrint Archive, Paper 2006/162. https://e
print.iacr.org/2006/162. 2006. url: https://eprint.iacr.org/2006/162.
 - Kazukuni Kobara and Hideki Imai. “Semantically Secure McEliece Public-Key Cryptosystems - Conversions for McEliece PKC”. In: International Conference on Theory
and Practice of Public Key Cryptography. 2001.
 - Nicolas Sendrier. “Encoding information into constant weight words”. In: Proceedings. International Symposium on Information Theory, 2005. ISIT 2005. 2005,
pp. 435–438. doi: 10.1109/ISIT.2005.1523371.
 - Alessandro Barenghi and Gerardo Pelosi. “Constant weight strings in constant time: a building block for code-based post-quantum cryptosystems”. In: Proceedings of the
17th ACM International Conference on Computing Frontiers (2020).
 - Pierre-Louis Cayrel, Gerhard Hoffmann, and Edoardo Persichetti. “Efficient Implementation of a CCA2-Secure Variant of McEliece Using Generalized Srivastava Codes”. 
 In: International Conference on Theory and Practice of Public Key Cryptography. 2012.
'''

from sage.all_cmdline import *   # import sage library
from math import ceil, floor, log2
import numpy as np
import timeit

import classic
import sendrier
import ideal_stc
import auxiliary

error_vec_list = []

#The Fujisaki-Okamoto transform with Sendrier's function for converting bitstrings to constant-weight vectors
def fujisaki_okamoto_encrypt_sendrier(m, n, k, pk):
    t = pk[1]
    #Generate r
    r = random_matrix(Integers(2), 1, k)
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    z1 = bin(auxiliary.H(in1, n, t))[2:] #BtoCW takes binary strings as input 
    z2 = sendrier.BtoCW(n, t, 0, z1, 0)
    z = auxiliary.positional_to_vector(z2, n)
    assert vector(z).hamming_weight() == t
    #z = Conv(z, n)
    #print('r', r, 'z', z)
    c1 = classic.encrypt(r, z, pk)
    in2 = auxiliary.vector_to_bytes(r)
    c2 = auxiliary.R(in2, k) + m
    return c1, c2

#The Fujisaki-Okamoto transform with Barenghi and Pelosi's function for converting bitstrings to constant-weight vectors
def fujisaki_okamoto_encrypt_ideal(m, n, k, pk):
    t = pk[1]
    l, d = ideal_stc.fix_l_d(n, t)
    #print(l, d)
    #Generate r
    r = random_matrix(Integers(2), 1, k)
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    B = bin(auxiliary.H1(in1, l))[2:] #StC takes binary strings as input 
    lv = ideal_stc.StC(B, d, n, t)
    z = auxiliary.positional_to_vector(lv, n)
    assert vector(z).hamming_weight() == t
    #z = Conv(z, n)
    #print('r', r, 'z', z)
    c1 = classic.encrypt(r, z, pk)
    in2 = auxiliary.vector_to_bytes(r)
    c2 = auxiliary.R(in2, k) + m
    return c1, c2

#The Fujisaki-Okamoto transform that does not use the conversion function (from Cayrel et al)
def alt_fujisaki_okamoto_encrypt(m, n, k, pk):
    t = pk[1]
    r = random_matrix(Integers(2), 1, n)
    classic.select_error(r, t, n)
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    z = matrix(Integers(2), [auxiliary.pad_as_list(auxiliary.H1(in1, k), k)])
    c1 = classic.encrypt(z, r, pk)
    in2 = auxiliary.vector_to_bytes(r)
    c2 = auxiliary.R(in2, k) + m
    return c1, c2

#The Kobara-Imai gamma transform, implemented with the Barenghi-Pelosi conversion
#Note that this does not work in its current form, as the Barenghi-Pelosi method requires bitstrings to be of length < log(C(n,t)) 
#and the gamma transform sends inputs of length exactly log(C(n,t)) for conversion    
def kobara_imai_gamma_encrypt(m, n, k, const, pk):
    r_len = 160
    const_len = 160
    m_len = m.ncols()
    assert const.ncols() == const_len
    t = pk[1]
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    lognct = int(floor(log(nct, 2).n()))

    r = random_matrix(Integers(2), 1, r_len)
    #l, d = ideal_stc.fix_l_d(n, t)
    l = lognct - 10
    lim = int(n - t) / l
    u = int(log2(lim))
    d = int(2 ** u)

    c1_len = m_len + const_len
    c2_len = r_len
    c3_len = k
    c4_len = lognct
    c5_len = m_len + const_len + r_len - c4_len - k
    #c6_len = c1_len + c2_len - lognct - k

    c1 = auxiliary.R(auxiliary.vector_to_bytes(r), c1_len) + auxiliary.concat_vectors(m, const)
    c2 = r + auxiliary.R(auxiliary.vector_to_bytes(c1), r_len)
    c2c1 = auxiliary.concat_vectors(c2, c1)
    assert c2c1.ncols() == (c1.ncols() + c2.ncols())
    #print(c2c1.ncols(), c3_len, lognct)
    c3 = auxiliary.LSB(c2c1, c3_len)
    c5c4 = auxiliary.MSB(c2c1, c5_len + c4_len)
    c4 = auxiliary.LSB(c5c4, c4_len)
    c5 = auxiliary.MSB(c5c4, c5_len)
    zpos = ideal_stc.StC(vector_to_bitstring(c4), d, n, t)
    z = auxiliary.positional_to_vector(zpos, n)
    if c5_len > 0:
        #c6 = MSB(c2c1, c6_len)
        e = classic.encrypt(c3, z, pk)
        c = auxiliary.concat_vectors(c5, e)
        return c
    else:
        c = classic.encrypt(c3, z, pk)
        return c

#The Kobara-Imai alpha protocol, implemented with the Barenghi-Pelosi conversion
def kobara_imai_alpha_encrypt(m, n, k, pk):
    t = pk[1]
    l, d = ideal_stc.fix_l_d(n, t)
    r_len = 160
    m_len = m.ncols()
    r = random_matrix(Integers(2), 1, r_len)
    out1 = auxiliary.H(auxiliary.concat_vectors_to_bytearray(r, m), n, t)
    zbarbin = auxiliary.pad_as_bitstring(out1, l)
    zbar = zbarbin[:l]
    zbar_bytes = auxiliary.bitstring_to_bytes(zbar)
    y1y2 = auxiliary.R(zbar_bytes, r_len + m_len) + auxiliary.concat_vectors(r, m)
    y1 = auxiliary.MSB(y1y2, k)
    y2 = auxiliary.LSB(y1y2, r_len + m_len - k)
    lv = ideal_stc.StC(zbar, d, n, t)
    z = auxiliary.positional_to_vector(lv, n)
    c1 = classic.encrypt(y1, z, pk)
    c2 = y2
    return c1, c2

#The combinadics approach to creating a conversion function            
def generate_all_error_vecs(n, t):
    limit = 1 << n
    val = (1 << t) - 1
    list_all = []
    while val < limit:
        yield int("{0:0{1}b}".format(val,n), 2)
        minbit=val&-val #rightmost 1 bit
        fillbit = (val+minbit)&~val  #rightmost 0 to the left of that bit
        val = val+minbit | (fillbit//(minbit<<1))-1

#Test of the Fujisaki-Okamoto protocol as implemented with the two different conversion functions        
def test_original_f_o(n, t, k):
    #n = 4096
    #t = 128
    #k = 2560
    print("Testing the Fujisaki-Okamoto transform with n=", n, "t=", t, "k=", k)
    m = random_matrix(Integers(2), 1, k)
    pk, sk = classic.keygen(n, t, k)
    print("Classic McEliece key generated...")
    
    num_iter = 10000
    duration_sendrier = 0
    duration_ideal = 0
    
    for i in range(num_iter):
        if i % (num_iter / 10) == 0:
            print(i, "iterations...")
        #Timing Sendrier
        start = timeit.default_timer()
        c1, c2 = fujisaki_okamoto_encrypt_sendrier(m, n, k, pk)
        stop = timeit.default_timer()
        duration_sendrier += (stop - start)
        #print(c1)
        #print(c2)
        
        #Timing Barenghi-Pelosi
        start = timeit.default_timer()
        c1, c2 = fujisaki_okamoto_encrypt_ideal(m, n, k, pk)
        stop = timeit.default_timer()
        duration_ideal += (stop - start)
        #print(c1)
        #print(c2)
    print("Average encryption time with Sendrier", duration_sendrier / num_iter)
    print("Average encryption time with Barenghi-Pelosi", duration_ideal / num_iter)

#Test of the conversion-free Fujisaki-Okamoto transform    
def test_alt_f_o():
    #n = 1024
    #t = 50
    #k = 524
    n = 4096
    t = 128
    k = 2560
    m = random_matrix(Integers(2), 1, k)
    pk, sk = classic.keygen(n, t, k)
    print("Classic McEliece key generated...")
    
    num_iter = 10
    duration = 0
    
    for i in range(num_iter):
        start = timeit.default_timer()
        c1, c2 = alt_fujisaki_okamoto_encrypt(m, n, k, pk)
        stop = timeit.default_timer()
        duration += stop - start 
    #print(c1)
    #print(c2)
    print("Without Conversion", duration / num_iter)

#Test of the Kobara-Imai alpha protocol
def test_kobara_imai(n, t, k):
    #n = 4096
    #t = 128
    #k = 2560
    
    print("Testing the Kobara-Imai alpha transform with n=", n, "t=", t, "k=", k)
    m = random_matrix(Integers(2), 1, k + 500)
    pk, sk = classic.keygen(n, t, k)
    const = random_matrix(Integers(2), 1, 160)
    
    num_iter = 10000
    duration = 0
    
    for i in range(num_iter):
        start = timeit.default_timer()
        c1, c2 = kobara_imai_alpha_encrypt(m, n, k, pk)
        stop = timeit.default_timer()
        duration += stop - start
    print("Average encryption time of Kobara-Imai alpha", duration / num_iter)

#test_original_f_o(1024, 38, 644)
#test_original_f_o(2048, 69, 1289)
#test_original_f_o(4096, 128, 2560)
#test_alt_f_o()
test_kobara_imai(1024, 28, 644)
test_kobara_imai(2048, 69, 1289)
test_kobara_imai(4096, 128, 2560)