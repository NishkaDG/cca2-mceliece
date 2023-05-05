'''
This file contains implementations of various IND-CCA2 secure transformations of the Classic McEliece protocol, using a variety of conversion functions. 
References:
 - D. Engelbert, R. Overbeck, and A. Schmidt. A Summary of McEliece-Type Cryptosystems and their Security. Cryptology ePrint Archive, Paper 2006/162. https://eprint.iacr.org/2006/162. 2006. url: https://eprint.iacr.org/2006/162.
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

#Encryption with the Fujisaki-Okamoto transform using Sendrier's function for converting bitstrings to constant-weight vectors
def fujisaki_okamoto_encrypt_sendrier(m, n, k, pk):
    t = pk[1]
    #Generate r
    r = random_matrix(GF(2), 1, k)
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    z1 = bin(auxiliary.H(in1, n, t))[2:] #BtoCW takes binary strings as input 
    z2 = sendrier.BtoCW(n, t, 0, z1, 0)
    z = auxiliary.positional_to_vector(z2, n)
    assert vector(z).hamming_weight() == t
    c1 = classic.encrypt(r, z, pk)
    in2 = auxiliary.vector_to_bytes(r)
    c2 = auxiliary.R(in2, k) + m
    return c1, c2
    
#Decryption with the Fujisaki-Okamoto transform using Sendrier's function for converting bitstrings to constant-weight vectors
def fujisaki_okamoto_decrypt_sendrier(c1, c2, pk, sk):
    k = pk[0].nrows()
    n = pk[0].ncols()
    t = pk[1]
    r, z = classic.decrypt(c1, sk, pk)
    in2 = auxiliary.vector_to_bytes(r)
    m = c2 + auxiliary.R(in2, k)
    
    #Now test 
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    z1 = bin(auxiliary.H(in1, n, t))[2:]
    z2 = sendrier.BtoCW(n, t, 0, z1, 0)
    expected_z = auxiliary.positional_to_vector(z2, n)
    expected_c1 = classic.encrypt(r, expected_z, pk)
    if c1 == expected_c1:
        return m
    else:
        return None

#Encryption with the Fujisaki-Okamoto transform using Barenghi and Pelosi's function for converting bitstrings to constant-weight vectors
def fujisaki_okamoto_encrypt_ideal(m, n, k, pk):
    t = pk[1]
    l, d = ideal_stc.fix_l_d(n, t)
    #Generate r
    r = random_matrix(GF(2), 1, k)
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    aux = auxiliary.H1(in1, l)
    B = auxiliary.H1(in1, l) #StC takes binary strings as input 
    lv = ideal_stc.StC(B, d, n, t)
    z = auxiliary.positional_to_vector(lv, n)
    assert vector(z).hamming_weight() == t
    c1 = classic.encrypt(r, z, pk)
    in2 = auxiliary.vector_to_bytes(r)
    c2 = auxiliary.R(in2, k) + m
    return c1, c2
    
#Decryption with the Fujisaki-Okamoto transform using Barenghi and Pelosi's function for converting bitstrings to constant-weight vectors
#Since Conv() in the forward direction in this protocol is one-to-many, we need to unconvert to check
def fujisaki_okamoto_decrypt_ideal(c1, c2, pk, sk):
    k = pk[0].nrows()
    n = pk[0].ncols()
    t = pk[1]
    l, d = ideal_stc.fix_l_d(n, t)
    
    r, z = classic.decrypt(c1, sk, pk)
    in2 = auxiliary.vector_to_bytes(r)
    m = c2 + auxiliary.R(in2, k)
        
    #Now test 
    lv = auxiliary.vector_to_positional(z)
    expected_B = ideal_stc.CtS(lv, d, t, n, l)[:l]
    
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    B = auxiliary.H1(in1, l)
    if B == expected_B:
        return m
    else:
        return None

#Encryption with the Fujisaki-Okamoto transform that does not use the conversion function (from Cayrel et al)
def alt_fujisaki_okamoto_encrypt(m, n, k, pk):
    t = pk[1]
    r = random_matrix(GF(2), 1, n)
    classic.select_error(r, t, n)
    assert vector(r).hamming_weight() == t
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    z = auxiliary.bitstring_to_vector(auxiliary.H1(in1, k))
    c1 = classic.encrypt(z, r, pk)
    in2 = auxiliary.vector_to_bytes(r)
    c2 = auxiliary.R(in2, k) + m
    return c1, c2
    
#Encryption with the Fujisaki-Okamoto transform that does not use the conversion function (from Cayrel et al)
def alt_fujisaki_okamoto_decrypt(c1, c2, pk, sk):
    k = pk[0].nrows()
    n = pk[0].ncols()
    t = pk[1]
    z, r = classic.decrypt(c1, sk, pk)
    in2 = auxiliary.vector_to_bytes(r)
    m = c2 + auxiliary.R(in2, k)
    
    #Now test 
    in1 = auxiliary.concat_vectors_to_bytearray(r, m)
    expected_z = auxiliary.bitstring_to_vector(auxiliary.H1(in1, k))
    expected_c1 = classic.encrypt(expected_z, r, pk)
    if c1 == expected_c1 and z == expected_z:
        return m
    else:
        return None 
    
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

    r = random_matrix(GF(2), 1, r_len)
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

#Encryption with the Kobara-Imai alpha protocol, implemented with the Barenghi-Pelosi conversion
def kobara_imai_alpha_encrypt(m, n, k, pk):
    t = pk[1]
    l, d = ideal_stc.fix_l_d(n, t)
    r_len = 160
    m_len = m.ncols()
    r = random_matrix(GF(2), 1, r_len)
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

#Decryption with the Kobara-Imai alpha protocol, implemented with the Barenghi-Pelosi conversion
def kobara_imai_alpha_decrypt(c1, c2, pk, sk):
    k = pk[0].nrows()
    n = pk[0].ncols()
    t = pk[1]
    l, d = ideal_stc.fix_l_d(n, t)
    
    y3, z = classic.decrypt(c1, sk, pk)
    y2 = c2
    c_len = y3.ncols() + y2.ncols()
    lv = auxiliary.vector_to_positional(z)
    zbar = ideal_stc.CtS(lv, d, t, n, l)
    rm = auxiliary.R(auxiliary.bitstring_to_bytes(zbar), c_len) + auxiliary.concat_vectors(y3, y2)
    out1 = auxiliary.H(auxiliary.vector_to_bytes(rm), n, t)
    expected_zbar = auxiliary.pad_as_bitstring(out1, l)[:l]
    if zbar == expected_zbar:
        m = auxiliary.LSB(rm, c_len - 160)
        return m 
    else:
        return None

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

#Timing of the Fujisaki-Okamoto protocol as implemented with the two different conversion functions        
def time_original_f_o(n, t, m):
    print("Testing the Fujisaki-Okamoto transform with n=", n, "t=", t, "m=", m)
    pk, sk = classic.keygen(n, t, m)
    print("Classic McEliece key generated...")
    k = pk[0].nrows()
    msg = random_matrix(GF(2), 1, k)
    
    num_iter = 10
    duration_sendrier_enc = 0
    duration_sendrier_dec = 0
    duration_ideal_enc = 0
    duration_ideal_dec = 0
    
    for i in range(num_iter):
        if i % (num_iter / 10) == 0:
            print(i, "iterations...")
        #Timing Sendrier
        start_enc = timeit.default_timer()
        c1, c2 = fujisaki_okamoto_encrypt_sendrier(msg, n, k, pk)
        stop_enc = timeit.default_timer()
        d = fujisaki_okamoto_decrypt_sendrier(c1, c2, pk, sk)
        assert d == msg
        stop_dec = timeit.default_timer()
        duration_sendrier_enc += (stop_enc - start_enc)
        duration_sendrier_dec += stop_dec - stop_enc
        
        #Timing Barenghi-Pelosi
        start_enc = timeit.default_timer()
        c1, c2 = fujisaki_okamoto_encrypt_ideal(msg, n, k, pk)
        stop_enc = timeit.default_timer()
        d = fujisaki_okamoto_decrypt_ideal(c1, c2, pk, sk)
        assert d == msg
        stop_dec = timeit.default_timer()
        duration_ideal_enc += (stop_enc - start_enc)
        duration_ideal_dec += (stop_dec - stop_enc)
    print("Average encryption time of Fujisaki-Okamoto with Sendrier's conversion", duration_sendrier_enc / num_iter)
    print("Average decryption time of Fujisaki-Okamoto with Sendrier's conversion", duration_sendrier_dec / num_iter)
    print("Average encryption time of Fujisaki-Okamoto with Barenghi-Pelosi's conversion", duration_ideal_enc / num_iter)
    print("Average decryption time of Fujisaki-Okamoto with Barenghi-Pelosi's conversion", duration_ideal_dec / num_iter)

#Timing of the conversion-free Fujisaki-Okamoto transform    
def time_alt_f_o(n, t, m):
    pk, sk = classic.keygen(n, t, m)
    print("Classic McEliece key generated...")
    k = pk[0].nrows()
    msg = random_matrix(GF(2), 1, k)

    num_iter = 1
    duration_enc = 0
    duration_dec = 0
    
    for i in range(num_iter):
        if i % (num_iter / 10) == 0:
            print(i, "iterations...")
        start_enc = timeit.default_timer()
        c1, c2 = alt_fujisaki_okamoto_encrypt(msg, n, k, pk)
        stop_enc = timeit.default_timer()
        d = alt_fujisaki_okamoto_decrypt(c1, c2, pk, sk)
        assert d == msg
        stop_dec = timeit.default_timer()
        duration_enc += stop_enc - start_enc 
        duration_dec += stop_dec - stop_enc
    print("Average encryption time of Fujisaki-Okamoto without Conversion", duration_enc / num_iter)
    print("Average decryption time of Fujisaki-Okamoto without Conversion", duration_dec / num_iter)

#Timing of the Kobara-Imai alpha protocol
def time_kobara_imai(n, t, m):
    print("Testing the Kobara-Imai alpha transform with n=", n, "t=", t, "m=", m)
    pk, sk = classic.keygen(n, t, m)
    k = pk[0].nrows()
    msg = random_matrix(GF(2), 1, k + 500)
    const = random_matrix(GF(2), 1, 160)
    
    num_iter = 10
    duration_enc = 0
    duration_dec = 0
    
    for i in range(num_iter):
        if i % (num_iter / 10) == 0:
            print(i, "iterations...")
        start_enc = timeit.default_timer()
        c1, c2 = kobara_imai_alpha_encrypt(msg, n, k, pk)
        stop_enc = timeit.default_timer()
        d = kobara_imai_alpha_decrypt(c1, c2, pk, sk)
        assert d == msg
        stop_dec = timeit.default_timer()
        duration_enc += stop_enc - start_enc
        duration_dec += stop_dec - stop_enc
    print("Average encryption time of Kobara-Imai alpha", duration_enc / num_iter)
    print("Average decryption time of Kobara-Imai alpha", duration_dec / num_iter)

#time_original_f_o(1024, 38, 10)
#time_original_f_o(2048, 69, 11)
#time_original_f_o(4096, 128, 12)
#time_alt_f_o(1024, 38, 10)
#time_kobara_imai(1024, 38, 10)
#time_kobara_imai(2048, 69, 11)
#time_kobara_imai(4096, 128, 12)