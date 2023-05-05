'''
This file contains an implementation of Classic McEliece PKC as well as Goppa Code initialization via Bernstein.
References:
 - Daniel J. Bernstein. Understanding binary-Goppa decoding. Cryptology ePrint Archive, Paper 2022/473. https://eprint.iacr.org/2022/473. 2022. url: https://eprint.iacr.org/2022/473
 - D. Engelbert, R. Overbeck, and A. Schmidt. A Summary of McEliece-Type Cryptosystems and their Security. Cryptology ePrint Archive, Paper 2006/162. https://eprint.iacr.org/2006/162. 2006. url: https://eprint.iacr.org/2006/162.
'''

from sage.all_cmdline import *   # import sage library
import bernstein
import timeit

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_38 = Integer(38); _sage_const_6 = Integer(6); _sage_const_5 = Integer(5); _sage_const_69 = Integer(69); _sage_const_128 = Integer(128); _sage_const_7 = Integer(7); _sage_const_0 = Integer(0); _sage_const_1024 = Integer(1024); _sage_const_10 = Integer(10); _sage_const_2048 = Integer(2048); _sage_const_11 = Integer(11); _sage_const_4096 = Integer(4096); _sage_const_12 = Integer(12)

def generate_P(n):
    R = GF(2)
    M = identity_matrix(R, n)
    perm = Permutations(n).random_element()
    for i in range(n):
        j = perm[i] - 1
        M.swap_rows(i, j)
    return M

def generate_S(k):
    R = GF(2)
    M = random_matrix(R, k, k)
    while M.is_singular():
        M = random_matrix(R, k, k)
    return M

def generate_G_squarefree(n, t, m):
    q = 2**m 
    F = GF(q)
    Fpoly = F['x']
    (x,) = Fpoly._first_ngens(1)
    a = list(F)
    while True:
        shuffle(a)
        L = a[:n]
        g = Fpoly([F.random_element() for j in range(t)] + [1])
        if g.is_squarefree():
            if all(g(aj) != 0 for aj in L):
                break
    C = codes.GoppaCode(g, L)
    G = C.generator_matrix()
    k = G.nrows()
    return (k, G, g, L, F) 
    
def generate_G_irreducible(n, t, m):
    Fp = GF(_sage_const_2 )
    Fpm = GF(_sage_const_2 **m)

    R = Fpm['x']; (x,) = R._first_ngens(1)
    g = x**t + _sage_const_1 
    if t == _sage_const_38 :
        g = g + x**_sage_const_6  + x**_sage_const_5  + x
    elif t == _sage_const_69 :
        g = g + x**_sage_const_6  + x**_sage_const_5  + x**_sage_const_2 
    elif t == _sage_const_128 :
        g = g + x**_sage_const_7  + x**_sage_const_2  + x
    else:
        print("Undefined behaviour!")

    L = []
    ctr = _sage_const_0 

    while ctr < n:
        y = Fpm.random_element()
        assert(g(y) != _sage_const_0 )
        if y not in L:
            L.append(y)
            ctr = ctr + _sage_const_1 

    C = codes.GoppaCode(g, L)
    G = C.generator_matrix()
    k = G.nrows()

    assert k >= (n - m*t)

    return (k, G, g, L, Fpm)
    
def keygen(n, t, m):
    goppa_info = generate_G_squarefree(n, t, m)
    k = goppa_info[0]
    G1 = goppa_info[1]
    decoding_info = (goppa_info[2], goppa_info[3], goppa_info[4])
    P = generate_P(n)
    S = generate_S(k)
    
    G = S * G1 * P
    pk = (G, t)
    sk = (S, P, decoding_info)
    return pk, sk
    
def encrypt(m, z, pk):
    G = pk[0]
    c = (m * G) + z
    return c
    
def decrypt(c, sk, pk):
    # c = mSGP + e 
    # do cP^{-1} = mSG + eP^{-1} (see Ivan mail on why this is important)
    # do Bernstein error correcting to remove eP^{-1} (P is a permutation matrix, so this term is also a vector of weight t)
    # now do SG.solve_left(cP^{-1}) to get m
    S = sk[0]
    P = sk[1]
    decoding_info = sk[2]
    SGP = pk[0] # this is SG'P where G is the generator matrix
    n = SGP.ncols()
    k = SGP.nrows()
    t = pk[1]
    g = decoding_info[0]
    alpha = decoding_info[1]
    F = decoding_info[2]
    
    P1 = P.inverse()
    SG = SGP * P1
    c = c * P1 #now we have c = mSG + eP^{-1}
    
    e_list = bernstein.goppa_errors(n, t, F, alpha, g, c[0])
    eP = matrix(GF(2), 1, n, [e_list]) #Remember that we multiplied with P^{-1} so the error that we corrected is not the original error e 
    e = eP * P
    
    c = c + eP #now we have c = mSG = (mS)(G)
    m = SG.solve_left(c)
    
    return m, e
    
def select_error(z, t, n):
    RR = Integers(n)
    wt = vector(z).hamming_weight()
    while not wt == t:
        pos_to_change = RR.random_element(n)
        if wt < t:
            z[0, pos_to_change] = 1
        elif wt > t:
            z[0, pos_to_change] = 0
        wt = vector(z).hamming_weight()
    
def test_keygen(n, t, m):
    num_iter = 10000
    duration = 0
    for i in range(num_iter):
        start = timeit.default_timer()
        pk, sk = keygen(n, t, m)
        stop = timeit.default_timer()
        duration += stop - start
    print("Average time for McEliece key generation is", duration / num_iter)
    
def test_encrypt(n, t, m):
    pk, sk = keygen(n, t, m)
    print("Keygen done.")
    k = pk[0].nrows()
    msg = random_matrix(GF(2), 1, k)
    
    num_iter = 10000
    duration = 0
    for i in range(num_iter):
        if (i % (num_iter / 10)) == 0:
            print(i, "iterations...")
        start = timeit.default_timer()
        z = random_matrix(GF(2), 1, n)
        select_error(z, t, n)
        c = encrypt(msg, z, pk)
        stop = timeit.default_timer()
        duration = duration + stop - start
    print("Average encryption time of classic McEliece (including error vector generation) is", duration / num_iter)
    
def test_decrypt(n, t, m):
    print("testing for n=", n, "t=", t, "m=", m)
    pk, sk = keygen(n, t, m)
    print("Keygen done")
    k = pk[0].nrows()
    num_iter = 10
    duration = 0
    
    for i in range(num_iter):
        if (i % (num_iter / 10)) == 0:
            print(i, "iterations...")
        msg = random_matrix(GF(2), 1, k)
        z = matrix(GF(2), 1, n)
        select_error(z, t, n)
        c = encrypt(msg, z, pk)
        start = timeit.default_timer()
        d, e = decrypt(c, sk, pk)
        assert d == msg
        stop = timeit.default_timer()
        duration += stop - start
    print("Average decryption time of classic McEliece is", duration / num_iter)
    
#test_decrypt(1024, 38, 10)
#test_decrypt(2048, 69, 11)
#test_decrypt(4096, 128, 12)