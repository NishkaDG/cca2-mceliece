'''
This file contains an implementation of Classic McEliece PKC
'''

#TODO: What is the relationship between t and k
#Actual goppa decryption

from sage.all_cmdline import *   # import sage library
import timeit

_sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_38 = Integer(38); _sage_const_6 = Integer(6); _sage_const_5 = Integer(5); _sage_const_69 = Integer(69); _sage_const_128 = Integer(128); _sage_const_7 = Integer(7); _sage_const_0 = Integer(0); _sage_const_1024 = Integer(1024); _sage_const_10 = Integer(10); _sage_const_2048 = Integer(2048); _sage_const_11 = Integer(11); _sage_const_4096 = Integer(4096); _sage_const_12 = Integer(12)

def generate_P(n):
    R = Integers(2)
    M = identity_matrix(R, n)
    perm = Permutations(n).random_element()
    for i in range(n):
        j = perm[i] - 1
        M.swap_rows(i, j)
    return M

def generate_S(k):
    R = Integers(2)
    M = random_matrix(R, k, k)
    while M.is_singular():
        M = random_matrix(R, k, k)
    return M

def generate_G_squarefree(n, t, m):
    q = 2**m 
    k = GF(q)
    kpoly = k['x']
    (x,) = kpoly._first_ngens(1)
    a = list(k)
    while True:
        shuffle(a)
        L = a[:n]
        coeffs = ([0] * (t - 1))
        flip = 0
        while flip < 3:
            pos = randrange(t - 1)
            if coeffs[pos] == 0:
                coeffs[pos] = k.random_element()
                flip = flip + 1
        g = kpoly([k.random_element() for j in range(t)] + [1])
        g = kpoly([1] + coeffs + [1])
        if g.is_squarefree():
            if all(g(aj) != 0 for aj in L):
                break
    #print(g, g.is_irreducible())
    C = codes.GoppaCode(g, L)
    G = C.generator_matrix()
    H = C.parity_check_matrix()
    k = G.nrows()
    #print(k)
    return G, k
    
def generate_G_irreducible(n, t, m):
    Fp = GF(_sage_const_2 )
    Fpm = GF(_sage_const_2 **m)

    R = Fpm['x']; (x,) = R._first_ngens(1)
    #K.<a> = GF(2**m, name = 'a', modulus="minimal_weight")
    g = x**t + _sage_const_1 
    if t == _sage_const_38 :
        g = g + x**_sage_const_6  + x**_sage_const_5  + x
    elif t == _sage_const_69 :
        g = g + x**_sage_const_6  + x**_sage_const_5  + x**_sage_const_2 
    elif t == _sage_const_128 :
        g = g + x**_sage_const_7  + x**_sage_const_2  + x
    else:
        print("Undefined behaviour!")
    #print(g, g.is_squarefree())

    L = []
    ctr = _sage_const_0 

    while ctr < n:
        y = Fpm.random_element()
        assert(g(y) != _sage_const_0 )
        if y not in L:
            L.append(y)
            ctr = ctr + _sage_const_1 

    #print(len(L))
    #print(type(g))

    C = codes.GoppaCode(g, L)
    G = C.generator_matrix()
    H = C.parity_check_matrix()
    #print(G.nrows(), 'x', G.ncols())
    #print(H.nrows(), 'x', H.ncols())

    k = G.nrows()

    assert k >= (n - m*t)

    #E = codes.encoders.GoppaCodeEncoder(C)
    #D = codes.decoders.GoppaCodeDecoder(C)

    #E = codes.encoders.GoppaCodeEncoder(C)
    #D = C.decoder()

    #msg = random_matrix(GF(2), 1, k)[0]
    #print(msg)

    #c = E.encode(msg)
    #w = D.decode_to_message(c) 

    return G, k
    
def keygen(n, t, m):
    G1, k = generate_G_squarefree(n, t, m)
    P = generate_P(n)
    S = generate_S(k)
    
    G = S * G1 * P
    pk = (G, t)
    sk = (S, P, k)
    return pk, sk
    
def encrypt(m, z, pk):
    #print(m, m.parent())
    G = pk[0]
    #print(G)
    c = (m * G) + z
    return c
    
def decrypt(c, sk):
    # c = mSGP + e 
    # do cP^{-1} = mSG + eP^{-1} (see Ivan mail on why this is important)
    # do Bernstein error correcting to remove eP^{-1} (P is a permutation matrix, so this term is also a vector of weight t)
    # now do SG.solve_left(cP^{-1}) to get m
    return None
    
def select_error(z, t, n):
    RR = Integers(n)
    wt = vector(z).hamming_weight()
    while not wt == t:
        #print("weight of z is", wt)
        pos_to_change = RR.random_element(n)
        if wt < t:
            z[0, pos_to_change] = 1
        elif wt > t:
            z[0, pos_to_change] = 0
        wt = vector(z).hamming_weight()
    #return z
    
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
    k = sk[2]
    m = random_matrix(Integers(2), 1, k)
    
    num_iter = 10000
    duration = 0
    for i in range(num_iter):
        start = timeit.default_timer()
        z = random_matrix(Integers(2), 1, n)
        #print(vector(z).hamming_weight())
        select_error(z, t, n)
        #print(vector(z).hamming_weight())
        #print("Error vector selected.")
        c = encrypt(m, z, pk)
        #print("Encryption done.")
        stop = timeit.default_timer()
        duration = duration + stop - start
    print("Average time for encryption (including error vector generation) is", duration / num_iter)
    
#generate_G_squarefree(1024, 38, 10)
#generate_G_irreducible(1024, 38, 10)
#generate_G_irreducible(2048, 69, 11)
#generate_G_irreducible(4096, 128, 12)
test_keygen(1024, 38, 10)
