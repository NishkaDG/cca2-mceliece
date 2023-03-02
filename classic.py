#TODO: What is the relationship between t and k
#Actual goppa decryption + generation (sage has a goppa library that can be used as placeholder)

from sage.all_cmdline import *   # import sage library

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
    
def generate_G(k, n):
    R = Integers(2)
    M = random_matrix(R, k, n)
    return M
    
def keygen(n, t, k):
    P = generate_P(n)
    S = generate_S(k)
    G1 = generate_G(k, n)
    
    G = S * G1 * P
    pk = (G, t)
    sk = (S, P)
    return pk, sk
    
def encrypt(m, z, pk):
    #print(m, m.parent())
    G = pk[0]
    #print(G)
    c = (m * G) + z
    return c
    
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
    
def test():
    n = 1024
    k = 524
    t = 50
    pk, sk = keygen(n, t, k)
    print("Keygen done.")
    m = random_matrix(Integers(2), 1, k)
    z = random_matrix(Integers(2), 1, n)
    print(vector(z).hamming_weight())
    select_error(z, t, n)
    print(vector(z).hamming_weight())
    print("Error vector selected.")
    c = encrypt(m, z, pk)
    print("Encryption done.")
    
#test()