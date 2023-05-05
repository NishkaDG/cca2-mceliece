'''
This file contains miscellaneous conversion functions, hash functions, and other helper functions for the other programs.
'''

from sage.all_cmdline import *   # import sage library
from Crypto.Hash import SHAKE256
from math import ceil, floor, log2
from random import randrange

#Takes an integer input num and a bitlength length
#Returns a vector (as a list) of the binary representation of num in exactly length bits
#Truncates or pads as needed
def pad_as_list(num, length):
    res = []
    bin_val = bin(num)[2:]
    if len(bin_val) < length:
        for i in range(length - len(bin_val)):
            res.append(0)
    else:
        bin_val = bin_val[:length]
    for i in bin_val:
        res.append(int(i))
    #print(res)
    return res

#Takes an integer input num and a bitlength length
#Returns a string of bits of the binary representation of num in exactly length bits 
#Truncates or pads as needed
def pad_as_bitstring(num, length):
    res = ''
    bin_val = bin(num)[2:]
    if len(bin_val) < length:
        for i in range(length - len(bin_val)):
            res = res + '0'
    else:
        bin_val = bin_val[:length]
    for i in bin_val:
        res = res + i 
    return res

#Takes inputs two vectors (sagemath 1 * ncols matrix) 
#Returns a bytearray of the concatenation of these vectors
def concat_vectors_to_bytearray(vec1, vec2):
    cc = bytearray(b'\x00')
    for i in vec1[0]:
        cc.append(i)
    for i in vec2[0]:
        cc.append(i)
    return cc

#Takes as input a single vector (sagemath 1 * ncols matrix) 
#Returns the bytearray conversion of the vector   
def vector_to_bytes(vec):
    cc = bytearray(b'\x00')
    for i in vec[0]:
        cc.append(i)
    return cc

#Takes as inputs two vectors (sagemath 1 * ncols matrices) 
#Returns their concatenation as a sagemath matrix
def concat_vectors(v1, v2):
    assert v1.nrows() == 1 and v2.nrows() == 1
    new_len = v1.ncols() + v2.ncols()
    res = matrix(GF(2), 1, new_len)
    len_1 = v1.ncols()
    for i in range(new_len):
        if i < len_1:
            res[0, i] = v1[0, i]
        else:
            res[0, i] = v2[0, i - len_1]
    return res

#Takes as input a vector (sagemath row matrix)
#Returns the bitstring (elements of the vector concatenated into a string)    
def vector_to_bitstring(vec):
    B = ''
    for bit in vec[0]:
        B = B + str(bit)
    return B
    
#Takes as input bitstring
#Returns the bitstring as a vector (sagemath row matrix)
def bitstring_to_vector(bitstring):
    n = len(bitstring)
    vec = matrix(GF(2), 1, n)
    for i in range(n):
        vec[0, i] = int(bitstring[i])
    return vec

#Takes an input bitstring 
#Returns the bitstring as a bytearray    
def bitstring_to_bytes(bitstring):
    cc = bytearray(b'\x00')
    for i in bitstring:
        cc.append(int(i))
    return cc

#Takes an input bitstring 
#Returns the run-length encoding of the bitstring, i.e, the number of consecutive 0s preceding each occurrence of 1, as a list
def bitstring_to_positional(bitstring):
    delta_lst = []
    ctr = 0
    for ele in bitstring:
        if int(ele) == 1:
            delta_lst.append(ctr)
            ctr = 0
        else:
            ctr = ctr + 1
    return delta_lst

#Takes an input sagemath row matrix vec 
#Returns the run-length encoding of the bitstring as above
def vector_to_positional(vec):
    delta_lst = []
    ctr = 0
    for ele in vec[0]:
        if int(ele) == 1:
            delta_lst.append(ctr)
            ctr = 0
        else:
            ctr = ctr + 1
    return delta_lst

#Takes as input a bitlength n and list delta_lst of run-length encodings i.e, the number of consecutive 0s preceding each occurrence of 1
#Returns the corresponding bitstring 
def positional_to_bitstring(delta_lst, n):
    print('res', delta_lst)
    bitstring = ''
    delpos = 0
    ctr = 0
    for i in range(n):
        if ctr == delta_lst[delpos]:
            bitstring = bitstring + '1'
            delpos = delpos + 1
            ctr = 0
        else:
            bitstring = bitstring + '0'
            ctr = ctr + 1
    assert len(bitstring) == n        
    return bitstring

#Takes as input a bitlength n and list delta_lst of run-length encodings as above 
#Returns a sagemath row matrix of the corresponding bitstring. 
def positional_to_vector(delta_lst, n):
    z = matrix(GF(2), 1, n)
    ctr = 0
    for d in delta_lst:
        ctr = ctr + d
        z[0, ctr] = 1
        ctr = ctr + 1
    return z        

#Takes as input a vector (sagemath row matrix) vec and an integer x 
#Returns a new vector (sagemath row matrix) consisting of the x least significant bits of vec    
def LSB(vec, x):
    assert vec.ncols() >= x
    res = matrix(GF(2), 1, x)
    start = vec.ncols() - x 
    for i in range(x):
        res[0, i] = vec[0, i + start]
    return res 

#Takes as input a vector (sagemath row matrix) vec and an integer x 
#Returns a new vector (sagemath row matrix) consisting of the x most significant bits of vec        
def MSB(vec, x):
    assert vec.ncols() >= x 
    res = matrix(GF(2), 1, x)
    for i in range(x):
        res[0, i] = vec[0, i]
    return res

#Takes as input a bitstring (as bytes), n, t 
#Returns the SHAKE256 XOF output of the bitstring in C(n,t) bits as an integer
def H(bitstring, n, t):
    nct = int(ceil(factorial(n)/(factorial(t) * factorial(n-t))))
    bitlength = int(ceil(log(nct, 2).n()))
    bytelength = int(ceil(bitlength / 8))
    
    shake = SHAKE256.new()
    shake.update(bitstring)
    h = shake.read(bytelength).hex()
    h = int(h, 16) % nct
    #print(nct)
    return h

#Takes as input a bitstring (as bytes), and a bitlength k to hash to 
#Returns the SHAKE256 XOF output of the bitstring in k bits as a bitstring 
def H1(bitstring, k):
    bytelength = int(ceil(k / 8))
    
    shake = SHAKE256.new()
    shake.update(bitstring)
    h = shake.read(bytelength).hex()
    h = int(h, 16) % (2 ** k)
    binh = pad_as_bitstring(h, k)
    return binh

#Takes as input a bitstring r (as bytes) and a bitlength k to hash to 
#Returns the SHAKE256 XOF output of the bitstring in k bits as a sagemath matrix 
def R(r, k):
    bytelength = int(ceil(k / 8))
    shake = SHAKE256.new()
    shake.update(r)
    h_hex = shake.read(bytelength).hex()
    h_list = pad_as_list(int(h_hex, 16), k)
    h_vec = matrix(GF(2), h_list)
    return h_vec

def test_positional_vector_interconversion():
    num_iter = 10000
    
    for i in range(num_iter):
        n = randrange(2, 200)
        t = randrange(n//2)
        lv = []
        lim = n-t
        for i in range(t):
            lv.append(randrange(lim - sum(lv)))
        v = positional_to_vector(lv, n)
        lv2 = vector_to_positional(v)
        if not (lv == lv2):
            print("Failure!", v, lv, lv2)
            break
    
#test_positional_vector_interconversion()