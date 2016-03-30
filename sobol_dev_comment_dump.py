import numpy
from _tgo import *
#d       s       a       m_i
#2       1       0       1
#3       2       1       1 3
d = 3
s = 2
a = 1
# k =  1  2
m_i = (1, 3)


#>>> bin(1)
#'0b1'
#>>> bin(1)[2:]
#'1'
#m_k_j
m_1_j = m_i[0]
m_2_j = m_i[1]
m_3_j = (2 ** (s - 1)) * m_i[1] ^ ((2 ** s)) * m_i[0]
m_4_j = (2 ** (s - 1)) * m_i[1] ^ ((2 ** s)) * m_i[0] ^ m_3_j

# m_k_j =   (2**(0        )) * m_i[k - 1]       * a_1
#         ^ (2**(2        )) * m_i[k - 2]       * a_2
#         ...
#         ^ (2**(s - 1    )) * m_i[k - s + 1]   * a_s-1
#         ^ (2**(s        )) * m_i[k - s]
#         ^ m[k - s]

# v_k_j = m_k_j / 2^k


# Example Kuo's notes # x^3 + 0 x^2 + 1 x + 1
s = 3
#  a_1 a_2
a = [0, 1]
# k =  0,    1  2  3
m_i = [None, 1, 3, 7]
a_i = [None, 0, 1]
# Find m_4
k = 4
m_4 =(  ((2**(s - 2)) * m_i[k - s + 2])  # k - s + 2 = k - 1 = 3
        * a_i[1]
      ^ ((2**(s - 1)) * m_i[k - s + 1])  # k - s + 1 = k - 2 = 2
        * a_i[2]                         # s - 1 = 2 --> a_s-1
      ^ ((2**(s    )) * m_i[k - s])
      ^                 m_i[k - s]
     )

m_i.append(m_4)

k = 5
m_5 =(  ((2**(s - 2)) * m_i[k - s + 2])
        * a_i[1]
      ^ ((2**(s - 1)) * m_i[k - s + 1])
        * a_i[2]
      ^ ((2**(s    )) * m_i[k - s])
      ^                 m_i[k - s]
     )

#NOTE: The form of the recurrance relation does not change for higher
# values of k

m_i.append(m_5)

print 'm_4 = '
print m_4

print 'm_5 = '
print m_5

for k in range(1, 6):
    pass
    print 'v_{}_int = '.format(k)
    v_int = m_i[k]# / (2.0 ** k)
    print v_int
    print 'v_{}_bin * (2.0 ** k) = '.format(k)
    print bin(v_int)
    #print m_i[k]
    #print bin( v_int )
    #print m_i[k] / (2.0 ** k)


# Test XOR operations with fractions
# (0.1)_2 XOR (0.11)_2 = (0.01)_2 = 0.25  (= 1/ 2.0**2)
# 1 XOR 3 = 2
# 2 / 2.0 ** 3 = 0.25


# (0.1)_2 XOR (0.111)_2 = (0.011)_2 = 0.25  (= 1/ 2.0**2)
# 1 XOR 7 = 6 / 2.0**4 = 0.375
# 2 / 2.0 ** 3 = 0.25


print '='*100
def frac_dec2bin(f):  # Convert frac decimals to their bin strings (w/ no dec)
    if f == 0:
        return '0'

    ig = 2**(-1)
    st = ""
    while f > 0 or ig >= 1:
        if f < 1:
            if len(st) >= 1000: # 10 fractional digits max
                   break
        if f >= ig:
            st += "1"
            f -= ig
        else:
            st += "0"
        ig /= 2
    return st


# Define x_j with x_0_j = 0
import fractions
#x = [0.0]
import numpy
N = 6
x  = numpy.zeros(N)
x[0] = 0.0
#x_int.append(fractions.Fraction(x[0]))
#x_int = [0, 0, 0, 0, 0, 0]
#x_int[0] = fractions.Fraction(x[0])

x_int = numpy.zeros(N)
#c_i = []
c_i = numpy.zeros(N)
for i in range(0, N):
    #print 'i = {}'.format(i)
    bin_string = bin(i)  # Binary string representation of i
    #print 'bin i = {}'.format(bin_string)
    #print 'bin i [2:]= {}'.format(bin_string[2:])
    try:
        c = bin_string[2:][::-1].index('0') + 1
    except(ValueError):  # No substring containing 0
        c = len(bin_string[2:]) + 1

    #print 'c_{} = {}'.format(i, c)
    c_i[i] = int(c) # TODO: Add check

    #print 'm_[c_({} - 1)] = {}'.format(i, m_i[c_i[i - 1]])

    if i > 0:
        k = int(c_i[i - 1])
        #x_new_int = x_int[i - 1].numerator ^ m_i[k]
        #x_new = ((x_new_int) / float(x_int[i].denominator)) / (2.0 ** k)
        #x_bin = bin(x_int[i - 1].numerator)
        #m_bin = bin(m_i[k])
        #x_bin ^ m_bin
        #x_bin ^ m_bin

        v_k_j =  m_i[k] / (2.0 ** k)
        #print 'v_{}_j = {}'.format(k, v_k_j)

        #print str(x[i - 1])
        #print len(str(x[i - 1])[2:])

    #    fracxl = len(str(x[i - 1])[2:])
    #    fracvl = len(str(v_k_j)[2:])
    #    x_10 = x[i - 1]*10.0**fracxl
        # print x_10
    #    v_10 = v_k_j*10.0**fracvl
        # print v_10
        #maxf = max(fracxl, fracvl)
        #print 'maxf = {}'.format(maxf)
        #x_10 = x[i - 1]*10.0**maxf
        #print 'x_10= {}'.format(x_10)
        #v_10 = v_k_j*10.0**maxf
        #print 'v_10= {}'.format(v_10)
    #    x_new_int = int(x_10) ^ int(v_10)
    #    x_new = float(x_new_int) / 10.0**(fracxl+fracvl)
        #x_new = float(x_new_int) / float(10.0**maxf)
    #    x[i] = x_new
        #x_int.append(fractions.Fraction(x[i]))

        #x[i] = m_i[k]


        v_bin_s = frac_dec2bin(v_k_j)
        #print 'v_bin_s = {}'.format(v_bin_s)
        x_bin_s = frac_dec2bin(x[i - 1])
        #print 'x_bin_s = {}'.format(x_bin_s)
        #print x_bin_s
        lv = len(v_bin_s)
        lx = len(x_bin_s)
        if lx > lv:
            v_bin_s += '0' * (lx - lv)
            lmax = lx
        elif lv > lx:
            x_bin_s += '0' * (lv - lx)
            lmax = lv
        else:
            lmax = lx

        #print 'v_bin_s fixed = {}'.format(v_bin_s)
        #print 'x_bin_s fixed = {}'.format(x_bin_s)

        v_bin = int(v_bin_s, 2)
        x_bin = int(x_bin_s, 2)
        x_xor = x_bin ^ v_bin
        print x_xor
        x[i] = x_xor/ 2.0**(len(v_bin_s))
        print 'x_{} = {} '.format(i, x[i])




import struct
def binary(num):
    return ''.join(bin(ord(c)).replace('0b', '').rjust(8, '0') for c in struct.pack('!f', num))




# Use 2 to find m_4 = 5 and m_5 = 7

# x_0_k = 0
# x_i_j = x_[i - 1]_j ^ v_[c_{i - 1}, j]
# See notes
# c_i is the first 0 digit from the right of the binary representation
#i = (...i_3 i_2 i_1 )_2

# Practice, shift indices back later
# m_i = [None, 1, 3]
#
#
# print m_3_j
# print m_4_j
#
# k = 1
# print m_i[0] / 2.0
# k = 2
# print m_i[1] / (2.0 ** k)
# k = 3
# print m_3_j / (2.0 ** k)
# k = 3
# print m_4_j / (2.0 ** k)

# Bratley and Fox test:
# x^3 + 0 x^2 + 1 x + 1
# NOTE SHIFT a_i one to the right for Kuo's convention
# a1 = 0
# a2 = 1
#m_1_j = 1
#m_2_j = 3
#m_3_j = 7
#m_k_j =