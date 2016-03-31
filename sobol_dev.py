import numpy
from _tgo import *
def frac_dec2bin(float):  # Convert frac decimals to their bin strings (w/ no dec)
    if float == 0:
        return '0'

    decp = 2**(-1)
    fstr = ""
    while float > 0 or decp >= 1:
        if float < 1:
            if len(fstr) >= 1000: # 10 fractional digits max
                   break
        if float >= decp:
            fstr += "1"
            float -= decp
        else:
            fstr += "0"
        decp /= 2
    return fstr





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

k =6
m_6 =(  ((2**(s - 2)) * m_i[k - s + 2])
        * a_i[1]
      ^ ((2**(s - 1)) * m_i[k - s + 1])
        * a_i[2]
      ^ ((2**(s    )) * m_i[k - s])
      ^                 m_i[k - s]
     )

#NOTE: The form of the recurrance relation does not change for higher
# values of k

m_i.append(m_6)

print 'm_4 = '
print m_4

print 'm_5 = '
print m_5

for k in range(1, 6):
    pass
    print 'v_{}_int = '.format(k)
    v_int = m_i[k] / (2.0 ** k)
    print v_int
    print 'v_{}_bin * (2.0 ** k) = '.format(k)
    #print bin(v_int)
    print frac_dec2bin(v_int)
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



import numpy
N = 10
x  = numpy.zeros(N)
x[0] = 0.0  # Define x_j with x_0_j = 0
c_i = numpy.zeros(N)

print 'x_{} = {} '.format(0, 0.0)


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
    if i > 0:
        k = int(c_i[i - 1])
        v_k_j =  m_i[k] / (2.0 ** k)
        v_bin_s = frac_dec2bin(v_k_j)
        x_bin_s = frac_dec2bin(x[i - 1])
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

        v_bin = int(v_bin_s, 2)
        x_bin = int(x_bin_s, 2)
    #    print 'i = {}'.format(i)
     #   print 'v_bin_s = {}'.format(v_bin_s)
     #   print 'x_bin_s = {}'.format(x_bin_s)
        x_xor = x_bin ^ v_bin
     #   print 'x_xor_b = {}'.format(bin(x_xor)[2:])
     #   print 'x_xor   = {}'.format(x_xor)
        #print x_xor
        x[i] =  x_xor/ 2.0**(lmax)
        x[i] = float(x[i])

        print 'x_{} = {} '.format(i, x[i])

