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

# Define x_j with x_0_j = 0
x = [0]
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