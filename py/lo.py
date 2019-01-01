import numpy as np
x = False
y = True
z = True
a = True

for i in range(0,10):
    xi = (not x) or y ^ x or x ^ a
    yi = (not y) or x ^ y or y ^ a
    zi = (not z) or not a ^ z
    ai = (not a) or a ^ x or a ^ y
    x = xi
    y = yi
    z = zi
    a = ai
    print(x,y,z,a)


s=[
        False, False, False,
        False, False, False,
        ]
b=[
        False, False, False,
        False, False, False,
        ]

for i in range(0,20):
    b[0] = (not s[0]) or s[0] ^ s[1] or s[0] ^ s[2] or (s[0] ^ (s[3]))
    b[1] = (not s[1]) or s[1] ^ s[0] or s[1] ^ s[2]
    b[2] = (not s[2]) or s[2] ^ s[1] or s[2] ^ s[0] or (s[2] ^ (s[5]))

    b[3] = (not s[3]) or s[3] ^ (s[4]) or s[3] ^ (s[5]) or (s[3] ^ (s[0]))
    b[4] = (not s[4]) or (s[4]) ^ s[3] or (s[4]) ^ (s[5])
    b[5] = (not s[5]) or (s[5]) ^ s[3] or (s[5]) ^ (s[4]) or (s[5] ^ (s[2]))


    for j in range(6):
        s[j] = b[j]

    r = (not s[0] or not s[1] or not s[2]) and (s[0] or not s[1] or s[2])
    print(not s[0],not s[1], not s[2], s[3], not s[4], s[5],  r)
    print("a = %s, b = %s, c = %s" %(s[0], s[1], s[2]))
    if r:
        break;



print('---------------------')


s=[
        False, False, False,
        False, False, False,
        False, False, False,
        False, False, False,
        ]
b=[
        False, False, False,
        False, False, False,
        False, False, False,
        False, False, False,
        ]

for i in range(0,10):
    b[0] = (not s[0]) or s[0] ^ s[1] or s[0] ^ s[2] or (s[0] ^ (not s[6])) # a
    b[1] = (not s[1]) or s[1] ^ s[0] or s[1] ^ s[2] or (s[1] ^ (not s[10])) or (s[1] ^ (not s[3])) # b
    b[2] = (not s[2]) or s[2] ^ s[1] or s[2] ^ s[0] or (s[2] ^ (not s[4])) # c

    b[3] = (not s[3]) or s[3] ^ (s[4]) or s[3] ^ (s[5]) or (s[3] ^ s[1]) # -b
    b[4] = (not s[4]) or (s[4]) ^ s[3] or (s[4]) ^ (s[5]) or (s[4] ^ s[2]) or (s[4] ^ s[7]) # -c
    b[5] = (not s[5]) or (s[5]) ^ s[3] or (s[5]) ^ (s[4]) or (s[5] ^ (s[8])) # -d

    b[6] = (not s[6]) or s[6] ^ (s[7]) or s[6] ^ (s[8]) or (s[6] ^ (s[0])) or (s[6] ^ s[9]) # -a
    b[7] = (not s[7]) or (s[7]) ^ s[8] or (s[7]) ^ (s[6]) or (s[7] ^ (not s[4]))  # c
    b[8] = (not s[8]) or (s[8]) ^ s[7] or (s[8]) ^ (s[6]) or (s[8] ^ (not s[5])) or (s[8] ^ (not s[11])) # d

    b[9] = (not s[9]) or s[9] ^ (s[10]) or s[9] ^ (s[11]) or (s[9] ^ (not s[6])) # a
    b[10] = (not s[10]) or (s[10]) ^ s[11] or (s[10]) ^ (s[9]) or (s[10] ^ s[1])  # -b
    b[11] = (not s[11]) or (s[11]) ^ s[10] or (s[11]) ^ (s[9]) or (s[11] ^ (s[8])) # -d


    for j in range(12):
        s[j] = b[j]

    #r = (s[0] or s[1] or s[2]) and (not s[1] or not s[2] or not s[3]) and (not s[0] or s[2] or s[3]) and (s[0] or not s[1] or not s[3])
    r = (s[0] or s[1] or s[2]) and (not s[3] or not s[4] or not s[5]) and (not s[6] or s[7] or s[8]) and (s[9] or not s[10] or not s[11])
    print(s[0],s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11],r)
    _a = s[0]
    _b = s[1]
    _c = s[2]
    _d = s[8]
    r = (_a or _b or _c) and (not _b or not _c or not _d) and (not _a or _c or _d) and (_a or not _b or not _d)
    print("a = %s, b = %s, c = %s, d = %s" %(s[0], s[1], s[2], s[8]))
    print("a = %s, b = %s, c = %s, d = %s" %(_a, _b, _c, _d))
    if r:
        break;
