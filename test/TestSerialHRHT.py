# Test block HRHT (serial implementation)
import matplotlib.pyplot as plt
from src.hrht import *
from src.testSingVHRHT import *
import numpy as np
from numpy.linalg import svd
plt.ion()

m = 2**11
n = 2**8
V1 = createTestMatrices(m, n, type = 1)
V2 = createTestMatrices(m, n, type = 2)
V3 = createTestMatrices(m, n, type = 3)

# Vary the sketch size l
ls = [2**8, 2**9, 2**10]
s = 2**4
p = 4

plt.figure()
for l in ls:
    B1, _, _ = blockHRHT_serial(V1, l = l, s = s, p = p)
    _, Ss1, _ = svd(B1)
    plt.plot( range(len(Ss1)), Ss1, label = "SVDG, l = " + str(l)  )
    B2, _, _ = blockHRHT_serial(V2, l = l, s = s, p = p)
    _, Ss2, _ = svd(B2)
    plt.plot( range(len(Ss2)), Ss2, label = "colHad, l = " + str(l)  )
    B3, _, _ = blockHRHT_serial(V3, l = l, s = s, p = p)
    _, Ss3, _ = svd(B3)
    plt.plot( range(len(Ss3)), Ss3, label = "rowHad, l = " + str(l)  )

plt.title("Singular value of sketches, l changes")
plt.xlabel("i")
plt.ylabel("$\sigma_{i}$")


# Vary s
l = 2**9
ss = [1, 2, 2**2, 2**3, 2**4]
p = 4

plt.figure()
for s in ss:
    B1, _, _ = blockHRHT_serial(V1, l = l, s = s, p = p)
    _, Ss1, _ = svd(B1)
    plt.plot( range(len(Ss1)), Ss1, label = "SVDG, s = " + str(s)  )
    B2, _, _ = blockHRHT_serial(V2, l = l, s = s, p = p)
    _, Ss2, _ = svd(B2)
    plt.plot( range(len(Ss2)), Ss2, label = "colHad, s = " + str(s)  )
    B3, _, _ = blockHRHT_serial(V3, l = l, s = s, p = p)
    _, Ss3, _ = svd(B3)
    plt.plot( range(len(Ss3)), Ss3, label = "rowHad, s = " + str(s)  )

plt.title("Singular value of sketches, s changes")
plt.xlabel("i")
plt.ylabel("$\sigma_{i}$")


# Vary p
l = 2**9
ss = 2**3
ps = [1, 2, 4, 8]

plt.figure()
for p in ps:
    B1, _, _ = blockHRHT_serial(V1, l = l, s = s, p = p)
    _, Ss1, _ = svd(B1)
    plt.plot( range(len(Ss1)), Ss1, label = "SVDG, p = " + str(p)  )
    B2, _, _ = blockHRHT_serial(V2, l = l, s = s, p = p)
    _, Ss2, _ = svd(B2)
    plt.plot( range(len(Ss2)), Ss2, label = "colHad, p = " + str(p)  )
    B3, _, _ = blockHRHT_serial(V3, l = l, s = s, p = p)
    _, Ss3, _ = svd(B3)
    plt.plot( range(len(Ss3)), Ss3, label = "rowHad, p = " + str(p)  )

plt.title("Singular value of sketches, s changes")
plt.xlabel("i")
plt.ylabel("$\sigma_{i}$")


