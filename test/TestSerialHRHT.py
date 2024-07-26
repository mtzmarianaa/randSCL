# Test block HRHT (serial implementation)
import matplotlib.pyplot as plt
from src.hrht import *
from src.hrhtTest import *
import numpy as np
from numpy.linalg import svd, norm
from math import floor
plt.ion()

m = 2**11
n = 2**8
A = createTestMatrices(m, n, type = 7)
nr = 50
skip = 2
rs = np.arange(1, nr, skip)
nyErr = np.empty_like(rs)
nyCond = np.empty_like(rs)

# Vary the sketch size l
ls = [2**8, 2**9, 2**10]
s = 2**4
p = 4

plt.figure()
fig, axs = plt.subplots(1, 3)

for l in ls:
    B1, _, _ = blockHRHT_serial(A, l = l, s = s, p = p)
    _, Ss1, _ = svd(B1)
    axs[0, 0].plot( range(len(Ss1)), Ss1, label = "l = " + str(l)  )
    for i in range(floor(nr/skip)):
        U, Sn = rand_nystrom_svd(A, rs[i], l = l, s = s, p = p)
        nyErr[i] = norm( A - U, 'fro')
        nyCond[i] = Sn[0]/Sn[-1]
    axs[0, 1].plot( rs, nyErr, label = "l = " + str(l) )
    axs[0, 2].plot( rs, nyCond, label = "l = " + str(l) )

fig.suptitle('Effect from changing l')
axs[0, 0].set_title('Singular values of right sketch')
axs[0, 0].set_ylabel('$\sigma_{i}$')
axs[0, 0].set_xlabel('i')
axs[0, 1].set_title("Frobenius error, Nystrom approximation")
axs[0, 0].set_ylabel('Frobenius error')
axs[0, 0].set_xlabel('Rank')
axs[0, 2].set_title("Condition number, Nystrom approximation")
axs[0, 0].set_ylabel('$\kappa$')
axs[0, 0].set_xlabel('Rank')



# Vary s
l = 2**9
ss = [1, 2, 2**2, 2**3, 2**4]
p = 4

plt.figure()
for s in ss:
    B1, _, _ = blockHRHT_serial(A, l = l, s = s, p = p)
    _, Ss1, _ = svd(B1)
    plt.plot( range(len(Ss1)), Ss1, label = "s = " + str(s)  )

plt.title("Singular value of sketches, s changes")
plt.xlabel("i")
plt.ylabel("$\sigma_{i}$")


# Vary p
l = 2**9
ss = 2**3
ps = [1, 2, 4, 8]

plt.figure()
for p in ps:
    B1, _, _ = blockHRHT_serial(A, l = l, s = s, p = p)
    _, Ss1, _ = svd(B1)
    plt.plot( range(len(Ss1)), Ss1, label = "p = " + str(p)  )

plt.title("Singular value of sketches, s changes")
plt.xlabel("i")
plt.ylabel("$\sigma_{i}$")


