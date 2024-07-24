from src.hrhtTest import *
from matplotlib.colors import ListedColormap
import colorcet as cc
colormap2 = "cet_linear_worb_100_25_c53_r"

# SCRIPT TO TEST THE SINGULAR VALUES OF BLOCK HRHT

testType = 2
n = 2**8
m = 2**12
p = 2**1
r = int(m/p)
V = createTestMatrices(m, n, type = testType)

lTest = np.arange(n, int(m/p), 12)
sTest = np.arange(1, n+1, 4)

imErr2 = np.zeros((len(lTest), len(sTest)))
imErrMax = np.zeros((len(lTest), len(sTest)))
imErrMin = np.zeros((len(lTest), len(sTest)))
condNum = np.zeros((len(lTest), len(sTest)))

for i in range(len(lTest)):
    for j in range(len(sTest)):
        sigmaV, sigmaSketch, condN = testSingularValues(V, p = p, plot = False, l = lTest[i], s = sTest[j])
        imErr2[i,j] = norm(sigmaV - sigmaSketch)/norm(sigmaV)
        imErrMax[i,j] = abs( sigmaV[0] - sigmaSketch[0] )/sigmaV[0]
        imErrMin[i,j] = abs( sigmaV[-1] - sigmaSketch[-1] )/sigmaV[-1]
        condNum[i,j] = condN
    print(i)


# Plot as image
plt.figure()
plt.imshow(imErr2, cmap = colormap2, extent = [n, int(m/p), 1, len(lTest)], 
           origin="lower", aspect=0.3, vmin = 0, vmax = 2)
plt.title("L2 error singular values, V from colHadamard, r="+ str(r))
plt.xscale('symlog')
plt.yscale('symlog')
plt.colorbar()
plt.xlabel("l")
plt.ylabel("s")
plt.savefig("/Users/mmartine/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Mariana’s "+
            "MacBook Pro/EPFL/Block HRHT/Figures/Experiments epsilon subspace/p2imErr2_colHadamardZoom.png")

plt.figure()
plt.imshow(imErrMax, cmap = colormap2, extent = [n, int(m/p), 1, len(lTest)], 
           origin="lower", aspect=0.3, vmin = 0, vmax = 2)
plt.title("$|\sigma_{\max}(V) - \sigma_{\max}(HV)|$, V from colHadamard, r="+ str(r))
plt.xscale('symlog')
plt.yscale('symlog')
plt.colorbar()
plt.xlabel("l")
plt.ylabel("s")
plt.savefig("/Users/mmartine/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Mariana’s "+
            "MacBook Pro/EPFL/Block HRHT/Figures/Experiments epsilon subspace/p2imErrMax_colHadamardZoom.png")

plt.figure()
plt.imshow(imErrMin, cmap = colormap2, extent = [n, int(m/p), 1, len(lTest)], 
           origin="lower", aspect=0.3, vmin = 0, vmax = 2)
plt.title("$|\sigma_{\min}(V) - \sigma_{\min}(HV)|$, V from colHadamard, r="+ str(r))
plt.xscale('symlog')
plt.yscale('symlog')
plt.colorbar()
plt.xlabel("l")
plt.ylabel("s")
plt.savefig("/Users/mmartine/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Mariana’s "+
            "MacBook Pro/EPFL/Block HRHT/Figures/Experiments epsilon subspace/p2imErrMin_colHadamardZoom.png")


plt.figure()
plt.imshow(condNum, cmap = colormap2, extent = [n, int(m/p), 1, len(lTest)], 
           origin="lower", aspect=0.3)
plt.title("$\kappa(\Omega V)$, V from colHadamard, r="+ str(r))
plt.xscale('symlog')
plt.yscale('symlog')
plt.colorbar()
plt.xlabel("l")
plt.ylabel("s")
plt.savefig("/Users/mmartine/Library/Mobile Documents/com~apple~CloudDocs/Documents/Documents - Mariana’s "+
            "MacBook Pro/EPFL/Block HRHT/Figures/Experiments epsilon subspace/p2condNum_colHadamard.png")

# Pickle save
import pickle
file = open('data/p2imErr2_colHadamard.obj', 'wb')
pickle.dump(imErr2, file)
file.close()

file = open('data/p2imErrMax_colHadamard.obj', 'wb')
pickle.dump(imErrMax, file)
file.close()

file = open('data/p2imErrMin_colHadamard.obj', 'wb')
pickle.dump(imErrMin, file)
file.close()

file = open('data/p2condNum_colHadamard.obj', 'wb')
pickle.dump(condNum, file)
file.close()

# To read
# file = open("data/imErrMin_colHadamard.obj",'rb')
# imErrMinRead = pickle.load(file)
# file.close()

