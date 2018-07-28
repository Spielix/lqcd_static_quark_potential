import math
import numpy as np
from scipy.stats import norm
from scipy.stats import describe
import sys

A = float(sys.argv[1])
dA = float(sys.argv[2])
B = float(sys.argv[3])
dB = float(sys.argv[4])
corAB = float(sys.argv[5])
C = float(sys.argv[6])
dC = float(sys.argv[7])
D = float(sys.argv[8])
dD = float(sys.argv[9])
corCD = float(sys.argv[10])
N = int(sys.argv[11])

phiAB = 0.5 * math.asin(corAB)
phiCD = 0.5 * math.asin(corCD)

mAB = np.array([[math.sin(phiAB), math.cos(phiAB)], [math.cos(phiAB), math.sin(phiAB)]])
mCD = np.array([[math.sin(phiCD), math.cos(phiCD)], [math.cos(phiCD), math.sin(phiCD)]])

x = np.zeros(N);
for i in range(N):
    ABsamp = np.array([norm.rvs(size=1, scale=dA), norm.rvs(size=1, scale=dB)])
    CDsamp = np.array([norm.rvs(size=1, scale=dC), norm.rvs(size=1, scale=dD)])
    ABdat_corr = np.array([[A], [B]]) + np.dot(mAB, ABsamp)
    CDdat_corr = np.array([[C], [D]]) + np.dot(mCD, CDsamp)
#    print(mAB)
#    print(ABsamp)
#    print(ABdat_corr)
    x[i] =  (CDdat_corr[1] - ABdat_corr[1]) / (ABdat_corr[0] - CDdat_corr[0]) 

desc = describe(x)
print((D-B)/(A-C))
print(desc.mean)
print(math.sqrt(desc.variance))
