import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style()

K = np.loadtxt("K_perm.csv", skiprows=12)
lu = np.loadtxt("LU_perm.csv", skiprows=12)
# print(sc.linalg.bandwidth(K), sc.linalg.bandwidth(lu))
plt.spy(K)
plt.show()
plt.spy(lu)
plt.show()