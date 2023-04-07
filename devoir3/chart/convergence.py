import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("convergence.csv", sep=";")

gp = df.groupby("meshsize").mean()

ax = plt.subplot()
ax.bar(df["meshsize"].unique()-0.01, gp["power_iteration"], align="center", width=0.02, label="Power iteration")
ax.bar(df["meshsize"].unique()+0.01, gp["rayleigh_quotient_iteration"], align="center", width=0.02, label="rayleigh_quotient_iteration")
ax.set_xlabel("Meshsize", fontsize=18)
ax.set_ylabel("Nombres d'itérations moyenne", fontsize=18)
plt.legend(fontsize=18)
plt.show()