import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("hybrid.csv", sep=";")

gp = df.groupby("meshsize").mean()

ax = plt.subplot()
ax.bar(df["meshsize"].unique()-0.01, gp["power_iteration"], align="center", width=0.02, label="Power iteration")
ax.bar(df["meshsize"].unique()+0.01, gp["hybrid_pi"], align="center", width=0.02, label="Hybrid : Power iteration")
ax.bar(df["meshsize"].unique()+0.01, gp["hybrid_rqi"], bottom=gp["hybrid_pi"], align="center", width=0.02, label="Hybrid : Rayleigh quotient iteration")
ax.set_xlabel("Meshsize", fontsize=18)
ax.set_ylabel("Nombres d'itérations", fontsize=18)
plt.legend(fontsize=18)
plt.show()