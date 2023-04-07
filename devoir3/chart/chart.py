import pandas as pd
import matplotlib.pyplot as plt

perm_matrix = pd.read_csv("tuning_fork_perm.csv")
non_perm_matrix = pd.read_csv("tuning_fork_non_perm.csv")

def s_to_time(x):
    m, s = x[:-1].split("m")
    return int(m)*60 + float(s)
perm_matrix["time"] = perm_matrix["time"].apply(s_to_time)
non_perm_matrix["time"] = non_perm_matrix["time"].apply(s_to_time)
time_perm = perm_matrix.groupby("meshSize").mean()

time_non_perm = non_perm_matrix.groupby("meshSize").mean()

plt.plot(perm_matrix["meshSize"].unique(), time_perm, label="matrix permuté")
plt.plot(non_perm_matrix["meshSize"].unique(), time_non_perm, label="matrix non permuté")
plt.xlabel("Meshsize", fontsize=18)
plt.ylabel("Temps (s)", fontsize=18)
plt.legend(fontsize=18)
plt.show()