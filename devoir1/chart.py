import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set_style()
plt.yticks(fontsize=18,)
plt.xticks(fontsize=18,)

df = pd.read_csv("time_lu_solve.csv", sep="\t")
df_fast = pd.read_csv("time_lu_solve_03.csv", sep="\t")

groups = df.groupby("function")
groups_fast = df_fast.groupby("function")

lu = groups.get_group("lu").groupby("size").mean()
solve = groups.get_group("solve").groupby("size").mean()
plt.title("Temps d'éxécution de la fonction lu et solve en fonction de la taille de la matrice", fontsize=25)
plt.loglog(range(1,1001),lu["time"], label="lu")
plt.loglog(range(1,1001),solve["time"], label="solve")
plt.xlabel("Taille Matrice NxN", fontsize=20)
plt.ylabel("Temps en seconde", fontsize=20)
plt.legend(fontsize=20)
plt.show()

lu_fast = groups_fast.get_group("lu").groupby("size").mean()
plt.title("Comparaison pour la méthode lu sans/avec les optimisations du compilateur", fontsize=25)
plt.loglog(range(1,1001),lu_fast["time"], label="Optimisations")
plt.loglog(range(1,1001), lu["time"], label="Sans optimisations")
plt.xlabel("Taille Matrice NxN", fontsize=20)
plt.ylabel("Temps en seconde", fontsize=20)
plt.legend(fontsize=20)
plt.show()

solve_fast = groups_fast.get_group("solve").groupby("size").mean()
plt.title("Comparaison pour la méthode solve sans/avec les optimisations du compilateur", fontsize=25)
plt.loglog(range(1,1001),solve_fast["time"], label="Optimisations")
plt.loglog(range(1,1001), solve["time"], label="Sans optimisations")
plt.xlabel("Taille Matrice NxN", fontsize=20)
plt.ylabel("Temps en seconde", fontsize=20)
plt.legend(fontsize=20)
plt.show()