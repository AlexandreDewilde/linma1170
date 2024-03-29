
# Localisation du SDK Gmsh ; pour l'évaluation il sera bien dans le dossier parent
GMSH_DIR := ../gmsh-sdk

# Choix du compilateur
CC := gcc

# Flags de compilation ; vous pouvez ajouter d'autres flags si vous le souhaitez
CFLAGS := -Wall -O3

# Chemins vers les dossiers `include` où se trouvent les headers externes
# on donne le chemin vers le dossier `include` du SDK
INC_DIR := -I $(GMSH_DIR)/include # specify include directories

# Chemins vers les dossiers `lib` où se trouvent les librairies externes
# on donne le chemin vers le dossier `lib` du SDK
LIB_DIR := -L $(GMSH_DIR)/lib # specify library directories for compilation

# Spécification du runtime path
LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib

# On indique au compilateur de linker la librairie Gmsh
LDLIBS := -lgmsh -lm

# Nom du programme à compiler
PROG := solve_eigenmodes

SRC := $(wildcard *.c)
# Liste des objets nécessaires pour compiler le programme
# à modifier si vous ajoutez d'autres modules !
# OBJS := solve_eigenmodes.o elasticity.o lu.o matrix.o vector.o rmck.o power_method.o rayleigh_quotient.o reduce_matrix.o
SRCS := functions/elasticity.c functions/lu.c functions/matrix.c functions/vector.c functions/rmck.c functions/power_method.c functions/rayleigh_quotient.c

# Règle de compilation
all: $(PROG)

# Règle de compilation : on produit un .o à partir d'un .c
# on fournit les flags de compilation CFLAGS et les chemins include INC_DIR
# `$<` signifie: "la première dépendance" (le .c)
%.o: %.c
	$(CC) -c $(CFLAGS) $(INC_DIR) $<

# Règle de link : on produit le programme PROG à partir des fichiers objets OBJS
# on fournit la librairie à linker et le runtime path
# `$^` signifie : "toutes les dépendances" (les OBJS)
$(PROG): solve_eigenmodes.c $(SRCS)
	@$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^

test_convergence: test_convergence.c $(SRCS)
	@$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^

test_hybrid: test_hybrid.c $(SRCS)
	@$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^
	
test_inversion: test_inversion.c $(SRCS)
	@$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^

# Règle de nettoyage : supprime PROG et OBJS
clean:
	rm -f solve_eigenmodes
	rm -f test_convergence
	rm -f test_inversion
	rm -f test_hybrid