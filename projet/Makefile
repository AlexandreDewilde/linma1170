
# Localisation du SDK Gmsh ; pour l'évaluation il sera bien dans le dossier parent
GMSH_DIR := ../gmsh-sdk

# Choix du compilateur
CC := gcc

# Flags de compilation ; vous pouvez ajouter d'autres flags si vous le souhaitez
CFLAGS := -Wall

# Chemins vers les dossiers `include` où se trouvent les headers externes
# on donne le chemin vers le dossier `include` du SDK
INC_DIR := -I $(GMSH_DIR)/include # specify include directories

# Chemins vers les dossiers `lib` où se trouvent les librairies externes
# on donne le chemin vers le dossier `lib` du SDK
LIB_DIR := -L $(GMSH_DIR)/lib # specify library directories for compilation

# Spécification du runtime path
LDFLAGS := -Wl,-rpath,$(GMSH_DIR)/lib

# On indique au compilateur de linker la librairie Gmsh
LDLIBS := -lgmsh -llapack -lm -lcblas -llapacke

# Nom du programme à compiler
PROG := project

# Liste des objets nécessaires pour compiler le programme
# à modifier si vous ajoutez d'autres modules !
OBJS := elasticity.o lu.o matrix.o design.o eigen.o rmck.o

# Règle de compilation
all: $(PROG)

python_functions: $(OBJS) python_functions.o
	$(CC) -fPIC $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3 -shared -o python_functions.so

gradient_descent: $(OBJS) gradient_descent.o
	$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3

bisect: $(OBJS) bisect.o
	$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3

constant_area_bisect: $(OBJS) constant_area_bisect.o
	$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3

ternary_search: $(OBJS) ternary_search.o
	$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3

multi_layer: $(OBJS) multi_layer.o
	$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3


# Règle de compilation : on produit un .o à partir d'un .c
# on fournit les flags de compilation CFLAGS et les chemins include INC_DIR
# `$<` signifie: "la première dépendance" (le .c)
%.o: %.c
	$(CC) -g -c $(CFLAGS) $(INC_DIR) $< -O3

# Règle de link : on produit le programme PROG à partir des fichiers objets OBJS
# on fournit la librairie à linker et le runtime path
# `$^` signifie : "toutes les dépendances" (les OBJS)
$(PROG): $(OBJS) project.o
	$(CC) -g -o $@ $(LIB_DIR) $(LDLIBS) $(LDFLAGS) $^ -O3

# Règle de nettoyage : supprime PROG et OBJS
clean:
	rm *.o -f
	rm $(PROG) $(OBJS) -f
	rm *.so -f
	rm -f project.o project bisect execute ternary_search execute_2_layer