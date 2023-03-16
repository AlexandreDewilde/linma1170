CC=gcc
FLAGS=-lm -lgmsh -Wall
FILES=lu.c datastructures/matrix.c datastructures/vector.c elasticity.c rmck.c


solve_deformation: $(FILES) solve_deformation.c
	$(CC) $(FILES) solve_deformation.c $(FLAGS) -o solve_deformation -g -O3

clean:
	rm solve_deformation 