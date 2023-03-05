CC=gcc

main: lu.c matrix.c main.c
	@$(CC) lu.c matrix.c main.c -o main -O3 -lm
	@./main

no_optimization: lu.c matrix.c main.c
	@$(CC) lu.c matrix.c main.c -o main -lm
	@./main


memory_test: lu.c matrix.c main.c
	$(CC) lu.c matrix.c main.c -o main -lm
	@valgrind ./main

clean:
	@rm main -f
	@rm vgcore.* -f