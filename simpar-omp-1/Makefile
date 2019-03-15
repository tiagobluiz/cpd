CC = gcc-8
CFLAGS = -fopenmp -std=c11 -lm
OBJ = simpar-omp.1.o

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)

simpar-omp.1: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)