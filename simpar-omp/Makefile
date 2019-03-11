CC = gcc
CFLAGS = -fopenmp -std=c11 -lm
OBJ = simpar-omp.o

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)

simpar-omp: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)