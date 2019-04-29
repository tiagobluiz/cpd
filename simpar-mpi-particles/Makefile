CC = mpicc
CFLAGS = -std=c11 -lm
OBJ = simpar-mpi.o

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)

simpar-mpi: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
