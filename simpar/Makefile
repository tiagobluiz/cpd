CC = gcc
CFLAGS = -fopenmp -std=c11 -lm
OBJ = simpar.o

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)

simpar: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)