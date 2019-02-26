#include "Header.h"

void main(char * args []) {
	long seed = strtol(args[0], NULL, 10);
    long ncside = strtol(args[1], NULL, 10);
	long long n_part = strtol(args[2], NULL, 10);
	long iterations = strtol(args[3], NULL, 10);

	const long number_of_particles = n_part;
	particle_t * particles = malloc(n_part * sizeof(particle_t));

	

	init_particles(seed, ncside, n_part, particles);
}