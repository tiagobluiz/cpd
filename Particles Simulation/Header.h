#pragma once

#include <stdio.h>
#include <stdlib.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

typedef struct {
	double x;
	double y;
	double vx;
	double vy;
	double m;
}particle_t;

void init_particles(long seed, long ncside, long long n_part, particle_t* par);