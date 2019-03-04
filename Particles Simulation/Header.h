#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01

typedef struct {
	double x;
	double y;
	double vx;
	double vy;
	long long cellX;
	long long cellY;
	double m;
}particle_t;

typedef struct {
	double x;
	double y;
	double m;
}cell;

void init_particles(long seed, long ncside, long long n_part, particle_t* par);