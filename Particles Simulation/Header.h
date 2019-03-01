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

typedef struct particle_list {
	struct particle_list * prev;
	struct particle_list * next;
	particle_t * particle;
}particle_list;

typedef struct {
	double x;
	double y;
	double m;
	particle_t * particles;
}cell;

void init_particles(long seed, long ncside, long long n_part, particle_t* par);