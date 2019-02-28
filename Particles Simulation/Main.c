#include "Header.h"

int main(int args_length, char* args[]) {
	if (args_length < 4) {
		printf("Incorrect number of arguments! It should be 4!");
		return -1;
	}

	//init
	long seed = strtol(args[0], NULL, 10);
    long ncside = strtol(args[1], NULL, 10);
	long long n_part = strtol(args[2], NULL, 10);
	long iterations = strtol(args[3], NULL, 10);

	particle_t * particles = malloc(n_part * sizeof(particle_t));

	init_particles(seed, ncside, n_part, particles);




	free(particles);
}

double computeMagnitudeForce(particle_t* a, particle_t* b) {
	return G * ( a->m * b->m ) / ( exp2 ( sqrt ( exp2 ( a->x - b->x ) + exp2 ( a->y - b->y ))));
}

double computeAcceleration(double force, double mass) {
	return force / mass;
}

double computeVelocity(double currVelocity, double acceleration, double delta){
    return currVelocity * acceleration * delta;
}

void computeParticlePosition(double * x, double * y, double velocity, double acceleration, double delta){
    *x = *x + velocity + 1/2 * acceleration * exp2(delta);
    *y = *y + velocity + 1/2 * acceleration * exp2(delta);
}

void computeCellMass(particle_t * particles, int length, double * X, double * Y) {

    for (int i = 0; i < length; i++){
        
    }


	double totalY = 0;
	double totalX = 0;
	double totalMax = 0;
	
	for (int i = 0; i < length; i++) {
		totalMax += particles[i].m;
	}


	for (int i = 0; i < length; i++) {
		totalX += particles[i].m * particles[i].x;
		totalY += particles[i].m * particles[i].y;
	}


	//Affect out variable with center mass coordinates
	*X = totalX / totalMax;
	*Y = totalY / totalMax;
}

