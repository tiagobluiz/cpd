#include "Header.h"

particle_list * rmvlist( particle_list * list, particle_t * particle){
    particle_list * curr = list;
    while(curr->particle != particle && curr != NULL){
        curr = curr->next;
    }

    curr->next->prev = curr->prev;
    curr->prev->next = curr->next;

    return curr;
}

particle_list * addlist( particle_list * list, particle_t * particle){
    particle_list *newOne = (particle_list *)malloc(sizeof(particle_list));
    newOne->particle = particle;

    list->prev = newOne;
    newOne->next = list;

    return newOne;
}



cell ** createGrid(particle_t * particles, long long length, long ncSide){
    //Allocation of the cells
    cell ** grid = (cell **)malloc(ncSide*sizeof(cell*));
    for(int i = 0; i < ncSide; i++)
        grid[i] = (cell*)malloc(ncSide* sizeof(cell));

    double sizeCell = 1.0/ncSide;
    printf("CELLS %f \n", sizeCell);

    for(int i = 0; i < length; i++){
        long cellX = particles[i].x/sizeCell;
        long cellY = particles[i].y/sizeCell;
        cell * cell = &grid[cellX][cellY];
        particle_list *newOne = (particle_list *)malloc(sizeof(particle_list));
        if(grid[cellX][cellY].particles==NULL){
            newOne->particle = &particles[i];
            cell->particles = newOne;
        }else {
            cell->particles = addlist(cell->particles, &particles[i]);
        }

        particles[i].cellX = cellX;
        particles[i].cellY = cellY;

        printf("True CellX=%ld , True CellY=%ld , partilceX= %f , particleY= %f \n", cellX, cellY, particles[i].x, particles[i].y);
    }

    return grid;
}

//1 3 10 1
int main(int args_length, char* args[]) {
	if (args_length < 4) {
		printf("Incorrect number of arguments! It should be 4!");
		return -1;
	}

	//init
	long seed = strtol(args[1], NULL, 10);
    long ncside = strtol(args[2], NULL, 10);
	long long n_part = strtol(args[3], NULL, 10);
	long iterations = strtol(args[4], NULL, 10);

	particle_t * particles = malloc(n_part * sizeof(particle_t));

	init_particles(seed, ncside, n_part, particles);

	printf("Seed %ld, Ncside %ld, N_part %lld, Iterations %ld \n", seed, ncside, n_part, iterations);

	for(int i = 0; i < n_part; i++){
		printf("Index : %d \n", i);
		printf("X %f , Y %f, VX %f, VY %f, M %f \n",particles[i].x, particles[i].y, particles[i].vx, particles[i].vy, particles[i].m);
	}

    cell ** cellMatrix = createGrid(particles, n_part, ncside);

    for(int x =0;x<ncside;x++){
        for(int y =0;y<ncside;y++){
            printf("IndeX = %d IndeY = %d \n", x, y);
            particle_list * curr = cellMatrix[x][y].particles;
            while(curr!=NULL){
                printf("\tParticle \n \t  X = %lf Y= %lf \n", curr->particle->x, curr->particle->y);
                curr = curr->next;
            }
        }
    }

    free(cellMatrix);
	free(particles);
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

/**
 * Computes the center mass for each existing cell in grid
 *
 * @param cells array of cells that compose the grid
 * @param length number of particles (must match @cells length)
 */
void computeCellCenterMass(cell * cells, int length) {
    for (int cellIndex = 0; cellIndex < length; cellIndex++){
        particle_t * cellParticles = cells[cellIndex].particles;

        double totalY = 0;
        double totalX = 0;
        double totalMax = 0;


        for (int i = 0; i < length; i++) totalMax += cellParticles[i].m;                                                // Computes the total mass of the PIC


        for (int i = 0; i < length; i++) {                                                                              // Compute the influence of each PIC on center mass
            totalX += cells[i].m * cells[i].x;
            totalY += cells[i].m * cells[i].y;
        }


        //Affect out variable with center mass coordinates
        cells[cellIndex].x = totalX / totalMax;
        cells[cellIndex].y = totalY / totalMax;
        cells[cellIndex].m = totalMax;
    }
}


