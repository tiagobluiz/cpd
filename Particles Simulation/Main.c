#include "Header.h"

/**
 * Update the matrix coordenates of each particle
 *
 * @param particle A pointer to a particle
 * @param ncSide The number of cells on each side of the matrix of cells
 */
void updateParticle(particle_t * particle, long ncSide){

    double sizeCell = 1.0/ncSide;

    particle->cellX = particle->x/sizeCell;
    particle->cellY = particle->y/sizeCell;
}
/**
 *
 * @param particles Array that contains all the particles
 * @param length The size of the particles array
 * @param ncSide The number of cells on each side of the matrix of cells
 * @return Matrix of cells and update the matrix coordenates of each particle
 */
cell ** createGrid(particle_t * particles, long long length, long ncSide, double * cell_dimension){
    //Allocation of the cells
    cell ** grid = (cell **) malloc(ncSide * ncSide * sizeof(cell*));
    for(int i = 0; i < ncSide; i++)
        grid[i] = (cell*)calloc(ncSide, sizeof(cell)); // review this alloc - TODO

    *cell_dimension = 1.0/ncSide;

    for(int i = 0; i < length; i++){
        updateParticle(&particles[i], ncSide);
    }

    return grid;
}

/**
 * This function finds the remainder after division of one number by another
 *
 * @param dividend  the number that will be divided
 * @param diviser   the number that divides
 * @return          the remainder (always positive) of the operation
 */
int mod (int dividend, int diviser){
    int result = dividend % diviser;
    return result < 0? result + diviser : result;
}

double computeAcceleration(double force, double mass) {
	return force / mass;
}

double computeVelocity(double currVelocity, double acceleration){
    return currVelocity * acceleration;
}

void computeParticlePosition(double * x, double * y, double velocity, double acceleration){
    *x = *x + velocity + (1/2) * acceleration;
    *y = *y + velocity + (1/2) * acceleration;
}

/**
 * Computes the center mass for each existing cell in grid
 *
 * @param particles array of particles that compose the grid
 * @param length number of particles (must match @particles length)
 * @param cells bidimensional array with the grid of cells
 * @param ncside sides of the grid (how many rows the grid has)
 */
void computeCellCenterMass(particle_t * particles, long length, cell ** cells, long ncside) {
    for (long particleIndex = 0; particleIndex < length; particleIndex++){
        particle_t particle = particles[particleIndex];
        cells[particle.cellX][particle.cellY].x += particle.m * particle.x;
        cells[particle.cellX][particle.cellY].y += particle.m * particle.y;
        cells[particle.cellX][particle.cellY].m += particle.m;
        cells[particle.cellX][particle.cellY].part += 1;
    }

    for (long cellRowIndex = 0; cellRowIndex < ncside; cellRowIndex++){
        for (long cellColumnIndex = 0; cellColumnIndex < ncside; cellColumnIndex++){
            cells[cellRowIndex][cellColumnIndex].x /= cells[cellRowIndex][cellColumnIndex].m;
            cells[cellRowIndex][cellColumnIndex].y /= cells[cellRowIndex][cellColumnIndex].m;
        }
    }
}

/**
 * Compute the magnitude force between a particle and a center mass
 * If the particle and the center of mass are to close, the force should be 0
 *
 * @param a  the particle
 * @param b  the cell
 * @return   a double that represents the force being applied to the particle
 */
double computeMagnitudeForce(particle_t* a, cell * b) {
    double distance_a_b = sqrt ( exp2 ( a->x - b->x ) + exp2 ( a->y - b->y ));
    if(distance_a_b < MINIMUM_DISTANCE)
        return 0;
    return G * ( a->m * b->m ) / ( exp2 (distance_a_b));
}

/**
 * Method that updates the new position of the particle based on the force being applied
 *
 * @param particle  particle being updated
 * @param Fx        Force in axis X
 * @param Fy        Force in axis Y
 */
void updateParticlePosition(particle_t * particle, double Fx, double Fy, long ncside){
    double acceleration_x = Fx/particle->m;
    double acceleration_y = Fy/particle->m;
    particle->vx += acceleration_x;
    particle->vy += acceleration_y;
    double x= particle->x + particle->vx + acceleration_x/2;
    double y= particle->y + particle->vy + acceleration_y/2;
    x = fmod(x, MAX_COORDINATES_VALUE);
    y = fmod(y, MAX_COORDINATES_VALUE);
    particle->x = x < 0 ? x + MAX_COORDINATES_VALUE : x;
    particle->y = y < 0 ? y + MAX_COORDINATES_VALUE : y;

    updateParticle(particle, ncside);
}

/**
 * Function that gets the cell with the coordinates (unbound_row, unbound_column). Those coordinates will always represent
 * a cell because the space is wrapped, that is, if a particle exits through the side of the space it enters on the corresponding
 * position on the opposite side of the space.
 *
 * @param unbound_row       indice in X axis
 * @param unbound_column    indice in Y axis
 * @param cells             matrix with the cells to search
 * @param return_cell       cell which we were trying to get
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
cell * getCell(long long unbounded_row, long long unbounded_column, cell ** cells, cell * return_cell, long ncside, double cell_dimension){
    int bounded_row = mod(unbounded_row, ncside);
    int bounded_column = mod(unbounded_column, ncside);

    if((unbounded_row >= ncside && unbounded_column >= ncside) || (unbounded_row < 0 && unbounded_column < 0)){
        return_cell->x = unbounded_column*cell_dimension + cells[bounded_row][bounded_column].x;
        return_cell->y = unbounded_row*cell_dimension + cells[bounded_row][bounded_column].y;
        return_cell->m = cells[bounded_row][bounded_column].m;
    }
    else if(unbounded_row >= ncside || unbounded_row < 0){
        return_cell->y = unbounded_row*cell_dimension + cells[bounded_row][bounded_column].y;
        return_cell->x = cells[bounded_row][bounded_column].x;
        return_cell->m = cells[bounded_row][bounded_column].m;
    }
    else if(unbounded_column >= ncside || unbounded_column < 0){
        return_cell->x = unbounded_column*cell_dimension + cells[bounded_row][bounded_column].x;
        return_cell->y = cells[bounded_row][bounded_column].y;
        return_cell->m = cells[bounded_row][bounded_column].m;
    }
    else {
        return &cells[bounded_row][bounded_column];
    }

    return return_cell;
}

/**
 * Function that computes the force being applied to all particles and updates their cell after being moved
 *
 * @param particles         list with the particles
 * @param particlesLength   number of particles in the list
 * @param cells             matrix that represents the cells
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void computeForceAndUpdateParticles(particle_t *particles, int particlesLength, cell **cells, long ncside,
                                    double cell_dimension){

    //iterate all particles
    for(int i = 0; i < particlesLength; i++ ){

        //get the coordinates of the cell where the particle is located
        long long cell_x = particles[i].cellX;
        long long cell_y = particles[i].cellY;

        //resultant force in X and Y
        double Fx = 0;
        double Fy = 0;

        // iterate all the 9 mass cells around the particle
        for(int row = -1; row < 2; row++){
            for(int column = -1; column < 2; column++){
                //get the neighbor cell in the coordinates (row, column)
                cell created_neighbor_cell;
                cell * neighbor_cell = getCell(row + cell_x, column + cell_y, cells, &created_neighbor_cell, ncside, cell_dimension);

                //if that cell doesn't have any particle, as a consequence no central mass will exist, so we ignore
                if(neighbor_cell->m == 0)
                    continue;

                //compute angle
                double delta_x = neighbor_cell->x - particles[i].x;
                double delta_y = neighbor_cell->y - particles[i].y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force
                double force = computeMagnitudeForce(&particles[i], neighbor_cell);

                Fx += force * cos(vector_angle);
                Fy += force * sin(vector_angle);
            }
        }

        //update the particle position
        updateParticlePosition(&particles[i], Fx, Fy, ncside);
    }
}

void computeOverallCenterOfMass(cell ** cells, long ncside){
    cell overallCenterMass;// = {x:0, y:0, m:0};
    for (long cellRowIndex = 0; cellRowIndex < ncside; cellRowIndex++){
        for (long cellColumnIndex = 0; cellColumnIndex < ncside; cellColumnIndex++){
            if(cells[cellRowIndex][cellColumnIndex].m == 0) continue;
            overallCenterMass.x += cells[cellRowIndex][cellColumnIndex].x * cells[cellRowIndex][cellColumnIndex].m;
            overallCenterMass.y += cells[cellRowIndex][cellColumnIndex].y * cells[cellRowIndex][cellColumnIndex].m;
            overallCenterMass.m += cells[cellRowIndex][cellColumnIndex].m;
        }
    }
    overallCenterMass.x /= overallCenterMass.m;
    overallCenterMass.y /= overallCenterMass.m;
    printf("%0.2f  %0.2f", overallCenterMass.x, overallCenterMass.y);

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

    particle_t * particles = calloc(n_part, sizeof(particle_t));

    init_particles(seed, ncside, n_part, particles);

    printf("Seed %ld, Ncside %ld, N_part %lld, Iterations %ld \n", seed, ncside, n_part, iterations);

    double cell_dimension = 0;
    cell ** cellMatrix = createGrid(particles, n_part, ncside, &cell_dimension);

    for(int i = 0; i < iterations; i++){
        computeCellCenterMass(particles, n_part, cellMatrix, ncside);
        computeForceAndUpdateParticles(particles, n_part, cellMatrix, ncside, cell_dimension);
    }

    printf("%0.2f %0.2f \n", particles[0].x, particles[0].y);
    computeOverallCenterOfMass(cellMatrix, ncside);
    free(cellMatrix);
    free(particles);
}







