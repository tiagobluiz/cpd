#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.01
#define MAX_COORDINATES_VALUE 1.0

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
    int part;
}cell;

void init_particles(long seed, long ncside, long long n_part, particle_t *par)
{
    long long i;

    srandom(seed);

    for(i = 0; i < n_part; i++)
    {
        par[i].x = RND0_1;
        par[i].y = RND0_1;
        par[i].vx = RND0_1 / ncside / 10.0;
        par[i].vy = RND0_1 / ncside / 10.0;

        par[i].m = RND0_1 * ncside / (G * 1e6 * n_part);
    }
}

//#include <omp.h>

/**
 * Utility Methods
 */

/**
 * Update the matrix coordinates of each particle
 *
 * @param particle  A pointer to a particle
 * @param ncside    The number of cells on each side of the matrix of cells
 */
void update_particle(particle_t *particle, long ncside){
    double sizeCell = MAX_COORDINATES_VALUE/ncside;
    particle->cellX = particle->x/sizeCell;
    particle->cellY = particle->y/sizeCell;
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
    return result < 0 ? result + diviser : result;
}

/**
 * Cleans all cells and resets all parameters to 0 so that the previous state doesn't affect the new one
 * @param cells     bidimensional array with the grid of cells
 * @param ncside    sides of the grid (how many rows the grid has)
 */
void clean_cells(cell ** cells, long ncside){
    for (long cellRowIndex = 0; cellRowIndex < ncside; cellRowIndex++){
        for (long cellColumnIndex = 0; cellColumnIndex < ncside; cellColumnIndex++){
            cells[cellRowIndex][cellColumnIndex].x = cells[cellRowIndex][cellColumnIndex].y =
            cells[cellRowIndex][cellColumnIndex].m = cells[cellRowIndex][cellColumnIndex].part = 0;
        }
    }
}

/**
 * Implementation
 */


/**
 *
 * @param particles Array that contains all the particles
 * @param length    The size of the particles array
 * @param ncside    The number of cells on each side of the matrix of cells
 * @return Matrix of cells and update the matrix coordenates of each particle
 */
cell ** create_grid(particle_t * particles, long long length, long ncside, double *cell_dimension){
    //Allocation of the cells
    cell ** grid = (cell **) malloc(sizeof(cell*) * ncside + sizeof(cell) * ncside * ncside);

    cell * ptr = grid + ncside;
    for(int i = 0; i < ncside; i++)
        grid[i] = (ptr + ncside * i);

    *cell_dimension = MAX_COORDINATES_VALUE/ncside;

    for(int i = 0; i < length; i++){
        update_particle(&particles[i], ncside);
    }

    return grid;
}

/**
 * Computes the center mass for each existing cell in grid
 *
 * @param particles array of particles that compose the grid
 * @param length    number of particles (must match @particles length)
 * @param cells     bidimensional array with the grid of cells
 * @param ncside    sides of the grid (how many rows the grid has)
 */
void compute_cell_center_mass(particle_t *particles, long length, cell ** cells, long ncside) {

  /*  #pragma omp parallel{
    /**
     * criar um array por thread, e cada um ter o seu. O array teria o length de ncside + ncside
     * dado que cada thread teria o seu array, a posição da celula a ser afetada seria dada pela soma do cellX + cellY
     * no final um outro for juntaria essas threads
     *//*
//        cell ** temporaryCells[#NUM_THREADS][ncside + ncside];
        #pragma omp for{
            for (long particleIndex = 0; particleIndex < length; particleIndex++) {
                particle_t particle = particles[particleIndex];
                cells[omp_get_thread_num()][particle.cellX + particle.cellY].x += particle.m * particle.x;
                cells[omp_get_thread_num()][particle.cellX + particle.cellY].y += particle.m * particle.y;
                cells[omp_get_thread_num()][particle.cellX + particle.cellY].m += particle.m;
                cells[omp_get_thread_num()][particle.cellX + particle.cellY].part += 1;
            }
        }
        #pragma omp for{
            for (int threadUsed = 0; threadUsed < #NUM_THREADS; threadUsed++){
                //ir buscar a cada array e juntar no final a cells
            }
        }
    };*/

    for (long particleIndex = 0; particleIndex < length; particleIndex++){
        particle_t particle = particles[particleIndex];
        cells[particle.cellX][particle.cellY].x += particle.m * particle.x;
        cells[particle.cellX][particle.cellY].y += particle.m * particle.y;
        cells[particle.cellX][particle.cellY].m += particle.m;
        cells[particle.cellX][particle.cellY].part += 1;
    }

    //#pragma omp parallel for
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
double compute_magnitude_force(particle_t * a, cell * b) {
    double distance_a_b = sqrt ( exp2 ( a->x - b->x ) + exp2 ( a->y - b->y ));
    if(distance_a_b < EPSLON)
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
void update_particle_position(particle_t * particle, double Fx, double Fy, long ncside){
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

    update_particle(particle, ncside);
}

/**
 * Function that gets the cell with the coordinates (unbound_row, unbound_column). Those coordinates will always represent
 * a cell because the space is wrapped, that is, if a particle exits through one side of the space, it enters on the corresponding
 * opposite side of the space.
 *
 * @param unbound_row       indice unbounded in X axis
 * @param unbound_column    indice unbounded in Y axis
 * @param cells             matrix with the cells to search
 * @param return_cell       cell used to represent the unbounded cell
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 * @return                  the neighbor cell
 */
cell * get_cell(long long unbounded_row, long long unbounded_column, cell **cells, cell * return_cell, long ncside,
                double cell_dimension){
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
 * @param particles_length  number of particles in the list
 * @param cells             matrix that represents the cells
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void compute_force_and_update_particles(particle_t *particles, int particles_length, cell **cells, long ncside,
                                        double cell_dimension){

    //iterate all particles
    for(int i = 0; i < particles_length; i++ ){

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
                cell * neighbor_cell = get_cell(row + cell_x, column + cell_y, cells, &created_neighbor_cell, ncside,
                                                cell_dimension);

                //if that cell doesn't have any particle, no central mass will exist, so we ignore
                if(neighbor_cell->m == 0)
                    continue;

                //compute angle
                double delta_x = neighbor_cell->x - particles[i].x;
                double delta_y = neighbor_cell->y - particles[i].y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force
                double force = compute_magnitude_force(&particles[i], neighbor_cell);
                Fx += force * cos(vector_angle);
                Fy += force * sin(vector_angle);
            }
        }

        //update the particle position
        update_particle_position(&particles[i], Fx, Fy, ncside);
    }
}

/**
 * Computes the center mass for each existing cell in grid
 *
 * @param particles     array of particles that compose the grid
 * @param length        number of particles (must match @particles length)
 */
void compute_overall_center_mass(particle_t * particles, long length){
    cell overallCenterMass = {x:0, y:0, m:0};
    for (long particleIndex = 0; particleIndex < length; particleIndex++){
        overallCenterMass.x += particles[particleIndex].x * particles[particleIndex].m;
        overallCenterMass.y += particles[particleIndex].y * particles[particleIndex].m;
        overallCenterMass.m += particles[particleIndex].m;
    }
    overallCenterMass.x /= overallCenterMass.m;
    overallCenterMass.y /= overallCenterMass.m;
    printf("%0.2f %0.2f \n", overallCenterMass.x, overallCenterMass.y);

}

int main(int args_length, char* args[]) {
    if (args_length < 4) {
        printf("Incorrect number of arguments! It should be 4!");
        return -1;
    }

    long seed = strtol(args[1], NULL, 10);
    long ncside = strtol(args[2], NULL, 10);
    long long n_part = strtol(args[3], NULL, 10);
    long iterations = strtol(args[4], NULL, 10);

    particle_t * particles = calloc(n_part, sizeof(particle_t));

    init_particles(seed, ncside, n_part, particles);

    double cell_dimension = 0;
    cell ** cellMatrix = create_grid(particles, n_part, ncside, &cell_dimension);

    for(int i = 0; i < iterations; i++){
        compute_cell_center_mass(particles, n_part, cellMatrix, ncside);
        compute_force_and_update_particles(particles, n_part, cellMatrix, ncside, cell_dimension);
        clean_cells(cellMatrix, ncside);
    }

    printf("%0.2f %0.2f \n", particles[0].x, particles[0].y);
    compute_overall_center_mass(particles, n_part);
    free(cellMatrix);
    free(particles);
}







