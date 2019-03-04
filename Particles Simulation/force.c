#include "Header.h"
#define MINIMUM_DISTANCE 0.01

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
    particle->x += particle->vx + acceleration_x/2;
    particle->y += particle->vy + acceleration_y/2;


    updateParticle(particle, ncside);

}

/**
 * Function that gets the cell with the coordinates (unbound_row, unbound_column). Those coordinates will always represent
 * a cell because the space is wrapped, that is, if a particle exits through the side of the space it enters on the corresponding
 * position on the opposite side of the space.
 *
 * @param unbound_row       indice in X axis
 * @param unbound_column    indice in Yaxis
 * @param cells             matrix with the cells to search
 * @param return_cell       cell which we were trying to get
 * @param ncside              number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
cell * getCell(long long unbounded_row, long long unbounded_column, cell ** cells, cell * return_cell, long ncside, double cell_dimension){
    long long bounded_row = unbounded_row%ncside;
    long long bounded_column = unbounded_column%ncside;

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
 * Function that calculates the force being applied to all particles
 *
 * @param particles         list with the particles
 * @param particlesLength   number of particles in the list
 * @param cells             matrix that represents the cells
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void calculateForce(particle_t * particles, int particlesLength, cell ** cells, long ncside, double cell_dimension){

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

                //compute angle
                double delta_x = neighbor_cell->x - particles[i].x;
                double delta_y = neighbor_cell->y - particles[i].y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force
                double force = computeMagnitudeForce(&particles[i], neighbor_cell);

                Fx += force*cos(vector_angle);
                Fy += force*sin(vector_angle);
            }
        }

        //update the particle position
        updateParticlePosition(&particles[i], Fx, Fy, ncside);
    }
}

