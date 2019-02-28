#include "Header.h"
#define MINIMUM_DISTANCE 0,01

//struct that represents the force applied to a particle
typedef struct {
    double force;       // double which contains the force applied to the particle in N
    double angle;       //double which contains the angle of the force where the X axis is the base
} force_vector;


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!! if the particle and the center of mass are to close, the force should be 0
double computeMagnitudeForce(particle_t* a, cell * b) {
    double distance_a_b = sqrt ( exp2 ( a->x - b->x ) + exp2 ( a->y - b->y ));
    if(distance_a_b < MINIMUM_DISTANCE)
        return 0;
    return G * ( a->m * b->m ) / ( exp2 (distance_a_b));
}

/**
 * Function that gets the cell with the coordinates (unbound_row, unbound_column). Those coordinates will always represent
 * a cell because the space is wrapped, that is, if a particle exits through the side of the space it enters on the corresponding
 * position on the opposite side of the space.
 * @param unbound_row       indice in X axis
 * @param unbound_column    indice in Yaxis
 * @param cells             matrix with the cells to search
 * @param return_cell       cell which we were trying to get
 * @param side              number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void getCell(int unbound_row, int unbound_column, cell ** cells, cell * return_cell int side, double cell_dimension){
    int bounded_row = unbounded_row%side;
    int bounded_column = unbounded_column%side;

    if((unbounded_row >= side && unbounded_column >= side) || unbounded_row < 0 && unbounded_column < 0){
        return_cell->x = unbounded_column*cell_dimension + cells[bounded_row][bounded_column].x;
        return_cell->y = unbounded_row*cell_dimension + cells[bounded_row][bounded_column].y;
        return_cell->m = cells[bounded_row][bounded_column].m;
    }
    else if(unbounded_row >= side || unbounded_row < 0){
        return_cell->y = unbounded_row*cell_dimension + cells[bounded_row][bounded_column].y;
        return_cell->x = cells[bounded_row][bounded_column].x;
        return_cell->m = cells[bounded_row][bounded_column].m;
    }
    else if(unbounded_column >= side || unbounded_column < 0){
        return_cell->x = unbounded_column*cell_dimension + cells[bounded_row][bounded_column].x;
        return_cell->y = cells[bounded_row][bounded_column].y;
        return_cell->m = cells[bounded_row][bounded_column].m;
    }
    else {
        return_cell = &cells[bounded_row][bounded_column];
    }
}

// we can create DUMMY CELLS around the grid!!
/**
 * Function that calculates the force being applied to all particles
 * @param particles         list with the particles
 * @param particlesLength   number of particles in the list
 * @param cells             matrix that represents the cells
 * @param forceVector       list to return the resultant force being applied to each particle
 * @param side              number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void calculateForce(particle_t * particles, int particlesLength, cell ** cells, force_vector * forceVector, int side, double cell_dimension){

    //iterate all cells
    for(int i = 0; i < particlesLength; i++ ){

        //change to get the coordinates of the cell where the particle is located !!!!!!!!!!!!!!!!!!!1
        double cell_x = particles[i].x;
        double cell_y = particles[i].y;

        //resultant force in X and Y
        double Fx = 0;
        double Fy = 0;

        // iterate all the 9 mass cells around the particle
        for(int row = -1; row < 2; row++){
            for(int column = -1; column < 2; column++){
                //get the neighbor cell in the coordinates (row, column)
                cell neighbor_cell;
                getCell(row + cell_x, column + cell_y, cells, &neighbor_cell, side, cell_dimension);

                //compute angle
                double delta_x = neighbor_cell.x - particles[i].x;
                double delta_y = neighbor_cell.y - particles[i].y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force !!!!!!!!!!!!!!!! verificar se é assim
                double force = computeMagnitudeForce(&particles[i], &neighbor_cell);

                //compute resultant force
                // it's in radians or degrees????????????
                Fx += force*cos(vector_angle);
                Fy += force*((M_PI/2) - vector_angle);
            }
        }

        //calculate the resultant force and the respective angle
        forceVector[i].force = sqrt( Fx*Fx + Fy*Fy );
        forceVector[i].angle = acos(Fx/resultant_force);
    }
}

