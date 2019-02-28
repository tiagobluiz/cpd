#include "Header.h"

typedef struct {
    double force;
    double angle;
} force_vector;

double computeMagnitudeForce(particle_t* a, cell * b) {
    return G * ( a->m * b->m ) / ( exp2 ( sqrt ( exp2 ( a->x - b->x ) + exp2 ( a->y - b->y ))));
}


// we can create DUMMY CELLS around the grid!!
void calculateForce(particle_t * particles, int particlesLength, cell ** cells, int cellsLength, force_vector * forceVector, int side, double dimension_cell){

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
                i_row = row + cell_x;
                i_column = column + cell_y;

                if(i_row > side){

                }
                //get respective cell
                cell * p_cell = &cells[()%side][()%side];

                //compute angle
                double delta_x = p_cell.x - particles[i].x;
                double delta_y = p_cell.y - particles[i].y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force !!!!!!!!!!!!!!!! verificar se é assim
                double force = computeMagnitudeForce(&particles[i], p_cell);

                //compute resultant force
                // it's in radians or degrees????????????
                Fx += force*cos(vector_angle);
                Fy += force*((M_PI/2) - vector_angle);
            }
        }
        double resultant_force = sqrt( Fx*Fx + Fy*Fy );
        double resultant_force_angle = acos(Fx/resultant_force);

        forceVector[i].force = resultant_force;
        forceVector[i].angle = resultant_force_angle;
    }
}

