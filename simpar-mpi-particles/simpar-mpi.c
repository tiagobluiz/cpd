#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005
#define MAX_COORDINATES_VALUE 1.0

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)-1)

int NUMBER_OF_PROCESSES;
long NCSIDE, MAX_BLOCK_SIZE=0;

MPI_Op reduceCellOp;
MPI_Datatype particleMPIType;


typedef struct {
    double x;
    double y;
    double vx;
    double vy;
    double m;
    long cellX;
    long cellY;  
}particle_t;

typedef struct {
    double x;
    double y;
    double m;
}cell;

/**
 * Utility Methods
 */

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
void clean_cells(cell * cells, long ncside){
    for (long cellIndex = 0; cellIndex < ncside * ncside; cellIndex++){
        cells[cellIndex].x = cells[cellIndex].y =
        cells[cellIndex].m = 0;
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
cell * create_grid(particle_t * particles, long long length, long ncside){
    for(int i = 0; i < length; i++){
        update_particle(&particles[i], ncside);
    }
}

/**
 * Reduces two arrays of cells into one (reduction function of redcell)
 *
 * @param a     accumulated array
 * @param b     incoming array, which will be deallocated from memory
 */
void reduceCellsMatrix (cell * in, cell * out, int *len , MPI_Datatype *datatype){
    for (long cellIndex = 0; cellIndex < NCSIDE * NCSIDE; cellIndex++){
            out[cellIndex].x += in[cellIndex].x;
            out[cellIndex].y += in[cellIndex].y;
            out[cellIndex].m += in[cellIndex].m;
    }
    free(in);
}


/**
 * Computes the center mass for each existing cell in grid
 *
 * @param particles     array of particles that compose the grid
 * @param length        number of particles (must match @particles length)
 * @param cells         bidimensional array with the grid of cells
 * @param ncside        sides of the grid (how many rows the grid has)
 * @param process_id    id of the process
 */
void compute_cell_center_mass(particle_t *particles, long length, cell * cells, long ncside, int process_id, MPI_Datatype datatype) { 
    clean_cells(cells, ncside);
    cell cellLocalMatrix[ncside * ncside];
    memcpy(cellLocalMatrix, cells, ncside * ncside);

    for(long particleIndex = BLOCK_LOW(process_id, NUMBER_OF_PROCESSES, length); 
            particleIndex < BLOCK_SIZE(process_id, NUMBER_OF_PROCESSES, length); 
            particleIndex++){
        particle_t particle = particles[particleIndex];
        cellLocalMatrix[particle.cellX * ncside + particle.cellY].x += particle.m * particle.x;
        cellLocalMatrix[particle.cellX * ncside + particle.cellY].y += particle.m * particle.y;
        cellLocalMatrix[particle.cellX * ncside + particle.cellY].m += particle.m;
    }

    MPI_Allreduce(cellLocalMatrix, cells, ncside * ncside, datatype, reduceCellOp, MPI_COMM_WORLD);

    for (long cellIndex = 0; cellIndex < ncside * ncside; cellIndex++){
        cells[cellIndex].x /= cells[cellIndex].m;
        cells[cellIndex].y /= cells[cellIndex].m;
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
cell * get_cell(long long unbounded_row, long long unbounded_column, cell *cells, cell * return_cell, long ncside){
    int bounded_row = mod(unbounded_row, ncside);
    int bounded_column = mod(unbounded_column, ncside);

    if (unbounded_column >= 0 && unbounded_column < ncside && unbounded_row >= 0 && unbounded_row < ncside){
        return &cells[bounded_row * ncside + bounded_column];
    } else {
        return_cell->x = cells[bounded_row * ncside + bounded_column].x;
        return_cell->y = cells[bounded_row * ncside + bounded_column].y;
        return_cell->m = cells[bounded_row * ncside + bounded_column].m;

        if (unbounded_column < 0)
            return_cell->x -= MAX_COORDINATES_VALUE;
        else if (unbounded_column >= ncside)
            return_cell->x += MAX_COORDINATES_VALUE;

        if (unbounded_row < 0)
            return_cell->y -= MAX_COORDINATES_VALUE;
        else if (unbounded_row >= ncside)
            return_cell->y += MAX_COORDINATES_VALUE;
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
void compute_force_and_update_particles(particle_t *particles, int particles_length, cell *cells, long ncside, int process_id){

    particle_t particleCopy[MAX_BLOCK_SIZE];
    long copyIndex =0;
    //iterate all particles
    for(long particleIndex = BLOCK_LOW(process_id, NUMBER_OF_PROCESSES, particles_length); 
            particleIndex < BLOCK_SIZE(process_id, NUMBER_OF_PROCESSES, particles_length); 
            particleIndex++){

        //get the coordinates of the cell where the particle is located
        long long cell_x = particles[particleIndex].cellX;
        long long cell_y = particles[particleIndex].cellY;

        //resultant force in X and Y
        double Fx = 0;
        double Fy = 0;

        // iterate all the 9 mass cells around the particle
        for(int row = -1; row < 2; row++){
            for(int column = -1; column < 2; column++){
                //get the neighbor cell in the coordinates cellMatrix = create_grid(particles, n_part, ncside);(row, column)
                cell created_neighbor_cell;
                cell * neighbor_cell = get_cell(row + cell_x, column + cell_y, cells, &created_neighbor_cell, ncside);

                //if that cell doesn't have any particle, no central mass will exist, so we ignore
                if(neighbor_cell->m == 0)
                    continue;

                //compute angle
                double delta_x = neighbor_cell->x - particles[particleIndex].x;
                double delta_y = neighbor_cell->y - particles[particleIndex].y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force
                double force = compute_magnitude_force(&particles[particleIndex], neighbor_cell);
                Fx += force * cos(vector_angle);
                Fy += force * sin(vector_angle);
            }
        }

        //update the particle position
        update_particle_position(&particles[particleIndex], Fx, Fy, ncside);   
        particleCopy[copyIndex++] = particles[particleIndex];     
    }
    particle_t particleReceive[MAX_BLOCK_SIZE];

    for(int i =0; i<NUMBER_OF_PROCESSES; i++){
        if(i!=process_id){
            MPI_Bcast(particleReceive, BLOCK_SIZE(process_id, NUMBER_OF_PROCESSES, particles_length), particleMPIType, process_id, MPI_COMM_WORLD);
            memcpy(&particles[BLOCK_LOW(i, NUMBER_OF_PROCESSES, particles_length)], particleReceive, BLOCK_SIZE(i, NUMBER_OF_PROCESSES, particles_length));
        }
        else
            MPI_Bcast(particleCopy, BLOCK_SIZE(process_id, NUMBER_OF_PROCESSES, particles_length), particleMPIType, process_id, MPI_COMM_WORLD);
    }
}

/**
 * Computes the center mass for each existing cell in grid
 *
 * @param particles     array of particles that compose the grid
 * @param length        number of particles (must match @particles length)
 */
void compute_overall_center_mass(particle_t * particles, long length){
    cell overallCenterMass = {.x=0, .y=0, .m=0};
    for (long particleIndex = 0; particleIndex < length; particleIndex++){
        overallCenterMass.x += particles[particleIndex].x * particles[particleIndex].m;
        overallCenterMass.y += particles[particleIndex].y * particles[particleIndex].m;
        overallCenterMass.m += particles[particleIndex].m;
    }
    overallCenterMass.x /= overallCenterMass.m;
    overallCenterMass.y /= overallCenterMass.m;
    printf("%0.2f %0.2f \n", overallCenterMass.x, overallCenterMass.y);
}

void mapCellToMPI(MPI_Datatype * newType){
    MPI_Type_contiguous(3, MPI_DOUBLE, newType);
    MPI_Type_commit(newType);
}

void mapParticleToMPI(MPI_Datatype * newType){
    int blocklens[] = {5 /*doubles*/,2 /*long long*/};
    MPI_Aint extent;
    MPI_Type_extent(MPI_DOUBLE, &extent);
    MPI_Aint indices[] = {0, 5 * extent /* we have 5 doubles */};
    MPI_Datatype oldTypes[] = {MPI_DOUBLE, MPI_LONG};
    MPI_Type_struct(1, blocklens, indices, oldTypes, newType);
    MPI_Type_commit(newType);
}


int main(int args_length, char* args[]) {
    if (args_length < 4) {
        printf("Incorrect number of arguments! It should be 4!");
        return -1;
    }

    long seed = strtol(args[1], NULL, 10);
    long ncside = NCSIDE = strtol(args[2], NULL, 10);
    long long n_part = strtol(args[3], NULL, 10);
    long iterations = strtol(args[4], NULL, 10);

    particle_t * particles = calloc(n_part, sizeof(particle_t));

    int me;
    MPI_Init( &args_length, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMBER_OF_PROCESSES);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    //Type & Operations Setup
    MPI_Op_create(reduceCellsMatrix, 1, &reduceCellOp);
    MPI_Datatype cellMPIType;
    mapCellToMPI(&cellMPIType);
    mapParticleToMPI(&particleMPIType);

    cell * cellMatrix = (cell*) calloc(ncside * ncside, sizeof(cell));

    if(me==0){
        init_particles(seed, ncside, n_part, particles); 
        create_grid(particles, n_part, ncside);
    }

    MPI_Bcast(cellMatrix, ncside * ncside, cellMPIType, 0, MPI_COMM_WORLD);
    MPI_Bcast(particles, n_part, particleMPIType, 0, MPI_COMM_WORLD);    

    for(int i = 0; i < NUMBER_OF_PROCESSES;i++){
        long blockSize = BLOCK_SIZE(i, NUMBER_OF_PROCESSES, n_part);
        if(blockSize>MAX_BLOCK_SIZE)
            MAX_BLOCK_SIZE = blockSize;
    }

    for(int i = 0; i < iterations; i++){
        compute_cell_center_mass(particles, n_part, cellMatrix, ncside, me, cellMPIType);
        compute_force_and_update_particles(particles, n_part, cellMatrix, ncside, me);  
    }

    MPI_Finalize();
    if(me==0){
        printf("%0.2f %0.2f \n", particles[0].x, particles[0].y);
        compute_overall_center_mass(particles, n_part);
    }
    free(cellMatrix);
    free(particles);
    
    
    
}
