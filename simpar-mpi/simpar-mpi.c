#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005
#define MAX_COORDINATES_VALUE 1.0
#define NUMBER_OF_GHOST_ROWS 2 //Não pode ser por define (e se numero de colunas/processadores <  NUMBER_OF_GHOST_ROWS ?? não permite alteração

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

// s = ncside
#define STARTING_ROW(id,p,n,s) BLOCK_LOW(id,p,n)/s
#define ENDING_ROW(id,p,n,s) STARTING_ROW((id+1),p,n,s)-1
#define ROWS_NUMBER(id,p,n) (ENDING_ROW(id,p,n)-STARTING_ROW(id,p,n)+1)


#define STARTING_ROW_IDX(id,p,n,s) BLOCK_LOW(id,p,n)
#define ENDING_ROW_IDX(id,p,n,s) BLOCK_HIGH(id,p,n)//ENDING_ROW((id),p,n,s) * s
#define NUMBER_OF_ELEMENTS(id,p,n,s) BLOCK_SIZE(id,p,n)//(ENDING_ROW_IDX(id,p,n,s)-STARTING_ROW_IDX(id,p,n,s)+1) // TODO #wrong


int NUMBER_OF_PROCESSES, NUMBER_OF_GHOST_ROWS_2;
long NCSIDE;

MPI_Op reduceOverallCellOp;
MPI_Datatype cellMPIType;
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

typedef struct particle_list {
    struct particle_list * prev;
    struct particle_list * next;
    particle_t * particle;
}particle_list;

typedef struct {
    double x;
    double y;
    double m;
    particle_list * particles;
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

particle_list * rmvList( particle_list * list, particle_t * particle){
    particle_list * curr = list;
    while(curr->particle != particle && curr != NULL){
        curr = curr->next;
    }

    curr->next->prev = curr->prev;
    curr->prev->next = curr->next;

    return curr;
}

particle_list * addList( particle_list * list, particle_t * particle){
    particle_list *newOne = (particle_list *)malloc(sizeof(particle_list));
    newOne->particle = particle;

    newOne->next = list->next;
    newOne->prev = list;
    list->next = newOne;
    if(newOne->next!=NULL)
        newOne->next->prev = newOne;    

    return newOne;
}

/**
 * Update the matrix coordinates of each particle
 *
 * @param particle  A pointer to a particle
 * @param ncside    The number of cells on each side of the matrix of cells
 */
void update_particle(particle_t *particle, cell * cells, long ncside, int processId){
    //Verifies if this cell belongs to this processor
    double sizeCell = MAX_COORDINATES_VALUE/ncside;
    long newCellX = particle->x/sizeCell;
    long newCellY = particle->y/sizeCell;
    long cellIndex = newCellX * ncside + newCellY;

    if((cellIndex >= BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside * ncside)) &&
            (cellIndex <= BLOCK_HIGH(processId, NUMBER_OF_PROCESSES, ncside * ncside))){
        particle->cellX = newCellX;
        particle->cellY = newCellY;
        addList( cells[cellIndex].particles, particle );
    }
}

/**
 * Update the matrix coordinates of each particle
 *
 * @param particle  A pointer to a particle
 * @param ncside    The number of cells on each side of the matrix of cells
 */
void move_particle(particle_t *particle, cell * cells, long ncside, int processId){
    //Verifies if this cell belongs to this processor
    long oldCellIndex = particle->cellX * ncside + particle->cellY;
    double sizeCell = MAX_COORDINATES_VALUE/NCSIDE;
    long newCellX = particle->x/sizeCell;
    long newCellY = particle->y/sizeCell;
    long cellIndex = newCellX * ncside + newCellY;

    if(!(cellIndex >= BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside * ncside) &&
         cellIndex <= BLOCK_HIGH(processId, NUMBER_OF_PROCESSES, ncside * ncside)) ||
         cellIndex == oldCellIndex)return;
    rmvList( cells[oldCellIndex].particles, particle );

    particle->cellX = newCellX;
    particle->cellY = newCellY;
    addList( cells[cellIndex].particles, particle );
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
void clean_cells(cell * cells, long ncside, int processId){
    for (long cellIndex = 0; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside); cellIndex++){
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
 * @param numberOfParticles    The size of the particles array
 * @param ncside    The number of cells on each side of the matrix of cells
 * @return Matrix of cells and update the matrix coordenates of each particle
 */
cell * create_grid(particle_t * particles, long long numberOfParticles, cell * cells, long ncside, int processId){
    for (long cellIndex = 0;
            cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside) + ncside * NUMBER_OF_GHOST_ROWS * 2; cellIndex++)
        cells[cellIndex].particles = (particle_list *) calloc(1, sizeof(particle_list));

    for(long long particleIndex = 0; particleIndex < numberOfParticles; particleIndex++){
        update_particle(&particles[particleIndex], &cells[ncside * NUMBER_OF_GHOST_ROWS], ncside, processId);
    }
}

/**
 * Sends a row to another process
 *
 * @param
 * @param
 * @param
 */
#define SEND_GHOST_ROW_TAG 1
#define TOP_PROCESS_GHOST_ROW_REQUEST_IDX 0
#define BOTTOM_PROCESS_GHOST_ROW_REQUEST_IDX 1
#define TOP_PROCESS_GHOST_ROW_SEND_IDX 2
#define BOTTOM_PROCESS_GHOST_ROW_SEND_IDX 3

/**
 *
 * @param row
 * @param ncside
 * @param topProcessRow
 * @param bottomProcessRow
 * @param senderProcessId
 * @param requests
 * @param statuses
 */
void exchangeGhostRows (cell * cells, long ncside, int senderProcessId, cell *  topProcessRow, cell *  bottomProcessRow,  MPI_Request * requests, MPI_Status * statuses) {
    if (NUMBER_OF_PROCESSES == 1) return; // If it is only one processor, there is no need for communication

    int topProcessId = (senderProcessId - 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;
    int bottomProcessId = (senderProcessId + 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;

    cell * topGhostRow = (cells - (ncside * NUMBER_OF_GHOST_ROWS));
    cell * bottomGhostRow = &cells[NUMBER_OF_ELEMENTS(bottomProcessId, NUMBER_OF_PROCESSES, ncside * ncside, ncside)];

    MPI_Irecv(topGhostRow, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, topProcessId,
              SEND_GHOST_ROW_TAG, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(bottomGhostRow, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, bottomProcessId,
              SEND_GHOST_ROW_TAG, MPI_COMM_WORLD, &requests[1]);

    MPI_Send(cells, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, topProcessId, SEND_GHOST_ROW_TAG, MPI_COMM_WORLD);

    MPI_Send(&cells[NUMBER_OF_ELEMENTS(senderProcessId, NUMBER_OF_PROCESSES, ncside * ncside, ncside) - ncside * NUMBER_OF_GHOST_ROWS],
              ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, bottomProcessId, SEND_GHOST_ROW_TAG, MPI_COMM_WORLD);

    // Ideally use the fact that a row is received to compute the movements on the received side (more parallelism = less time)
    MPI_Waitall(2, requests, statuses);

    // Check the status for a possible error TODO use the count and cancelled to check errors
    for (int statusesIndex = 0; statusesIndex < 2; statusesIndex++)
        printf("Sender Id: %d | Status index: %d | Cancelled: %d | Count: %d\n",
               senderProcessId, statusesIndex, statuses[statusesIndex]._cancelled, statuses[statusesIndex]._ucount);

    printf("SID %d TOP||  X:%0.2f; Y:%0.2f; M:%0.2f\n", senderProcessId, topGhostRow[0].x,topGhostRow[0].y,topGhostRow[0].m);
    printf("SID %d BOT||  X:%0.2f; Y:%0.2f; M:%0.2f\n", senderProcessId, bottomGhostRow[0].x,bottomGhostRow[0].y,bottomGhostRow[0].m);
}

/**
 * Computes the center mass for each existing cell in grid
 *
 * @param particles     array of particles that compose the grid
 * @param length        number of particles (must match @particles length)
 * @param cells         bidimensional array with the grid of cells
 * @param ncside        sides of the grid (how many rows the grid has)
 * @param processId    id of the process
 */
void compute_cell_center_mass(particle_t *particles, long length, cell * cells, long ncside, int processId) {
    clean_cells(cells, ncside, processId);

    for (long cellIndex = 0; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside);
    cellIndex++){
        particle_list * currentParticleList = cells[cellIndex].particles->next;
        if(currentParticleList != NULL){ //avoid div by 0 after while
            do {
                particle_t * currentParticle = currentParticleList->particle;
                cells[cellIndex].x += currentParticle->x * currentParticle->m;
                cells[cellIndex].y += currentParticle->y * currentParticle->m;
                cells[cellIndex].m += currentParticle->m;
                currentParticleList = currentParticleList->next;
            } while (currentParticleList != NULL);
            cells[cellIndex].x /= cells[cellIndex].m;
            cells[cellIndex].y /= cells[cellIndex].m;
        }
    }

    printf("PID %d TOP|| X:%0.2f; Y:%0.2f; M:%0.2f\n", processId, cells[0].x,cells[0].y,cells[0].m);
    printf("PID %d BOT|| X:%0.2f; Y:%0.2f; M:%0.2f\n", processId,
            cells[NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside)-1].x,
            cells[NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside)-1].y,
            cells[NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside)-1].m);

    cell topProcessRow, bottomProcessRow;
    MPI_Status statuses[4];
    MPI_Request requests[4];

    exchangeGhostRows(cells, ncside, processId,
            &topProcessRow, &bottomProcessRow, &requests[0], &statuses[0]);
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
 * @param cells     Cells matrix
 * @param ncside    number of rows/columns of the matrix
 * @param processId ID of the current process
 */
void update_particle_position(particle_t * particle, double Fx, double Fy, cell * cells, long ncside, int processId){
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

    move_particle(particle, cells, ncside, processId);
}

/**
 * Function that gets the cell with the coordinates (unbound_row, unbound_column). Those coordinates will always represent
 * a cell because the space is wrapped, that is, if a particle exits through one side of the space, it enters on the corresponding
 * opposite side of the space.
 *
 * @param unbound_row       index unbounded in X axis
 * @param unbound_column    index unbounded in Y axis
 * @param cells             matrix with the cells to search
 * @param return_cell       cell used to represent the unbounded cell
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 * @return                  the neighbor cell
 */
cell * get_cell(long long unbounded_row, long long unbounded_column, cell *cells, cell * return_cell, long ncside){

    if (unbounded_column >= 0 && unbounded_column < ncside && unbounded_row >= 0 && unbounded_row < ncside){
        return &cells[unbounded_row * ncside + unbounded_column];
    } else {
        cell * matrixWithGhostRows = cells - (ncside * NUMBER_OF_GHOST_ROWS);
        long bounded_row = (unbounded_row + ncside) % ncside;
        /*
         *  -------
         * |  TOP  | Since top row is BEFORE the data, any reference to it will result in a negative column index
         * |-------|
         * |  DATA |
         * |-------|
         * |  BOT  | Since bot row is AFTER the data, any reference to it will result in a column index higher than its possible
         *  -------
         */
        long bounded_column = (unbounded_column < 0) ? NUMBER_OF_GHOST_ROWS + unbounded_column : unbounded_column;
        long cellIndex = bounded_row * ncside + bounded_column;

        return_cell->x = matrixWithGhostRows[cellIndex].x;
        return_cell->y = matrixWithGhostRows[cellIndex].y;
        return_cell->m = matrixWithGhostRows[cellIndex].m;

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
 * @param particlesLength  number of particles in the list
 * @param cells             matrix that represents the cells
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void compute_force_and_update_particles(cell *cells, long ncside, int processId){
    for(long cellIndex = 0;
            cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside); cellIndex++) {
        particle_list *cellParticles = cells[cellIndex].particles;
        cell neighboursList[9];

        long cellX = ceil(cellIndex / ncside);
        long cellY = cellIndex % ncside;

        //resultant force in X and Y
        double Fx = 0;
        double Fy = 0;

        int neighboursListIndex = 0;
        for (int row = -1; row < 2; row++) {
            for (int column = -1; column < 2; column++) {
                get_cell(row + cellX, column + cellY, cells, &neighboursList[neighboursListIndex++], ncside);
            }

            while (cellParticles->next != 0) {
                particle_t *currentParticle = cellParticles->particle;
                for (neighboursListIndex = 0; neighboursListIndex < 9; neighboursListIndex++) {
                    //if that cell doesn't have any particle, no central mass will exist, so we ignore
                    if (neighboursList[neighboursListIndex].m == 0)
                        continue;

                    //compute angle
                    double delta_x = neighboursList[neighboursListIndex].x - currentParticle->x;
                    double delta_y = neighboursList[neighboursListIndex].y - currentParticle->y;
                    double vector_angle = atan2(delta_y, delta_x);

                    //compute force
                    double force = compute_magnitude_force(currentParticle, &neighboursList[neighboursListIndex]);
                    Fx += force * cos(vector_angle);
                    Fy += force * sin(vector_angle);
                }

                update_particle_position(currentParticle, Fx, Fy, cells, ncside, processId);
                cellParticles = cellParticles->next;
            }
        }
    }



   /* //iterate all particles
    for(long particleIndex = BLOCK_LOW(processId, NUMBER_OF_PROCESSES, particlesLength);
            particleIndex <= BLOCK_HIGH(processId, NUMBER_OF_PROCESSES, particlesLength);
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
        update_particle_position(&particles[particleIndex], Fx, Fy, cells, ncside, processId);
    }*/
}

void reduceOverallCellsMatrix (cell * in, cell * out, int *len , MPI_Datatype *datatype){
    out->x += in->x;
    out->y += in->y;
    out->m += in->m;
}


/**
 * Computes the center mass for each existing cell in grid
 *
 * @param cells     array of particles that compose the grid
 * @param ncside        number of particles (must match @particles length)
 */
void compute_overall_center_mass(cell * cells, long ncside, int processId){
    cell overallCenterMass = {.x=0, .y=0, .m=0};

    for(long cellIndex = 0; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside * ncside, ncside);
    cellIndex++) {
        particle_list * currentParticleList = cells[cellIndex].particles;
        while (currentParticleList->particle != 0){
            particle_t * currentParticle = currentParticleList->particle;
            overallCenterMass.x += currentParticle->x * currentParticle->m;
            overallCenterMass.y += currentParticle->y * currentParticle->m;
            overallCenterMass.m += currentParticle->m;
        }
    }

    cell outOverallCenterMass = {.x=0, .y=0, .m=0};
    printf("Process id: %d  |  %0.2f %0.2f %0.2f \n", processId, overallCenterMass.x, overallCenterMass.y, overallCenterMass.m);

    MPI_Reduce(&overallCenterMass, &outOverallCenterMass, 1, cellMPIType, reduceOverallCellOp, 0, MPI_COMM_WORLD);

    outOverallCenterMass.x /= outOverallCenterMass.m;
    outOverallCenterMass.y /= outOverallCenterMass.m;

    if(processId == 0)
        printf("%0.2f %0.2f \n", outOverallCenterMass.x, outOverallCenterMass.y);
}


//TODO alterar o metodo para a nova celula
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
    MPI_Type_struct(2, blocklens, indices, oldTypes, newType);
    MPI_Type_commit(newType);
}


int main(int args_length, char* args[]) {
    if (args_length < 4) {
        printf("Incorrect number of arguments! It should be 4!");
        return -1;
    }

    int factor = 0;
    if(args_length == 8) factor = 3;

    long seed = strtol(args[1 + factor], NULL, 10);
    long ncside = NCSIDE = strtol(args[2 + factor], NULL, 10);
    long long n_part = strtol(args[3 + factor], NULL, 10);
    long iterations = strtol(args[4 + factor], NULL, 10);

    particle_t * particles = calloc(n_part, sizeof(particle_t));

    int rank;
    MPI_Init( &args_length, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMBER_OF_PROCESSES);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //Type & Operations Setup
    mapCellToMPI(&cellMPIType);
    mapParticleToMPI(&particleMPIType);
    MPI_Op_create(reduceOverallCellsMatrix, 1, &reduceOverallCellOp);

    /*
     * CellMatrix organization:
     *    ----------------------------------------------------------
     *   |  Top Ghost Row  |  Process's cells  |  Bottom Ghost Row  |
     *    ----------------------------------------------------------
     */
    printf("Allocating for %d, %d cells\n", rank, NUMBER_OF_ELEMENTS(rank, NUMBER_OF_PROCESSES, ncside * ncside, ncside) + (ncside * NUMBER_OF_GHOST_ROWS * 2));
    cell * cellMatrix = (cell*) calloc(NUMBER_OF_ELEMENTS(rank, NUMBER_OF_PROCESSES, ncside * ncside, ncside) + (ncside * NUMBER_OF_GHOST_ROWS * 2), sizeof(cell));

    if(rank==0){
        init_particles(seed, ncside, n_part, particles);
    }


    MPI_Bcast(particles, n_part, particleMPIType, 0, MPI_COMM_WORLD);

    create_grid(particles, n_part, cellMatrix, ncside, rank);


//    for(int i = 0; i < n_part; i++)
//        printf("Process id: %d | Particula %d | %0.2f %0.2f %0.2f \n", rank, i, particles[i].x, particles[i].y, particles[i].m);
    /*
    for(long particleIndex = BLOCK_LOW(me, NUMBER_OF_PROCESSES, n_part); 
            particleIndex <= BLOCK_HIGH(me, NUMBER_OF_PROCESSES, n_part); 
            particleIndex++){
        printf("Process id: %d | Particula %d | %0.2f %0.2f %0.2f \n", me, particleIndex, particles[particleIndex].x, particles[particleIndex].y, particles[particleIndex].m);
    }*/

    for(int i = 0; i < iterations; i++){
        //To simplify the index treatment (to start with 0) the cells matrix that goes through parameter omits the top ghost row
        compute_cell_center_mass(particles, n_part, &cellMatrix[ncside * NUMBER_OF_GHOST_ROWS], ncside, rank);
//        compute_force_and_update_particles(cellMatrix, ncside, rank);
    }


    if(rank==0){
        printf("%0.2f %0.2f \n", particles[0].x, particles[0].y);
    }
//    compute_overall_center_mass(cellMatrix, ncside, rank);

    MPI_Finalize();

/**
    printf("END\n");    
    for(int i = 0; i < n_part; i++){
        printf("%d Process | Particle %d %0.2f %0.2f\n", me, i, particles[0].x, particles[0].y);
    }

    printf("END\n");
    */
   
    free(cellMatrix);
    free(particles);
    
    
    
}
