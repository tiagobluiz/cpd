#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define RND0_1 ((double) random() / ((long long)1<<31))
#define G 6.67408e-11
#define EPSLON 0.0005
#define MAX_COORDINATES_VALUE 1.0
#define NUMBER_OF_GHOST_ROWS 1

#define BLOCK_LOW(id,p,n)  ceil(((id)*(n))/(p))*n
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define NUMBER_OF_ELEMENTS(id,p,n) (long)(BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)


int NUMBER_OF_PROCESSES;
long TOTAL_ELEMENTS;

MPI_Op reduceOverallCellOp;
MPI_Datatype cellMPIType;
MPI_Datatype particleMPIType;


typedef struct {
    double x;
    double y;
    double vx;
    double vy;
    double m;
    long arrayIndex;
    long alreadyMoved; //Flag to check if a particle was already verified
    long creationIndex;
}particle_t;

typedef struct {
    double x;
    double y;
    double m;
    long long nParticles;
    long long allocatedSpace;
    particle_t * particles;
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

        //our
        par[i].creationIndex = i;
    }
}

void rmvList( cell * cell, particle_t * particle, long oldParticleArrayIndex){
    if(particle->creationIndex != cell->particles[oldParticleArrayIndex].creationIndex)
        printf("\t\t/!\\ A particle was not removed as it should be\n");
    if(oldParticleArrayIndex + 1 < cell->nParticles){
        memcpy(&cell->particles[oldParticleArrayIndex], &cell->particles[cell->nParticles - 1], sizeof(particle_t));
        cell->particles[oldParticleArrayIndex].arrayIndex = oldParticleArrayIndex;
    }

    cell->nParticles--;

    // If we are only using 1/3 of the space, we can realloc to use less space (just cut at half to allow some movement)
    if(cell->nParticles < cell->allocatedSpace/3){
        cell->allocatedSpace /= 2;
        particle_t * buffer = (particle_t * )realloc(cell->particles, cell->allocatedSpace * sizeof(particle_t));
        if(buffer == NULL){
            printf("/!\\ There was an error on realloc, the execution is corrupt\n");
            return;
        }
        cell->particles = buffer;
    }
}

void addList(cell * cell, particle_t * particle){
    if( cell->nParticles + 1 > cell->allocatedSpace ){
        if(cell->allocatedSpace > 0){
            cell->allocatedSpace *= 2;
            particle_t * buffer = (particle_t * )realloc(cell->particles, cell->allocatedSpace * sizeof(particle_t));
            if(buffer == NULL){
                printf("/!\\ There was an error on realloc, the execution is corrupt\n");
                return;
            }
            cell->particles = buffer;
        } else {
            cell->allocatedSpace = 1;
            cell->particles = (particle_t * )calloc(cell->allocatedSpace, sizeof(particle_t));
        }
    }

    particle->arrayIndex = cell->nParticles;
    cell->particles[cell->nParticles] = *particle;
    cell->nParticles++;
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
    if(newCellX > ncside || newCellY > ncside)
        printf("A particle was positioned out of matrix (particle coordinates >= 1.0)\n");

    long globalCellIndex = newCellY * ncside + newCellX;

    if(!(globalCellIndex >= BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside) &&
         globalCellIndex <= BLOCK_HIGH(processId, NUMBER_OF_PROCESSES, ncside))) return;

    int localCellIndex = globalCellIndex - BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside);
    addList( &cells[localCellIndex], particle );
}

/**
 * Update the matrix coordinates of each particle
 *
 * @param particle  A pointer to a particle
 * @param ncside    The number of cells on each side of the matrix of cells
 */
int move_particle(particle_t *particle, cell * cells, long ncside, long oldCellIndex, int particleIndex, int processId){

//    if(particle->y < 0 || particle->y >= MAX_COORDINATES_VALUE)
//        printf("--> %d ADDING %d to CI %d from CI %d\n", processId, particle->creationIndex, 0, oldCellIndex);

    //Verifies if this cell belongs to this processor
    double sizeCell = MAX_COORDINATES_VALUE/ncside;
    int y = floor(particle->y / sizeCell);
    int x = floor(particle->x / sizeCell);
    long globalCellIndex = y*ncside + x;




    long localCellIndex = globalCellIndex - BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside) + (ncside * NUMBER_OF_GHOST_ROWS);


    if(localCellIndex == oldCellIndex){
        return 0;
    }
    
    if(particle->y >= MAX_COORDINATES_VALUE) particle->y -= MAX_COORDINATES_VALUE;
    else if (particle->y < 0) particle->y += MAX_COORDINATES_VALUE;

    addList(&cells[localCellIndex], particle);

    rmvList( &cells[oldCellIndex], particle, particleIndex);

    return 1;
}

/**
 * Cleans all cells and resets all parameters to 0 so that the previous state doesn't affect the new one
 * @param cells     bidimensional array with the grid of cells
 * @param ncside    sides of the grid (how many rows the grid has)
 */
void clean_cells(cell * cells, long ncside, int processId){
    #pragma omp parallel for
    for (long cellIndex = 0; cellIndex < TOTAL_ELEMENTS; cellIndex++){
        cells[cellIndex].x = cells[cellIndex].y =
        cells[cellIndex].m = 0;
        if(cellIndex < ncside * NUMBER_OF_GHOST_ROWS || cellIndex >= TOTAL_ELEMENTS - (ncside * NUMBER_OF_GHOST_ROWS)){
            if(cells[cellIndex].nParticles > 0) free(cells[cellIndex].particles);
            cells[cellIndex].nParticles = cells[cellIndex].allocatedSpace = 0;
        }
        #pragma omp parallel for
        for (long particleIndex = 0; particleIndex < cells[cellIndex].nParticles; particleIndex++)
            cells[cellIndex].particles[particleIndex].alreadyMoved = 0;
    }
}

/**
 * Implementation
 */


/**
 * Sends rows to top & bottom processes (whose id is -1 and +1 of senderProcessId)
 *
 * @param cells             matrix of this process' cells (w/o ghost rows)
 * @param ncside            number of sides of the cell matrix
 * @param senderProcessId   id of the process performing the send operation
 */

#define SEND_TOP_GHOST_ROW_TAG 1
#define SEND_BOT_GHOST_ROW_TAG 2
#define SEND_TOP_GHOST_PARTICLES_TAG 3
#define SEND_BOT_GHOST_PARTICLES_TAG 4
void exchangeGhostRowsCells(cell *cells, long ncside, int senderProcessId) {
    if (NUMBER_OF_PROCESSES == 1) return; // If it is only one processor, there is no need for communication

    int topProcessId = (senderProcessId - 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;
    int bottomProcessId = (senderProcessId + 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;

    MPI_Status statuses[2];
    MPI_Request requests[2];

    cell *topGhostRow = (cells - (ncside * NUMBER_OF_GHOST_ROWS));
    cell *bottomGhostRow = cells + NUMBER_OF_ELEMENTS(senderProcessId, NUMBER_OF_PROCESSES, ncside);

    MPI_Irecv(topGhostRow, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, topProcessId,
              SEND_BOT_GHOST_ROW_TAG, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(bottomGhostRow, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, bottomProcessId,
              SEND_TOP_GHOST_ROW_TAG, MPI_COMM_WORLD, &requests[1]);


    MPI_Send(cells, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, topProcessId,
             SEND_TOP_GHOST_ROW_TAG, MPI_COMM_WORLD);
    MPI_Send(bottomGhostRow - (ncside * NUMBER_OF_GHOST_ROWS), ncside * NUMBER_OF_GHOST_ROWS, cellMPIType,
             bottomProcessId,
             SEND_BOT_GHOST_ROW_TAG, MPI_COMM_WORLD);


    MPI_Waitall(2, requests, statuses);

    #pragma omp parallel for
    for(long cellIndex = 0; cellIndex < ncside * NUMBER_OF_GHOST_ROWS; cellIndex++){
        topGhostRow[cellIndex].allocatedSpace = topGhostRow[cellIndex].nParticles =
                bottomGhostRow[cellIndex].allocatedSpace = bottomGhostRow[cellIndex].nParticles = 0;

        if (senderProcessId == 0) topGhostRow[cellIndex].y -= MAX_COORDINATES_VALUE;
        else if (senderProcessId == NUMBER_OF_PROCESSES - 1) bottomGhostRow[cellIndex].y += MAX_COORDINATES_VALUE;
    }
}

void exchangeGhostRowsParticles(cell * cells, long ncside, int senderProcessId){
//    printf("PID %d exchanging particles\n", senderProcessId);
    if (NUMBER_OF_PROCESSES == 1) return; // If it is only one processor, there is no need for communication

    int topProcessId = (senderProcessId - 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;
    int bottomProcessId = (senderProcessId + 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;

    cell *topGhostRow = cells;
    cell *bottomGhostRow = cells + (ncside * NUMBER_OF_GHOST_ROWS) + NUMBER_OF_ELEMENTS(senderProcessId, NUMBER_OF_PROCESSES, ncside);

    MPI_Status statuses[2];
    MPI_Request requests[2];


//    printf("PID %d before counting\n", senderProcessId);
    //count all particles that exists in frontier rows
    long long countsTopParticlesToSendByCell[ncside + 1], countsBotParticlesToSendByCell[ncside + 1];
    countsBotParticlesToSendByCell[0] = countsTopParticlesToSendByCell[0] = 0; //Ensure that position 0 starts with 0
    cell * topParticlesToSend = cells ;//+ (ncside * (NUMBER_OF_GHOST_ROWS - 1));
    cell * botParticlesToSend = bottomGhostRow;
    for(long cellIndex = 0; cellIndex < ncside; cellIndex++){
//        printf("%d -> TOP %d BOT %d\n", senderProcessId, topParticlesToSend[cellIndex].nParticles, botParticlesToSend[cellIndex].nParticles);
        countsTopParticlesToSendByCell[0] += topParticlesToSend[cellIndex].nParticles;
        countsBotParticlesToSendByCell[0] += botParticlesToSend[cellIndex].nParticles;

        countsTopParticlesToSendByCell[cellIndex + 1] = topParticlesToSend[cellIndex].nParticles;
        countsBotParticlesToSendByCell[cellIndex + 1] = botParticlesToSend[cellIndex].nParticles;
    }

    particle_t *topParticlesToSendBuffer = calloc(countsTopParticlesToSendByCell[0], sizeof(particle_t));
    particle_t *downParticlesToSendBuffer = calloc(countsBotParticlesToSendByCell[0], sizeof(particle_t));
    long long accumulatorTopParticlesToSend = 0;
    long long accumulatorBotParticlesToSend = 0;

//    printf("PID %d before memcopys\n", senderProcessId);
    for(int cellIndex = 0; cellIndex < ncside; cellIndex++) {
//        printf("%d || TOPAC %d NP %d\n", senderProcessId, accumulatorTopParticlesToSend, topParticlesToSend[cellIndex].nParticles);
        memcpy(topParticlesToSendBuffer + accumulatorTopParticlesToSend,
               topParticlesToSend[cellIndex].particles, topParticlesToSend[cellIndex].nParticles * sizeof(particle_t));
        accumulatorTopParticlesToSend += topParticlesToSend[cellIndex].nParticles;

//        printf("%d || BOTAC %d NP %d\n", senderProcessId, accumulatorBotParticlesToSend, botParticlesToSend[cellIndex].nParticles);
        memcpy(downParticlesToSendBuffer + accumulatorBotParticlesToSend,
               botParticlesToSend[cellIndex].particles,
               botParticlesToSend[cellIndex].nParticles * sizeof(particle_t));
        accumulatorBotParticlesToSend += botParticlesToSend[cellIndex].nParticles;

    }


    long long * countsTopParticlesToReceiveByCell = calloc(ncside + 1, sizeof(int));
    long long * countsBotParticlesToReceiveByCell = calloc(ncside + 1, sizeof(int));// [ncside + 1];

//    printf("PID %d before counts\n", senderProcessId);

//    for (int x = 0; x < ncside + 1; x++)
//        printf("%d || %d ; %d\n", senderProcessId, countsTopParticlesToSendByCell[x], countsBotParticlesToSendByCell[x]);

    //Send counters
    MPI_Irecv(countsTopParticlesToReceiveByCell, ncside + 1, MPI_INT, topProcessId,
              SEND_BOT_GHOST_PARTICLES_TAG, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(countsBotParticlesToReceiveByCell, ncside + 1, MPI_INT, bottomProcessId,
              SEND_TOP_GHOST_PARTICLES_TAG, MPI_COMM_WORLD, &requests[1]);

    MPI_Send(countsTopParticlesToSendByCell, ncside + 1, MPI_INT, topProcessId,
             SEND_TOP_GHOST_PARTICLES_TAG, MPI_COMM_WORLD);
    MPI_Send(countsBotParticlesToSendByCell, ncside + 1, MPI_INT, bottomProcessId,
             SEND_BOT_GHOST_PARTICLES_TAG, MPI_COMM_WORLD);

    MPI_Waitall(2, requests, statuses);

//    for (int x = 0; x < ncside + 1; x++)
//        printf("  %d || %d ; %d\n", senderProcessId, countsTopParticlesToReceiveByCell[x], countsBotParticlesToReceiveByCell[x]);
//
//    printf("PID %d after counts\n", senderProcessId);

    particle_t *topParticlesToReceiveBuffer = calloc(countsTopParticlesToReceiveByCell[0], sizeof(particle_t));
    particle_t *downParticlesToReceiveBuffer = calloc(countsBotParticlesToReceiveByCell[0], sizeof(particle_t));

    //Send data (particles)
    if(countsTopParticlesToReceiveByCell[0] > 0)
        MPI_Irecv(topParticlesToReceiveBuffer, countsTopParticlesToReceiveByCell[0], particleMPIType, topProcessId,
             SEND_BOT_GHOST_PARTICLES_TAG, MPI_COMM_WORLD, &requests[0]);
    if(countsBotParticlesToReceiveByCell[0] > 0)
        MPI_Irecv(downParticlesToReceiveBuffer, countsBotParticlesToReceiveByCell[0], particleMPIType, bottomProcessId,
             SEND_TOP_GHOST_PARTICLES_TAG, MPI_COMM_WORLD, &requests[1]);

    if(countsTopParticlesToReceiveByCell[0] > 0)
        MPI_Send(topParticlesToSendBuffer, countsTopParticlesToSendByCell[0], particleMPIType, topProcessId,
              SEND_TOP_GHOST_PARTICLES_TAG, MPI_COMM_WORLD);
    if(countsBotParticlesToSendByCell[0] > 0)
        MPI_Send(downParticlesToSendBuffer, countsBotParticlesToSendByCell[0], particleMPIType, bottomProcessId,
              SEND_BOT_GHOST_PARTICLES_TAG, MPI_COMM_WORLD);

    MPI_Waitall(2, requests, statuses);


//    printf("PID %d before particles\n", senderProcessId);
    long long topParticlesReceived = 0;
    long long botParticlesReceived = 0;
    cell * topParticlesToReceive = cells + (ncside * NUMBER_OF_GHOST_ROWS);
    cell * botParticlesToReceive = bottomGhostRow - ncside;

    for(long cellIndex = 0; cellIndex < ncside; cellIndex++){
        long long currCellTopParticlesToReceive = countsTopParticlesToReceiveByCell[cellIndex + 1];
        long long currCellBotParticlesToReceive = countsBotParticlesToReceiveByCell[cellIndex + 1];

//        printf("CTP %d CBP %d || ACT %d ACB %d\n", currCellTopParticlesToReceive, currCellBotParticlesToReceive, topParticlesReceived, botParticlesReceived);
        if(countsTopParticlesToReceiveByCell[0] > 0){
            for(int particleIndex = topParticlesReceived; particleIndex < currCellTopParticlesToReceive + topParticlesReceived; particleIndex++)
                addList(&topParticlesToReceive[cellIndex], &downParticlesToReceiveBuffer[particleIndex]);
            topParticlesReceived += currCellTopParticlesToReceive;
        }


        if(countsBotParticlesToReceiveByCell[0] > 0){
            for(int particleIndex = botParticlesReceived; particleIndex < currCellBotParticlesToReceive + botParticlesReceived; particleIndex++)
                addList(&botParticlesToReceive[cellIndex], &downParticlesToReceiveBuffer[particleIndex]);
            botParticlesReceived += currCellBotParticlesToReceive;
        }
    }

    free(topParticlesToReceiveBuffer);
    free(downParticlesToReceiveBuffer);
    free(topParticlesToSendBuffer);
    free(downParticlesToSendBuffer);
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
void compute_cell_center_mass(cell * cells, long ncside, int processId) {
//    clean_cells(cells, ncside, processId);

    for (long cellIndex = 0; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside);
         cellIndex++){
        cell currCell = cells[cellIndex];
        if(currCell.nParticles > 0){ //avoid div by 0
            for(long long particleIndex = 0; particleIndex < currCell.nParticles; particleIndex++) {
                particle_t * currentParticle = &currCell.particles[particleIndex];
                cells[cellIndex].x += currentParticle->x * currentParticle->m;
                cells[cellIndex].y += currentParticle->y * currentParticle->m;
                cells[cellIndex].m += currentParticle->m;
            }
            cells[cellIndex].x /= cells[cellIndex].m;
            cells[cellIndex].y /= cells[cellIndex].m;
        }
    }

    exchangeGhostRowsCells(cells, ncside, processId);
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
int update_particle_position(particle_t * particle, double Fx, double Fy, cell * cells, long ncside, int processId, long cellIndex,
        long long particleIndex){

    particle->alreadyMoved = 1;
    double acceleration_x = Fx/particle->m;
    double acceleration_y = Fy/particle->m;
    particle->vx += acceleration_x;
    particle->vy += acceleration_y;
    double x= particle->x + particle->vx + acceleration_x/2;
    particle->y = particle->y + particle->vy + acceleration_y/2;
    x = fmod(x, MAX_COORDINATES_VALUE);
    particle->x = x < 0 ? x + MAX_COORDINATES_VALUE : x;
    //if(cellIndex >= ncside * NUMBER_OF_GHOST_ROWS && cellIndex < (ncside * NUMBER_OF_GHOST_ROWS + ncside) && (particle->y < oldy)) printf("A %d PCI %d | %0.2f %0.2f | oLDY %0.2f\n", processId, particle->creationIndex, particle->x, particle->y, oldy);
    //Since the particles may come from another process
    return move_particle(particle, cells, ncside, cellIndex, particleIndex, processId);
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
void get_cell(long long unbounded_row, long long unbounded_column, cell *cells, cell * return_cell, long ncside){
    long cellIndex;
    if (unbounded_column >= 0 && unbounded_column < ncside){
        cellIndex = unbounded_row * ncside + unbounded_column;
    } else {
        long bounded_column = (unbounded_column + ncside) % ncside;
        cellIndex = unbounded_row * ncside + bounded_column;
    }

    return_cell->x = cells[cellIndex].x;
    return_cell->y = cells[cellIndex].y;
    return_cell->m = cells[cellIndex].m;
    return_cell->nParticles = cells[cellIndex].nParticles;
    return_cell->particles = cells[cellIndex].particles;
    return_cell->allocatedSpace = cells[cellIndex].allocatedSpace;

    if (unbounded_column < 0)
        return_cell->x -= MAX_COORDINATES_VALUE;
    else if (unbounded_column >= ncside)
        return_cell->x += MAX_COORDINATES_VALUE;
}

/**§
 * Function that computes the force being applied to all particles and updates their cell after being moved
 *
 * @param particles         list with the particles
 * @param particlesLength  number of particles in the list
 * @param cells             matrix that represents the cells
 * @param ncside            number of cells in each side
 * @param cell_dimension    dimension of each cell
 */
void compute_force_and_update_particles(cell *cells, long ncside, int processId){
    for(long cellIndex = NUMBER_OF_GHOST_ROWS * ncside; cellIndex < TOTAL_ELEMENTS - (ncside * NUMBER_OF_GHOST_ROWS); cellIndex++) {
        cell * currCell = &cells[cellIndex];

        long cellY = ceil(cellIndex / ncside);
        long cellX = cellIndex % ncside;

        //printf("Evaluating %d %d,%d in %d\n", cellIndex, cellX, cellY, processId);
        //Obtain the neighbours of the particles of this cell
        cell neighboursList[9];
        int neighboursListIndex = 0;
        for (int row = -1; row < 2; row++) {
            for (int column = -1; column < 2; column++) {
                get_cell(row + cellY, column + cellX, cells, &neighboursList[neighboursListIndex++], ncside);
            }
        }

        for(long long particleIndex = 0; particleIndex < currCell->nParticles; particleIndex++) {
            particle_t * currentParticle = &(currCell->particles[particleIndex]);
            if (currentParticle->alreadyMoved == 1) continue;

            //resultant force in X and Y
            double Fx = 0;
            double Fy = 0;
            for (neighboursListIndex = 0; neighboursListIndex < 9; neighboursListIndex++) {
                //if that cell doesn't have any particle, no mass will exist, so we ignore
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
            particleIndex -= update_particle_position(currentParticle, Fx, Fy, cells, ncside, processId, cellIndex, particleIndex);
        }
    }

    exchangeGhostRowsParticles(cells, ncside, processId);
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

    
    for(long cellIndex = ncside * NUMBER_OF_GHOST_ROWS;
        cellIndex < TOTAL_ELEMENTS - (ncside * NUMBER_OF_GHOST_ROWS);
        cellIndex++) {
        particle_t * currentParticleList = cells[cellIndex].particles;
         #pragma omp declare reduction \
        (cntmassrdc:cell:omp_out=reduceCell(omp_out,omp_in)) \
        initializer(omp_priv=initCell(omp_priv))

        #pragma omp parallel for reduction(cntmassrdc:overallCenterMass)
        for(long long particleIndex = 0; particleIndex < cells[cellIndex].nParticles; particleIndex++) {
            particle_t currentParticle = currentParticleList[particleIndex];
            overallCenterMass.x += currentParticle.x * currentParticle.m;
            overallCenterMass.y += currentParticle.y * currentParticle.m;
            overallCenterMass.m += currentParticle.m;
        }
    }

    cell outOverallCenterMass = {.x=0, .y=0, .m=0};

    MPI_Reduce(&overallCenterMass, &outOverallCenterMass, 1, cellMPIType, reduceOverallCellOp, 0, MPI_COMM_WORLD);

    outOverallCenterMass.x /= outOverallCenterMass.m;
    outOverallCenterMass.y /= outOverallCenterMass.m;

    if(processId == 0)
        printf("%0.2f %0.2f \n", outOverallCenterMass.x, outOverallCenterMass.y);
}

void mapCellToMPI(MPI_Datatype * newType){
    int blocklens[] = {3 /*doubles*/,3 /*long long*/};
    MPI_Aint extent;
    MPI_Type_extent(MPI_DOUBLE, &extent);
    MPI_Aint indices[] = {0, 3 * extent /* we have 3 doubles */};
    MPI_Datatype oldTypes[] = {MPI_DOUBLE, MPI_LONG_LONG};
    MPI_Type_struct(2, blocklens, indices, oldTypes, newType);
    MPI_Type_commit(newType);
}

void mapParticleToMPI(MPI_Datatype * newType){
    int blocklens[] = {5 /*doubles*/,3 /*long*/};
    MPI_Aint extentDouble, extentLong;
    MPI_Type_extent(MPI_DOUBLE, &extentDouble);
    MPI_Aint indices[] = {0, 5 * extentDouble /* we have 5 doubles */};
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
    long ncside = strtol(args[2 + factor], NULL, 10);
    long long n_part = strtol(args[3 + factor], NULL, 10);
    long iterations = strtol(args[4 + factor], NULL, 10);

    particle_t * particles = calloc(n_part, sizeof(particle_t));

    int rank;
    MPI_Init( &args_length, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &NUMBER_OF_PROCESSES);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    TOTAL_ELEMENTS = NUMBER_OF_ELEMENTS(rank, NUMBER_OF_PROCESSES, ncside) + (ncside * NUMBER_OF_GHOST_ROWS * 2);

    printf("PID %d NOE %d\n",rank, TOTAL_ELEMENTS);
    //Type & Operations Setup
    mapCellToMPI(&cellMPIType);
    mapParticleToMPI(&particleMPIType);
    MPI_Op_create((MPI_User_function *) reduceOverallCellsMatrix, 1, &reduceOverallCellOp);

    /*
     * CellMatrix organization:
     *    ----------------------------------------------------------
     *   |  Top Ghost Row  |  Process's cells  |  Bottom Ghost Row  |
     *    ----------------------------------------------------------
     */
    cell * cellMatrix = (cell*) calloc(TOTAL_ELEMENTS, sizeof(cell));
    init_particles(seed, ncside, n_part, particles);


//    MPI_Bcast(particles, n_part, particleMPIType, 0, MPI_COMM_WORLD);

    // Populates the cell with their particles
    for(long long particleIndex = 0; particleIndex < n_part; particleIndex++){
        update_particle(&particles[particleIndex], &cellMatrix[ncside * NUMBER_OF_GHOST_ROWS], ncside, rank);
    }



    for(int i = 0; i < iterations; i++){
        //To simplify the index treatment (to start with 0) the cells matrix that goes through parameter omits the top ghost row
        printf("%d ITERATION %d\n", rank,i);
        compute_cell_center_mass(&cellMatrix[ncside * NUMBER_OF_GHOST_ROWS], ncside, rank);
        
        compute_force_and_update_particles(cellMatrix, ncside, rank);

        clean_cells(cellMatrix, ncside, rank);

        MPI_Barrier(MPI_COMM_WORLD);
    }


    int found = 0;
    for (long cellIndex = ncside * NUMBER_OF_GHOST_ROWS; !found && cellIndex < TOTAL_ELEMENTS - (ncside * NUMBER_OF_GHOST_ROWS); cellIndex++){
        for(long long particleIndex = 0;!found && particleIndex < cellMatrix[cellIndex].nParticles; particleIndex++){
            if(cellMatrix[cellIndex].particles[particleIndex].creationIndex == 0){
                printf("%0.2f %0.2f \n", cellMatrix[cellIndex].particles[particleIndex].x, cellMatrix[cellIndex].particles[particleIndex].y);
                found = 1;
            }
        }
    }

    compute_overall_center_mass(cellMatrix, ncside, rank);

    MPI_Finalize();

    free(cellMatrix);
    free(particles);
}
