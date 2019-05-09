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
        par[i].index = i;
    }
}

void rmvList( cell * cell, long oldParticleArrayIndex){
    if(oldParticleArrayIndex + 1 < cell->nParticles){
        memcpy(&cell->particles[oldParticleArrayIndex], &cell->particles[cell->nParticles - 1], sizeof(particle_t));
        cell->particles[oldParticleArrayIndex].arrayIndex = oldParticleArrayIndex;
    }

    cell->nParticles--;
    //cell->allocatedSpace--;
//    for(long long particleIndex = oldParticleArrayIndex; particleIndex < cell->nParticles; particleIndex++)
//        cell->particles[particleIndex].arrayIndex -= 1;


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
        cell->allocatedSpace *= 2;
        if(cell->allocatedSpace > 0){
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
//    printf("Particle m: %0.2f | x: %0.2f | %0.2f\n", particle->m, particle->x, particle->y);
    //for(int i = 0; i < cell->nParticles; i++)
     //   printf("")

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
    long globalCellIndex = newCellY * ncside + newCellX;

    if(!(globalCellIndex >= BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside) &&
            globalCellIndex <= BLOCK_HIGH(processId, NUMBER_OF_PROCESSES, ncside))) return;


    int localCellIndex = globalCellIndex - BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside);
//    printf("Update Particle on Start | PID %d | GI %d | LI %d | X: %0.2f | Y: %0.2f\n", processId, globalCellIndex, localCellIndex, particle->x, particle->y);
    addList( &cells[localCellIndex], particle );
}

/**
 * Update the matrix coordinates of each particle
 *
 * @param particle  A pointer to a particle
 * @param ncside    The number of cells on each side of the matrix of cells
 */
void move_particle(particle_t *particle, cell * cells, long ncside, long oldCellIndex, int processId, long ci){
        //Verifies if this cell belongs to this processor
    double sizeCell = MAX_COORDINATES_VALUE/ncside;
    int y = particle->y/sizeCell;
    int x = particle->x/sizeCell;
    long globalCellIndex =  y*ncside + x;

    long localCellIndex = globalCellIndex - BLOCK_LOW(processId, NUMBER_OF_PROCESSES, ncside) + ncside * NUMBER_OF_GHOST_ROWS;
    if(localCellIndex == oldCellIndex) return;

    long oldParticleArrayIndex = particle->arrayIndex;
    if((globalCellIndex >= BLOCK_LOW(processId,NUMBER_OF_PROCESSES, ncside)) &&
    (globalCellIndex <= BLOCK_HIGH( processId, NUMBER_OF_PROCESSES, ncside))) {
//        printf("ENTROU NO ADD | PID %d || lets add/remove x %0.2f y %0.2f on ci%d (%0.2f , %0.2f) gi %d li %d \n",processId, particle->x, particle->y, oldCellIndex, lc, hc, globalCellIndex, localCellIndex);
        addList(&cells[localCellIndex], particle);
    }
//    printf("%d || lets remove x %0.2f y %0.2f on ci%d (rci %d) gi %d li %d \n",processId, particle->x, particle->y, oldCellIndex, ci, globalCellIndex, localCellIndex);
    rmvList( &cells[oldCellIndex], oldParticleArrayIndex);
}

/**
 * Cleans all cells and resets all parameters to 0 so that the previous state doesn't affect the new one
 * @param cells     bidimensional array with the grid of cells
 * @param ncside    sides of the grid (how many rows the grid has)
 */
void clean_cells(cell * cells, long ncside, int processId){
    for (long cellIndex = 0; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside); cellIndex++){
        cells[cellIndex].x = cells[cellIndex].y =
        cells[cellIndex].m = 0;
        for (long particleIndex = 0; particleIndex < cells[cellIndex].nParticles; particleIndex++)
            cells[cellIndex].particles[particleIndex].alreadyMoved = 0;
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

    /*for (long cellIndex = 0; cellIndex < TOTAL_ELEMENTS; cellIndex++){
        //printf("CREATE GRID | PID: %d | Cell Index: %d | Number of Particles: %d\n", processId, cellIndex, cells[cellIndex].nParticles);
        cells[cellIndex].particles = (particle_t *)calloc(1, sizeof(particle_t));
        cells[cellIndex].allocatedSpace = 1;
    }*/

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
#define SEND_TOP_GHOST_ROW_TAG 1
#define SEND_BOT_GHOST_ROW_TAG 2
#define SEND_TOP_GHOST_PARTICLES_TAG 3
#define SEND_BOT_GHOST_PARTICLES_TAG 4

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
void exchangeGhostRows (cell * cells, long ncside, int senderProcessId) {

    if (NUMBER_OF_PROCESSES == 1) return; // If it is only one processor, there is no need for communication
 /*  for(int i = 0 ; i < ncside * NUMBER_OF_GHOST_ROWS; i ++){
        printf("PID %d  | TOP GHOST SENT BY  |  CI %d  | NP: %d \n", senderProcessId, i, cells[i].nParticles);
    }

    cell * pointer = &cells[NUMBER_OF_ELEMENTS(senderProcessId, NUMBER_OF_PROCESSES, ncside) -
                            (ncside * NUMBER_OF_GHOST_ROWS)];
    for(int i = 0 ; i < ncside * NUMBER_OF_GHOST_ROWS; i ++){
        printf("PID %d  | BOT GHOST SENT BY  |  CI %d  | NP: %d \n", senderProcessId, i, pointer[i].nParticles);
    }*/

    int topProcessId = (senderProcessId - 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;
    int bottomProcessId = (senderProcessId + 1 + NUMBER_OF_PROCESSES) % NUMBER_OF_PROCESSES;

    MPI_Status statuses[2];
    MPI_Request requests[2];

    cell * topGhostRow = cells - (ncside * NUMBER_OF_GHOST_ROWS);
    cell * bottomGhostRow = cells + NUMBER_OF_ELEMENTS(senderProcessId, NUMBER_OF_PROCESSES, ncside);
//    printf("sPID %d bPID %d NOE %d || p1 %x p2 %x\n", senderProcessId, bottomProcessId, NUMBER_OF_ELEMENTS(bottomProcessId, NUMBER_OF_PROCESSES, ncside), bottomGhostRow - (ncside * NUMBER_OF_GHOST_ROWS), pointer);

//    printf("PID %d  | BOT | NP: %d NOF %d\n", senderProcessId, bottomGhostRow->nParticles, NUMBER_OF_ELEMENTS(bottomProcessId, NUMBER_OF_PROCESSES, ncside));

    MPI_Isend(cells, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, topProcessId,
             SEND_TOP_GHOST_ROW_TAG, MPI_COMM_WORLD, &requests[0]);

    MPI_Isend(bottomGhostRow - (ncside * NUMBER_OF_GHOST_ROWS),ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, bottomProcessId,
             SEND_BOT_GHOST_ROW_TAG, MPI_COMM_WORLD, &requests[1]);

    MPI_Recv(bottomGhostRow, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, bottomProcessId,
             SEND_TOP_GHOST_ROW_TAG, MPI_COMM_WORLD, &statuses[1]);
    MPI_Recv(topGhostRow, ncside * NUMBER_OF_GHOST_ROWS, cellMPIType, topProcessId,
              SEND_BOT_GHOST_ROW_TAG, MPI_COMM_WORLD, &statuses[0]);



    MPI_Waitall(2, requests, statuses);
    for (int statusesIndex = 0; statusesIndex < 2; statusesIndex++)
    printf("Sender Id: %d | Status index: %d | Cancelled: %d | Count: %d\n",
           senderProcessId, statusesIndex, statuses[statusesIndex]._cancelled, statuses[statusesIndex]._ucount);


/*    printf("---------------- after wait all--------\n");
   for(int i = 0 ; i < ncside * NUMBER_OF_GHOST_ROWS; i ++){
        printf("PID %d  | TOP GHOST RECV BY %d |  CI %d  | NP: %d \n", senderProcessId, topProcessId, i, topGhostRow[i].nParticles);
    }

    for(int i = 0 ; i < ncside * NUMBER_OF_GHOST_ROWS; i ++){
        printf("PID %d  | BOT GHOST RECV BY %d |  CI %d  | NP: %d \n", senderProcessId, bottomProcessId, i, bottomGhostRow[i].nParticles);
    }*/

    //count all particles that exists in frontier rows
    long long countTopParticlesToSend, countDownParticlesToSend, countTopParticlesToReceive, countDownParticlesToReceive;
    countTopParticlesToSend = countDownParticlesToSend = countTopParticlesToReceive = countDownParticlesToReceive = 0;
    cell * topParticlesToSend = cells;
    cell * downParticlesToSend = bottomGhostRow - ncside;
    cell * topParticlesToReceive = cells - ncside;
    cell * downParticlesToReceive = bottomGhostRow;
    for(long cellIndex = 0; cellIndex < ncside; cellIndex++){
        countTopParticlesToSend += topParticlesToSend[cellIndex].nParticles;
        countDownParticlesToSend += downParticlesToSend[cellIndex].nParticles;
        countTopParticlesToReceive += topParticlesToReceive[cellIndex].nParticles;
        countDownParticlesToReceive += downParticlesToReceive[cellIndex].nParticles;
    }
    printf("PID %d || CountTopPaRecv %d, CountTopPaSend %d || CountDownPaRecv %d CountDownPaSend %d\n", senderProcessId,
           countTopParticlesToReceive, countTopParticlesToSend,
           countDownParticlesToReceive, countDownParticlesToSend);

    for (long cellIndex = 0; cellIndex < ncside * NUMBER_OF_GHOST_ROWS; cellIndex++){
        if (senderProcessId == 0) topGhostRow[cellIndex].y -= MAX_COORDINATES_VALUE;
        else if (senderProcessId == NUMBER_OF_PROCESSES - 1) bottomGhostRow[cellIndex].y += MAX_COORDINATES_VALUE;
    }
//    printf("here1 by %d\n", senderProcessId);
    particle_t *topParticlesToSendBuffer = calloc(countTopParticlesToSend, sizeof(particle_t));//[countTopParticlesToSend];
    particle_t *downParticlesToSendBuffer = calloc(countDownParticlesToSend, sizeof(particle_t));//[countDownParticlesToSend];
    long long accumulatorTopParticlesToSend = 0;
    long long accumulatorDownParticlesToSend = 0;

//    printf("here2 by %d\n", senderProcessId);
    for(int cellIndex = 0; cellIndex < ncside; cellIndex++) {
//        printf("PID %d | Started top memcpy - copying %d particles\n", senderProcessId, topParticlesToSend[cellIndex].nParticles);
        memcpy(topParticlesToSendBuffer + accumulatorTopParticlesToSend,
               topParticlesToSend[cellIndex].particles, topParticlesToSend[cellIndex].nParticles * sizeof(particle_t));
        accumulatorTopParticlesToSend += topParticlesToSend[cellIndex].nParticles;

//        printf("PID %d | Started bot memcpy - copying %d particles\n", senderProcessId, downParticlesToSend[cellIndex].nParticles);
        memcpy(downParticlesToSendBuffer + accumulatorDownParticlesToSend,
               downParticlesToSend[cellIndex].particles,
               downParticlesToSend[cellIndex].nParticles * sizeof(particle_t));
        accumulatorDownParticlesToSend += downParticlesToSend[cellIndex].nParticles;

    }

//    printf("here4 by %d\n", senderProcessId);
    particle_t *topParticlesToReceiveBuffer = calloc(countTopParticlesToReceive, sizeof(particle_t));//[countTopParticlesToReceive];
    particle_t *downParticlesToReceiveBuffer = calloc(countDownParticlesToReceive, sizeof(particle_t)); //[countDownParticlesToReceive];
//    printf("here4 sID %d tID %d count %d\n", senderProcessId, topProcessId, countTopParticlesToSend);




    printf("pid %d after 1st wait toppid %d\n", senderProcessId, topProcessId);
    MPI_Irecv(topParticlesToReceiveBuffer, countTopParticlesToReceive, particleMPIType, topProcessId,
              SEND_BOT_GHOST_PARTICLES_TAG, MPI_COMM_WORLD, &requests[1]);

    printf("pid %d after 2nd wait botpid %d\n", senderProcessId, bottomProcessId);
    MPI_Irecv(downParticlesToReceiveBuffer, countDownParticlesToReceive, particleMPIType, bottomProcessId,
             SEND_TOP_GHOST_PARTICLES_TAG, MPI_COMM_WORLD, &requests[0]);
    printf("pid %d after Receives Particles wait \n", senderProcessId);

    printf("PId %d 1st send to top %d\n",senderProcessId, topProcessId);
    MPI_Send(topParticlesToSendBuffer, countTopParticlesToSend, particleMPIType, topProcessId,
             SEND_TOP_GHOST_PARTICLES_TAG, MPI_COMM_WORLD);//, &requests[0]);

    printf("PId %d 2nd send to bot %d\n",senderProcessId, bottomProcessId);
    MPI_Send(downParticlesToSendBuffer, countDownParticlesToSend, particleMPIType, bottomProcessId,
             SEND_BOT_GHOST_PARTICLES_TAG, MPI_COMM_WORLD);//, &requests[1]);

    printf("here4.1 sID %d bID %d count %d\n", senderProcessId, bottomProcessId, countDownParticlesToSend);

    MPI_Waitall(2, requests, statuses);
    printf("here6 by %d\n", senderProcessId);


    /*
    particle_t * aux = &topParticlesToSendBuffer;
    for(int i = 0; i < countTopParticlesToSend; i++){
        printf("ENVIAR UP | PID %d | X:%0.2f Y:%0.2f\n", senderProcessId, aux[i].x, aux[i].y);
    }

    aux = &downParticlesToReceiveBuffer;
    for(int i = 0; i < countDownParticlesToReceive; i++){
        printf("RECEBER DOWN | PID %d | X:%0.2f Y:%0.2f\n", senderProcessId, aux[i].x, aux[i].y);
    }

    aux = &downParticlesToSendBuffer;
    for(int i = 0; i < countDownParticlesToSend; i++){
        printf("ENVIAR DOWN | PID %d | X:%0.2f Y:%0.2f\n", senderProcessId, aux[i].x, aux[i].y);
    }

    aux = &topParticlesToReceiveBuffer;
    for(int i = 0; i < countTopParticlesToReceive; i++){
        printf("RECEBER UP | PID %d | X:%0.2f Y:%0.2f\n", senderProcessId, aux[i].x, aux[i].y);
    }*/

    particle_t * pointerTopParticlesToReceiveBuffer = topParticlesToReceiveBuffer;
    particle_t * pointerDownParticlesToReceiveBuffer = downParticlesToReceiveBuffer;
    
    for(long cellIndex = 0; cellIndex < ncside; cellIndex++){
//        topGhostRow[cellIndex].allocatedSpace = 1;
//        topGhostRow[cellIndex].nParticles = 0;
//        topGhostRow[cellIndex].particles = (particle_t *)calloc(topGhostRow[cellIndex].allocatedSpace * sizeof(particle_t));


        topParticlesToReceive[cellIndex].particles = (particle_t *)malloc(topParticlesToReceive[cellIndex].allocatedSpace * sizeof(particle_t));
        memcpy(topParticlesToReceive[cellIndex].particles, pointerTopParticlesToReceiveBuffer, topParticlesToReceive[cellIndex].nParticles * sizeof(particle_t));
        pointerTopParticlesToReceiveBuffer += topParticlesToReceive[cellIndex].nParticles;

        for(long long particleIndex = 0;
            senderProcessId == 0 && particleIndex < topParticlesToReceive[cellIndex].nParticles; particleIndex++){
            topParticlesToReceive[cellIndex].particles[particleIndex].y -= MAX_COORDINATES_VALUE;
        }

//        bottomGhostRow[cellIndex+ncside].allocatedSpace = 1;
//        bottomGhostRow[cellIndex+ncside].nParticles = 0;
//        bottomGhostRow[cellIndex+ncside].particles = (particle_t *)calloc(bottomGhostRow[cellIndex+ncside].allocatedSpace, sizeof(particle_t));

        downParticlesToReceive[cellIndex].particles = (particle_t *)malloc(downParticlesToReceive[cellIndex].allocatedSpace * sizeof(particle_t));
        memcpy(downParticlesToReceive[cellIndex].particles, pointerDownParticlesToReceiveBuffer, downParticlesToReceive[cellIndex].nParticles * sizeof(particle_t));
        pointerDownParticlesToReceiveBuffer += downParticlesToReceive[cellIndex].nParticles;
        for(long long particleIndex = 0;
        (senderProcessId == NUMBER_OF_PROCESSES - 1) && particleIndex < downParticlesToReceive[cellIndex].nParticles; particleIndex++){
            downParticlesToReceive[cellIndex].particles[particleIndex].y += MAX_COORDINATES_VALUE;
        }
    }

    free(topParticlesToReceiveBuffer);
    free(downParticlesToReceiveBuffer);
    free(topParticlesToSendBuffer);
    free(downParticlesToSendBuffer);

    // Check the status for a possible error TODO use the count and cancelled to check errors
//    for (int statusesIndex = 0; statusesIndex < 2; statusesIndex++)
//        printf("Sender Id: %d | Status index: %d | Cancelled: %d | Count: %d\n",
//               senderProcessId, statusesIndex, statuses[statusesIndex]._cancelled, statuses[statusesIndex]._ucount);

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
//    printf("PID %d started center mass\n", processId);
    clean_cells(cells, ncside, processId);

    for (long cellIndex = 0; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside);
    cellIndex++){
        cell currCell = cells[cellIndex];
        if(currCell.nParticles > 0){ //avoid div by 0
            for(long long particleIndex = 0; particleIndex < currCell.nParticles; particleIndex++) {
                particle_t * currentParticle = &currCell.particles[particleIndex];
//                printf("PID %d CI %d PI %d X %0.2f Y %0.2f \n", processId, cellIndex, particleIndex, currentParticle->x, currentParticle->y);
                cells[cellIndex].x += currentParticle->x * currentParticle->m;
                cells[cellIndex].y += currentParticle->y * currentParticle->m;
                cells[cellIndex].m += currentParticle->m;
            }
            cells[cellIndex].x /= cells[cellIndex].m;
            cells[cellIndex].y /= cells[cellIndex].m;
        }
//        printf("PID %d CI %d X %0.2f\n", processId, cellIndex, cells[cellIndex].x);
    }

    exchangeGhostRows(cells, ncside, processId);

    for (long cellIndex = -ncside; cellIndex < NUMBER_OF_ELEMENTS(processId, NUMBER_OF_PROCESSES, ncside)+ncside;
         cellIndex++){
//        printf("PID %d CI %d X %0.2f\n", processId, cellIndex, cells[cellIndex].x);
    }
//    printf("PID %d center mass\n", processId);
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
void update_particle_position(particle_t * particle, double Fx, double Fy, cell * cells, long ncside, int processId, long cellIndex){

//    printf("%d B|| px %0.2f py %0.2f m %0.2f ci %d\n", processId, particle->x, particle->y, particle->m, cellIndex);
    particle->alreadyMoved = 1;
    double acceleration_x = Fx/particle->m;
    double acceleration_y = Fy/particle->m;
    particle->vx += acceleration_x;
    particle->vy += acceleration_y;
    double x= particle->x + particle->vx + acceleration_x/2;
    particle->y = particle->y + particle->vy + acceleration_y/2;
    x = fmod(x, MAX_COORDINATES_VALUE);
    particle->x = x < 0 ? x + MAX_COORDINATES_VALUE : x;

 //   printf("%d A|| px %0.2f py %0.2f ci %d\n", processId, particle->x, particle->y, cellIndex);
    //Since the particles may come from another process
     move_particle(particle, cells, ncside, cellIndex, processId, cellIndex);
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
void get_cell(int pid, long long unbounded_row, long long unbounded_column, cell *cells, cell * return_cell, long ncside){
    long cellIndex;
    if (unbounded_column >= 0 && unbounded_column < ncside){
        cellIndex = unbounded_row * ncside + unbounded_column;
    } else {
        long bounded_column = (unbounded_column + ncside) % ncside;
//        printf("PID %d || UR %d UC %d BC %d \n", pid, unbounded_row, unbounded_column, bounded_column);
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
//    printf("PID %d started moved TOTal elem %d\n", processId, TOTAL_ELEMENTS);
    for(long cellIndex = ncside;cellIndex < TOTAL_ELEMENTS - ncside; cellIndex++) {
//        printf("PID % started for move in %d\n", processId, cellIndex);
        cell * currCell = &cells[cellIndex];

        long cellX = ceil(cellIndex / ncside);
        long cellY = cellIndex % ncside;

        //Obtain the neighbours of the particles of this cell
        cell neighboursList[9];
        int neighboursListIndex = 0;
        for (int row = -1; row < 2; row++) {
            for (int column = -1; column < 2; column++) {
                get_cell(processId, row + cellY, column + cellX, cells, &neighboursList[neighboursListIndex++], ncside);
//                if(neighboursList[neighboursListIndex-1].nParticles > 0) neighboursList[neighboursListIndex-1].particles[0].x = 1;
            }
        }
//        for (neighboursListIndex = 0; neighboursListIndex < 9; neighboursListIndex++) {
//            printf("PID %d | NI %d NP %d\n", processId, neighboursListIndex, neighboursList[neighboursListIndex].nParticles);
//        }
//        printf("PID %d obtained list of neihbours\n", processId);
        for(long long particleIndex = 0; particleIndex < currCell->nParticles; particleIndex++) {
//            printf("CHECK NUMBER PART | PID %d | CELL INDEX: %d | NUMBER OF PART: %d\n", processId, cellIndex, currCell->nParticles);
            particle_t * currentParticle = &(currCell->particles[particleIndex]);
//            printf("CHECK2 NUMBER PART | PID %d | CELL INDEX: %d | NUMBER OF PART: %d\n", processId, cellIndex, currCell->nParticles);
//            printf("CHECK3 NUMBER PART | PID %d | CELL INDEX: %d | NUMBER OF PART: %d | AM %d\n", processId, cellIndex, currCell->nParticles,currCell->particles[particleIndex].alreadyMoved);
            if (currentParticle->alreadyMoved) continue;
//            printf("CHECK4 NUMBER PART | PID %d | CELL INDEX: %d | NUMBER OF PART: %d\n", processId, cellIndex, currCell->nParticles);
            //resultant force in X and Y
            double Fx = 0;
            double Fy = 0;
            for (neighboursListIndex = 0; neighboursListIndex < 9; neighboursListIndex++) {
                //if that cell doesn't have any particle, no central mass will exist, so we ignore
               // printf("PID %d  ||  NI %d  ||  NP %d\n", processId, neighboursListIndex,
               //         neighboursList[neighboursListIndex].nParticles);
                if (neighboursList[neighboursListIndex].m == 0)
                    continue;
//                if(currentParticle->index == 0)
//                    printf("PID %d CI %d PI %d NI%d(%0.2f, %0.2f, m = %0.2f) || PX %0.2f PY %0.2f\n",processId, cellIndex, particleIndex, neighboursListIndex,neighboursList[neighboursListIndex].x, neighboursList[neighboursListIndex].y, neighboursList[neighboursListIndex].m,
//                           currentParticle->x, currentParticle->y);

                //compute angle
                double delta_x = neighboursList[neighboursListIndex].x - currentParticle->x;
                double delta_y = neighboursList[neighboursListIndex].y - currentParticle->y;
                double vector_angle = atan2(delta_y, delta_x);

                //compute force
                double force = compute_magnitude_force(currentParticle, &neighboursList[neighboursListIndex]);
                Fx += force * cos(vector_angle);
                Fy += force * sin(vector_angle);
            }
//            printf("update particle PID %d ||| CI %d\n", processId, cellIndex);

            update_particle_position(currentParticle, Fx, Fy, cells, ncside, processId, cellIndex);
        }
        printf("PID %d CI %d\n", processId, cellIndex);
    }
    printf("PID %d moved\n", processId);
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
        for(long long particleIndex = 0; particleIndex < cells[cellIndex].nParticles; particleIndex++) {
            particle_t currentParticle = currentParticleList[particleIndex];
            overallCenterMass.x += currentParticle.x * currentParticle.m;
            overallCenterMass.y += currentParticle.y * currentParticle.m;
            overallCenterMass.m += currentParticle.m;
        }
    }

    cell outOverallCenterMass = {.x=0, .y=0, .m=0};
   // printf("Process id: %d  |  %0.2f %0.2f %0.2f \n", processId, overallCenterMass.x, overallCenterMass.y, overallCenterMass.m);

    MPI_Reduce(&overallCenterMass, &outOverallCenterMass, 1, cellMPIType, reduceOverallCellOp, 0, MPI_COMM_WORLD);


    if(processId == 0){
        outOverallCenterMass.x /= outOverallCenterMass.m;
        outOverallCenterMass.y /= outOverallCenterMass.m;
        printf("%0.2f %0.2f \n", outOverallCenterMass.x, outOverallCenterMass.y);
    }
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
    int blocklens[] = {5 /*doubles*/,4 /*long*/};
    MPI_Aint extentDouble, extentLong;
    MPI_Type_extent(MPI_DOUBLE, &extentDouble);
    MPI_Aint indices[] = {0, 5 * extentDouble /* we have 5 doubles */};
    MPI_Datatype oldTypes[] = {MPI_DOUBLE, MPI_LONG};
    MPI_Type_struct(2, blocklens, indices, oldTypes, newType);
    MPI_Type_commit(newType);
}

int main(int args_length, char* args[]) {
    /*
     * LIST OF TODO
     *  MALLOC space for a particle only if there is one (use allocatedSpace to verify the need)
     *
     * */
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

    printf("PID %d NP %d| LB %0.2f HB %0.2f\n", rank,NUMBER_OF_PROCESSES ,(BLOCK_LOW(rank, NUMBER_OF_PROCESSES, ncside)), (BLOCK_HIGH(rank, NUMBER_OF_PROCESSES, ncside)));
    TOTAL_ELEMENTS = NUMBER_OF_ELEMENTS(rank, NUMBER_OF_PROCESSES, ncside) + (ncside * NUMBER_OF_GHOST_ROWS * 2);

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
    cell * cellMatrix = (cell*) calloc(TOTAL_ELEMENTS, sizeof(cell));
    if(rank==0){
        init_particles(seed, ncside, n_part, particles);
    }


    MPI_Bcast(particles, n_part, particleMPIType, 0, MPI_COMM_WORLD);

    create_grid(particles, n_part, cellMatrix, ncside, rank);

//    for (long cellIndex = ncside; cellIndex < TOTAL_ELEMENTS - ncside; cellIndex++){
//        for(long long particleIndex = 0; particleIndex < cellMatrix[cellIndex].nParticles; particleIndex++){
//            printf("PID %d CI %d PI %d||| px %0.2f  py %0.2f m %0.2f\n", rank, cellIndex, particleIndex, cellMatrix[cellIndex].particles[particleIndex].x, cellMatrix[cellIndex].particles[particleIndex].y, cellMatrix[cellIndex].particles[particleIndex].m);
//        }
//    }
//    printf("-------------\n");
//    for(int i = 0; i < n_part; i++)
//        printf("Process id: %d | Particula %d | %0.2f %0.2f %0.2f \n", rank, i, particles[i].x, particles[i].y, particles[i].m);
    /*
    for(long particleIndex = BLOCK_LOW(me, NUMBER_OF_PROCESSES, n_part); 
            particleIndex <= BLOCK_HIGH(me, NUMBER_OF_PROCESSES, n_part); 
            particleIndex++){
        printf("Process id: %d | Particula %d | %0.2f %0.2f %0.2f \n", me, particleIndex, particles[particleIndex].x, particles[particleIndex].y, particles[particleIndex].m);
    }*/

    for(int i = 0; i < iterations; i++){
//        for (long cellIndex = ncside; cellIndex < TOTAL_ELEMENTS - ncside; cellIndex++){
//            for(long long particleIndex = 0; particleIndex < cellMatrix[cellIndex].nParticles; particleIndex++){
//                printf("PID %d CI %d PI %d PRI %d||| px %0.2f  py %0.2f\n", rank, cellIndex, particleIndex, cellMatrix[cellIndex].particles[particleIndex].index, cellMatrix[cellIndex].particles[particleIndex].x, cellMatrix[cellIndex].particles[particleIndex].y);
//            }
//        }

        //To simplify the index treatment (to start with 0) the cells matrix that goes through parameter omits the top ghost row
        compute_cell_center_mass(&cellMatrix[ncside * NUMBER_OF_GHOST_ROWS], ncside, rank);

        compute_force_and_update_particles(cellMatrix, ncside, rank);
    }

    int found = 0;
    for (long cellIndex = ncside * NUMBER_OF_GHOST_ROWS; !found && cellIndex < TOTAL_ELEMENTS - (ncside * NUMBER_OF_GHOST_ROWS); cellIndex++){
        for(long long particleIndex = 0;!found && particleIndex < cellMatrix[cellIndex].nParticles; particleIndex++){
            if(cellMatrix[cellIndex].particles[particleIndex].index == 0){
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
