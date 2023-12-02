#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <omp.h>
#include <sys/time.h>
#include <time.h>

/*
 *	(x,y) = (1,1) for Glider.
 *	(x,y) = (10,30) for R-pentomino
 */

#define GLIDER_X 1
#define GLIDER_Y 1
#define R_PENTOMINO_X 10
#define R_PENTOMINO_Y 30

#define GENERATIONS 2000
#define GRID_SIZE 2048

/*
 *	Programacao Concorrente e Distribuida 2023 - Rainbow Game of Life
 *	Autor: Gustavo dos Santos Ferreira 140731
 */

float time_diff(struct timeval *start, struct timeval *end)
{
    return (end->tv_sec - start->tv_sec) + 1e-6*(end->tv_usec - start->tv_usec);
}

void initialShape(float*** grid) {
	/* GLIDER */
	(*grid)[GLIDER_X][GLIDER_Y+1] = 1.0;
	(*grid)[GLIDER_X+1][GLIDER_Y+2] = 1.0;
	(*grid)[GLIDER_X+2][GLIDER_Y] = 1.0;
	(*grid)[GLIDER_X+2][GLIDER_Y+1] = 1.0;
	(*grid)[GLIDER_X+2][GLIDER_Y+2] = 1.0;
	/* R-PENTOMINO */
	(*grid)[R_PENTOMINO_X][R_PENTOMINO_Y+1] = 1.0;
	(*grid)[R_PENTOMINO_X][R_PENTOMINO_Y+2] = 1.0;
	(*grid)[R_PENTOMINO_X+1][R_PENTOMINO_Y] = 1.0;
	(*grid)[R_PENTOMINO_X+1][R_PENTOMINO_Y+1] = 1.0;
	(*grid)[R_PENTOMINO_X+2][R_PENTOMINO_Y+1] = 1.0;
}

float calcMeanNeighbours(float*** grid, float** upperBorder, float** lowerBorder, int x, int y, int lineStart, int lineCount, int world_rank, int world_size) {
	int i;
	int n;
	float val;
	float posVal;
	float temp;
	int getUpper = 0;
	int getLower = 0;
	int proc;
	n = 0;
	val = 0;
	temp = (*grid)[x][y];
	for (i = 0; i < 9; i++) {
		int posX, posY;
		getUpper = 0;
		getLower = 0;
		if (i == 4) continue;
		posX = x - 1 + (i%3);
		posY = y - 1 + (i/3);
		// Infinite grid
		if (world_size > 1) {
			if (posX < 0) getUpper = 1;
			else if (posX > lineCount - 1) getLower = 1;
		}
		else {
			if (posX < 0) posX = GRID_SIZE - 1;
			else if (posX > GRID_SIZE - 1) posX = 0;
		}
		if (posY < 0) posY = GRID_SIZE - 1;
		else if (posY > GRID_SIZE - 1) posY = 0;
		if ((!getUpper && !getLower) || world_size == 1) {
			if ((*grid)[posX][posY]) {
				n++;
				val += (*grid)[posX][posY];
			}
		}
		else if (getUpper) {
			if ((*upperBorder)[posY]) {
				n++;
				val += (*upperBorder)[posY];
				//printf("[%d: %d, %d]\n", world_rank, lineStart-1, posY);
			}
		}
		else if (getLower) {
			if ((*lowerBorder)[posY]) {
				n++;
				val += (*lowerBorder)[posY];
				//printf("[%d: %d, %d]\n", world_rank, lineStart+lineCount, posY);
			}
		}
	}
	if (n != 2 && n != 3) {
		val = 0;
	}
	else if (n == 2) {
		val = temp;
	}
	else if (n == 3) {
		val = val / 8;
		if (val == 0) {
			val = 0.001;
		}
	}
	return val;
}

int countLivingCells(float*** grid, int world_rank, int lineStart, int lineCount) {
	int i;
	int j;
	int count;
	count = 0;
	for (i = 0; i < lineCount; i++) {
		for (j = 0; j < GRID_SIZE; j++) {
			if ((*grid)[i][j]) {
				count++;
				//printf("[%d: %d, %d]\n", i+lineStart, j);
			}
		}
	}
	//printf("Cell count for rank %d:%d\n", world_rank, count);
	return count;
}

int exchangeWithUpper(float*** grid, float** upperBorder, int lineCount, int world_rank, int world_size) {
	float val;
	int proc;
	proc = world_rank - 1;
	if (proc < 0) proc = world_size - 1;
	for (int i = 0; i < GRID_SIZE; i++) {
		MPI_Recv(&val, 1, MPI_FLOAT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		(*upperBorder)[i] = val;
		val = (*grid)[0][i];
		MPI_Send(&val, 1, MPI_FLOAT, proc, 0, MPI_COMM_WORLD);
	}
	return 0;
}

int exchangeWithLower(float*** grid, float** lowerBorder, int lineCount, int world_rank, int world_size) {
	float val;
	int proc;
	proc = (world_rank + 1) % world_size;
	for (int i = 0; i < GRID_SIZE; i++) {
		val = (*grid)[lineCount-1][i];
		MPI_Send(&val, 1, MPI_FLOAT, proc, 0, MPI_COMM_WORLD);
		MPI_Recv(&val, 1, MPI_FLOAT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		(*lowerBorder)[i] = val;
	}
	return 0;
}

int main (int argc, char *argv[]) {
	int lineStart, lineCount;
	int nProc, thisProc, otherProc;
	int ierr;
	int world_rank;
	int world_size;
	int local_sum;
	int global_sum;
	
	struct timeval start;
	struct timeval end;
	int i, j, k;
	float** cells;
	float** nextCells;

	float* upperBorder;
	float* lowerBorder;

	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	lineCount = (GRID_SIZE / world_size);
	lineStart = world_rank * lineCount;
	/* Initialize the arrays */
	cells = (float**) malloc(lineCount * sizeof(float*));
	for (i = 0; i < lineCount; i++) {
		cells[i] = (float*) malloc(GRID_SIZE * sizeof(float));
		for (j = 0; j < GRID_SIZE; j++) {
			cells[i][j] = 0.0;
		}
	}
	if (world_rank == 0) {
		initialShape(&cells);
	}
	nextCells = (float**) malloc(lineCount * sizeof(float*));
	for (i = 0; i < lineCount; i++) {
		nextCells[i] = (float*) malloc(GRID_SIZE * sizeof(float));
	}

	upperBorder = (float*) malloc(GRID_SIZE * sizeof(float));
	lowerBorder = (float*) malloc(GRID_SIZE * sizeof(float));

	/* Execution goes here */
	gettimeofday(&start, NULL);

	for (int k = 0; k < GENERATIONS; k++) {
		//printf("Rank %d at gen %d\n", world_rank, k);
		// Even proc
		if (world_rank % 2 == 0) {
			if (world_size > 1) {
				exchangeWithUpper(&cells, &upperBorder, lineCount, world_rank, world_size);
				exchangeWithLower(&cells, &lowerBorder, lineCount, world_rank, world_size);
			}
			for (i = 0; i < lineCount; i++) {
				for (j = 0; j < GRID_SIZE; j++) {
					//printf("[%d: %d, %d]\n", world_rank, i + lineStart, j);
					nextCells[i][j] = calcMeanNeighbours(&cells, &upperBorder, &lowerBorder, i, j, lineStart, lineCount, world_rank, world_size);
				}
			}
		}
		// Odd proc
		else {
			if (world_size > 1) {
				exchangeWithLower(&cells, &lowerBorder, lineCount, world_rank, world_size);
				exchangeWithUpper(&cells, &upperBorder, lineCount, world_rank, world_size);
			}
			for (i = lineCount-1; i >= 0; i--) {
				for (j = 0; j < GRID_SIZE; j++) {
					//printf("[%d: %d, %d]\n", world_rank, i + lineStart, j);
					nextCells[i][j] = calcMeanNeighbours(&cells, &upperBorder, &lowerBorder, i, j, lineStart, lineCount, world_rank, world_size);
				}
			}
		}

		// Finalize
		for (i = 0; i < lineCount; i++) {
			for (j = 0; j < GRID_SIZE; j++) {
				cells[i][j] = nextCells[i][j];
			}
		}
		countLivingCells(&cells, world_rank, lineStart, lineCount);
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	/* Execution ends here */
	local_sum = countLivingCells(&cells, world_rank, lineStart, lineCount);
	MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (world_rank == 0) {
		printf("[FINAL CELL COUNT FOR RANK %d: %d\n] ", world_rank, global_sum);
	}
	gettimeofday(&end, NULL);
	printf("Time taken: %.6f seconds.\n", time_diff(&start, &end));

	/* Free the arrays */
	for (i = 0; i < lineCount; i++) {
		free(cells[i]);
	}
	free(cells);

	for (i = 0; i < lineCount; i++) {
		free(nextCells[i]);
	}
	free(nextCells);

	MPI_Barrier(MPI_COMM_WORLD);
    ierr = MPI_Finalize();

    return 0;
}
