 /* Gaussian elimination without pivoting.
 * Compile with "gcc gauss.c"
 */

/* ****** ADD YOUR CODE AT THE END OF THIS FILE. ******
 * You need not submit the provided code.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
//#include <omp.h>
#include <mpi.h>

/* Program Parameters */
#define MAXN 2000  /* Max value of N */
int N;  /* Matrix size */
int p = 4;	//number of processors created. Default = 4. Command line argument can modify it.
int id;
/* Matrices and vectors */
//volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];
//MPI does not like 'volatile' matrix
/* A * X = B, solve for X */

double A[MAXN][MAXN], B[MAXN], X[MAXN];

/* junk */
#define randm() 4|2[uid]&3

void backSubstitution();
/* Prototype */
void gauss_mpi();  /* The function you will provide.
		* It is this routine that is timed.
		* It is called only on the parent.
		*/

/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}

/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
	int seed = 0;  /* Random seed */
	char uid[32]; /*User name */

	/* Read command-line arguments */
	srand(time_seed());  /* Randomize */

	if (argc == 3) {
		seed = atoi(argv[2]);
		srand(seed);
		p=atoi(argv[3]);	//get number of processors
		printf("Random seed = %i\n", seed);
	}
	if (argc >= 2) {
		N = atoi(argv[1]);
		if (N < 1 || N > MAXN) {
			printf("N = %i is out of range.\n", N);
			exit(0);
		}
	}
	else {
		printf("Usage: %s <matrix_dimension> [random seed]\n",
				argv[0]);
		exit(0);
	}

	/* Print parameters */
	printf("\nMatrix dimension N = %i.\n", N);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {
  int row, col;

  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}

/* Print input matrices */
void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
	printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 100) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}

int main(int argc, char **argv) {

	//MPI implementation

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


	/* Timing variables */
	struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
	struct timezone tzdummy;
	clock_t etstart2, etstop2;  /* Elapsed times using times() */
	unsigned long long usecstart, usecstop;
	struct tms cputstart, cputstop;  /* CPU times for my processes */

	/* Process program parameters */
	parameters(argc, argv);

	//if rank = 0 do main
	if(myid == 0){
		/* Initialize A and B */
		initialize_inputs();

		/* Print input matrices */
		print_inputs();

		/* Start Clock */
		printf("\nStarting clock.\n");
		gettimeofday(&etstart, &tzdummy);
		etstart2 = times(&cputstart);
	}
	/* Gaussian Elimination */
	gauss_mpi();

	if(myid == 0){
		/* Stop Clock */
		backSubstitution();
		gettimeofday(&etstop, &tzdummy);
		etstop2 = times(&cputstop);
		printf("Stopped clock.\n");
		usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
		usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

		/* Display output */
		print_X();

		/* Display timing results */
		printf("\nElapsed time = %g ms.\n",
		 (float)(usecstop - usecstart)/(float)1000);

		printf("(CPU times are accurate to the nearest %g ms)\n",
		 1.0/(float)CLOCKS_PER_SEC * 1000.0);
		printf("My total CPU time for parent = %g ms.\n",
		 (float)( (cputstop.tms_utime + cputstop.tms_stime) -
			  (cputstart.tms_utime + cputstart.tms_stime) ) /
		 (float)CLOCKS_PER_SEC * 1000);
		printf("My system CPU time for parent = %g ms.\n",
		 (float)(cputstop.tms_stime - cputstart.tms_stime) /
		 (float)CLOCKS_PER_SEC * 1000);
		printf("My total CPU time for child processes = %g ms.\n",
		 (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
			  (cputstart.tms_cutime + cputstart.tms_cstime) ) /
		 (float)CLOCKS_PER_SEC * 1000);
		  /* Contrary to the man pages, this appears not to include the parent */
		printf("--------------------------------------------\n");
	}
	MPI_Finalize();	//end MPI
	exit(0);
}

/* ------------------ Above Was Provided --------------------- */

/****** You will replace this routine with your own parallel version *******/
/* Provided global variables are MAXN, N, A[][], B[], and X[],
 * defined in the beginning of this code.  X[] is initialized to zeros.
 */
void gauss_mpi() {
	int norm, row, col;  /* Normalization row, and zeroing
			* element row and col */
	float multiplier;
	printf("Computing in parallel in MPI.\n");

	double endtime;
	double starttime = 0;
	int proc;

	MPI_Request request;
	MPI_Status status;

	MPI_Barrier(MPI_COMM_WORLD);	//waiting for p processors

	//processor 0 starts timer
	if(myid == 0){
		starttime = MPI_Wtime();
	}

	/* Gaussian elimination */
	for (norm = 0; norm < N - 1; norm++) {
		// bcast A and B matrix from 0th processor to other processors
		MPI_Bcast(&A[norm][0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&B[norm], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		//0th rank job
		if(myid == 0){
			//sending  A and B matrix to other processors
			for(proc=1;proc<p;proc++){
				for(row = norm + 1 + proc; row < N; row+=p){
					MPI_Isend(&A[row], N, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
					MPI_Isend(&B[row], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
				}
			}

			//doing gauss elimination on processor 0
			for (row = norm + 1; row < N; row+=p) {
				multiplier = A[row][norm] / A[norm][norm];
				for (col = norm; col < N; col++) {
					A[row][col] -= A[norm][col] * multiplier;
				}
				B[row] -= B[norm] * multiplier;
			}

			//receive values from other processors
			for(proc = 1; proc<p ; proc++){
				for(row = norm + 1 + proc; row < N; row+=p){
					MPI_Recv(&A[row], N, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, &status);
					MPI_Recv(&B[row], 1, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, &status);
				}
			}


			//stop timer
			if(norm == N - 2){
				endtime = MPI_Wtime();
				printf("MPI elapsed time = %f\n", endtime - starttime);
			}
		}
		//other processors
		else{
			for(row = norm + 1 + myid; row < N; row+=p){
				//receive A and B matrix from processor 0
				MPI_Recv(&A[row], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&B[row], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
				multiplier = A[row][norm] / A[norm][norm];
				for (col = norm; col < N; col++) {
					A[row][col] -= A[norm][col] * multiplier;
				}
				B[row] -= B[norm] * multiplier;

				//send back values to processor 0
				MPI_Isend(&A[row], N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &status);
				MPI_Isend(&B[row], 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
				MPI_Wait(&request, &status);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);	//wait for p processors
	}
}

void backSubstitution(){
	int row, col;
	for (row = N - 1; row >= 0; row--) {
      X[row] = B[row];
      for (col = N-1; col > row; col--) {
        X[row] -= A[row][col] * X[col];
      }
      X[row] /= A[row][row];
    }
}
