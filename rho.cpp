/******************************************************************************
* FILE: master.c
* DESCRIPTION:
*   Processor Farm method for Mandelbrot set
*   The row buffer contains the pixel information with the last index containing
*   the row the worker just processed.
* AUTHOR: Sanjay Agravat
******************************************************************************/
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "h5_read.h"
#include "h5_write.h"

#define WORKTAG 1
#define DIETAG 2

#define BUFFER_SIZE 383

#define LINE_SIZE 32
#define NUM_ROWS 20004
#define NUM_COLS 78996
//#define NUM_ROWS 2
//#define NUM_COLS 3
#define WORK_SIZE NUM_ROWS*NUM_COLS


double mean(double *v, int n) {

    double sum = 0;
    //printf("mean\n");
    for (int i = 0; i < n; i++) {
        sum += v[i];
        //printf("%f,%d, %f\n", v[i], n, sum);

    }

    return sum/n;
}

double stddev(double *v, double avg, int n) {

    double sum = 0;
    for (int i = 0; i < n; i++) 
        sum += (v[i] - avg)*(v[i] - avg);

    return sqrt(sum/(n-1));
}

void sub(double *v, double avg, int n, double *result) {

    //double * result = (double*) calloc(n, sizeof(double));
    //double *out = NULL;
    //double result[n];

    for (int i = 0; i < n; i++) {
        result[i] = v[i] - avg;
        //printf("%f,%f, %f\n", v[i], avg, result[i]);
    }

    //out = result;
    //return out;
    //return result;
}

void prod(double *x, double *y, int n, double *result) {

    //double * result = (double*) calloc(n, sizeof(double));
    //double result[n];
    //double *out = NULL;

    for (int i = 0; i < n; i++) {
        result[i] = x[i]*y[i];
        //printf("%f,%f, %f\n", x[i], y[i], result[i]);
    }

    //out = result;
    //return out;
    //return result;
}

double sum(double *v, int n) {

    double result = 0;

    for (int i = 0; i < n; i++) {
        result += v[i];
    }

    return result;

}

double pvalue(double r, int n) {

    return r * sqrt( (n-2)/( 1 - (r*r)));

}

void getrowvec(double **m, int row, int nc, double *result) {

    //double *result = (double*) calloc(nc+1, sizeof(double));
    //double result[nc+1];
    //double *out = NULL;

    int i;
    for (int i = 0; i < nc; i++) {
        result[i] = m[row][i];
    }
    result[nc] = row;

    //out = result;
    //return out;
    //return result;
}

void getcolvec(double **m, int col, int nr, double *result) {

    //double *result = (double*) calloc(nr+1, sizeof(double));
    //double result[nr+1];
    //double *out = NULL;

    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }
    result[nr] = col;

    //out = result;
    //return out;
    //return result;
}



void master(void);
void worker(void);
int get_next_row(void);
int get_next_col(void);
void do_work(double *x, double *y, double *result);

int current_row, current_col, work_completed;
double starttime, endtime;
double *xcentered, *ycentered, *prod_result;

enum {
    tag_work,
    tag_free
};
int
main(int argc, char **argv)
{
  starttime = MPI_Wtime();
  int myrank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  xcentered   = (double*)calloc(BUFFER_SIZE,sizeof(double));;
  ycentered   = (double*)calloc(BUFFER_SIZE,sizeof(double));;
  prod_result = (double*)calloc(BUFFER_SIZE,sizeof(double));;

  //printf("starting up thread %d\n", myrank);
  if (myrank == 0) {
    master();
  } else {
    worker();
  }

  /* Shut down MPI */

  MPI_Finalize();
  return 0;
}


void
master(void)
{
  double *result = (double*)calloc(3,sizeof(double));
  double **RHO = NULL;
  double *x = (double*) calloc(BUFFER_SIZE+1, sizeof(double));
  double *y = (double*) calloc(BUFFER_SIZE+1, sizeof(double));

  int ntasks, rank;
  int row, col;

  MPI_Status status;
  current_row = 0, current_col = 0, work_completed = 0;

  /* Find out how many processes there are in the default
     communicator */
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

        FILE *f;
        char line[LINE_SIZE];
        int *probe_ids = NULL;
        int num_lines = 0;
        int i;
        f = fopen("probe_ids.txt", "r");
        while(fgets(line, LINE_SIZE, f)) {
                num_lines++;
                probe_ids = (int*)realloc(probe_ids, sizeof(int)*num_lines);
                probe_ids[num_lines-1] = atoi(strdup(line));
        }
        fclose(f);
        //MPI_File log_file;
        //MPI_File_open(MPI_COMM_WORLD, "/nethome/sagrava/CopyNumberVariation/checkpoint.txt",
                //MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND,
                //MPI_INFO_NULL, &log_file);

        //MPI_File_set_view(log_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR,
                                //MPI_CHAR, "external32", MPI_INFO_NULL);

        unsigned long x_nr, x_nc, y_nr, y_nc;

        double **X = h5_read("ge_cnv.h5", 1, "/X",         &x_nr, &x_nc);
        double **Y = h5_read("filtered_probes.h5", 1, "/FilteredProbes", &y_nr, &y_nc);

        printf("loaded X, num rows = %d, num cols = %d\n", x_nr, x_nc);
        printf("loaded Y, num rows = %d, num cols = %d\n", y_nr, y_nc);

        unsigned long total_mem = (x_nr * y_nc);
        RHO = (double**) malloc(x_nr*sizeof(double*));
        if (RHO == NULL) {
            printf("error allcoating array\n");
            exit(EXIT_FAILURE);
        }
        RHO[0] = (double*)malloc( total_mem*sizeof(double) );
        if (RHO[0] == NULL) {
            printf("error allcoating array\n");
            exit(EXIT_FAILURE);
        }
        for (unsigned long i=1; i < x_nr; i++) RHO[i] = RHO[0]+i*y_nc;

  /* Seed the workers; send one unit of work to each worker. */
  for (rank = 1; rank < ntasks; ++rank) {

    /* Find the next item of work to do */
    row = get_next_row();
    col = get_next_col();
    //printf("row = %d, col = %d\n", row, col);
    getrowvec(X, row, x_nc, x);
    //if (rank == 1)
        //for (int i = 0; i < BUFFER_SIZE; i++)
                //printf("x[%d] = %f\n", i, x[i]);
    getcolvec(Y, col, y_nr, y);
    //MPI_Barrier( MPI_COMM_WORLD ) ;
   
    //printf("************** got vecs %f, %f\n", x[BUFFER_SIZE], y[BUFFER_SIZE]);
    //printf("[master] got next work item %d\n", work);

    /* Send it to each rank */
    MPI_Send(x,                /* message buffer */
             BUFFER_SIZE+1,       /* one data item */
             MPI_DOUBLE,        /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    MPI_Send(y,                /* message buffer */
             BUFFER_SIZE+1,       /* one data item */
             MPI_DOUBLE,        /* data item is an integer */
             rank,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

     free(x);
     free(y);
  }

  /* Loop over getting new work requests until there is no more work
     to be done */
  //row = get_next_row();
  //col = get_next_col();
  //while (work != NULL) {
  while (work_completed < WORK_SIZE) {

    /* Receive results from a worker */
    MPI_Recv(result,           /* message buffer */
             3,                 /* one data item */
             MPI_DOUBLE,        /* of type double real */
             MPI_ANY_SOURCE,    /* receive from any sender */
             MPI_ANY_TAG,       /* any type of message */
             MPI_COMM_WORLD,    /* default communicator */
             &status);          /* info about the received message */

            int ix = (int)result[1];
            int iy = (int)result[2];
            RHO[ix][iy] = result[0];
        if (work_completed % 10000 == 0) {
            printf("work = %d, result[%d,%d] from %d = %f\n", work_completed,
                (int)result[1], (int)result[2], status.MPI_SOURCE, result[0]);
            //MPI_File_write_shared(log_file, log_str, strlen(log_str),
                                        //MPI_CHAR, MPI_STATUS_IGNORE);

        }

    /* Send the worker a new work unit */

    /* Get the next unit of work to be done */
    //printf("[master] send row %d to rank %d\n", work, rank);
    row = get_next_row();
    col = get_next_col();
    //printf("row = %d, col = %d\n", row, col);
    getrowvec(X, row, x_nc, x);
    getcolvec(Y, col, y_nr, y);
    //MPI_Barrier( MPI_COMM_WORLD ) ;
    MPI_Send(x,                /* message buffer */
             BUFFER_SIZE+1,       /* one data item */
             MPI_DOUBLE,        /* data item is an integer */
             status.MPI_SOURCE, /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */

    MPI_Send(y,                /* message buffer */
             BUFFER_SIZE+1,       /* one data item */
             MPI_DOUBLE,        /* data item is an integer */
             status.MPI_SOURCE,              /* destination process rank */
             WORKTAG,           /* user chosen message tag */
             MPI_COMM_WORLD);   /* default communicator */
    //printf("[master] got next work item%d\n", work);
  }

  /* There's no more work to be done, so receive all the outstanding
     results from the workers. */
  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Recv(result, 3, MPI_DOUBLE, rank,
             MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            int ix = (int)result[1];
            int iy = (int)result[2];
            RHO[ix][iy] = result[0];

        if (work_completed % 10000 == 0) 
            printf("work = %d, result[%d,%d] from %d = %f\n", work_completed,
                (int)result[1], (int)result[2], status.MPI_SOURCE, result[0]);
        //printf("result from %d = %f\n", status.MPI_SOURCE, result);
  }

  /* Tell all the workers to exit by sending an empty message with the
     DIETAG. */
  for (rank = 1; rank < ntasks; ++rank) {
    //printf("[master] send die tag to %d\n", rank);
    MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  }

  endtime   = MPI_Wtime();
  free(x);
  free(y);
  free(X[0]);
  free(X);
  free(Y[0]);
  free(Y);
  h5_write(RHO, NUM_ROWS, NUM_COLS, "rho.h5", "/rho");
  free(RHO[0]);
  free(RHO);
  fclose(f);
  printf("work completed = %d\n", work_completed);
  printf("%d - %f\n",ntasks, endtime-starttime);
}


void 
worker(void)
{
  double *x = (double*)calloc(BUFFER_SIZE+1,sizeof(double));
  double *y = (double*)calloc(BUFFER_SIZE+1,sizeof(double));
  double *work_out= (double*)calloc(3,sizeof(double));;

  MPI_Status status;
  int rank;

  memset(x, 0, BUFFER_SIZE);
  while (1) {

    /* Receive a message from the master */
    MPI_Recv(x, BUFFER_SIZE+1, MPI_DOUBLE, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    if (status.MPI_TAG == DIETAG) {
      return;
    }
    MPI_Recv(y, BUFFER_SIZE+1, MPI_DOUBLE, 0, MPI_ANY_TAG,
             MPI_COMM_WORLD, &status);

    do_work(x, y, work_out);

    /* Send the result back */
    MPI_Send(work_out, 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  free(x);
  free(y);
  free(work_out);
}


int
get_next_row(void)
{
    int row = current_row;
    current_row++;
    work_completed++;
    //if (work_completed % NUM_ROWS == 0)
    //printf("current_row = %d\n", current_row);
    if (row >= NUM_ROWS-1) {
        current_row = 0;
    }
    return row;
}

int
get_next_col(void)
{
    int col = current_col;
    if (current_row == 0)
        current_col++;
    return col;
}

void
do_work(double *x, double *y, double *work_out)
{

    int x_nc = BUFFER_SIZE;
    int y_nr = BUFFER_SIZE;
    double avgx         = mean(x, x_nc);
    sub(x, avgx, x_nc, xcentered);

    double avgy         = mean(y, y_nr);
    sub(y, avgy, y_nr, ycentered);

    prod(xcentered, ycentered, x_nc, prod_result);
    double sum_prod = sum(prod_result, x_nc);

    double stdX = stddev(x, avgx, x_nc);
    double stdY = stddev(y, avgy, y_nr);
    double rho  = sum_prod/((x_nc-1)*(stdX*stdY));

    work_out[0] = rho;
    work_out[1] = x[BUFFER_SIZE];
    work_out[2] = y[BUFFER_SIZE];

}
