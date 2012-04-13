/******************************************************************************
* FILE: master.c
* DESCRIPTION:
*   Processor Farm method for Mandelbrot set
*   The row buffer contains the pixel information with the last index containing
*   the row the worker just processed.
* AUTHOR: Sanjay Agravat
******************************************************************************/
#include <mpi.h>
#include <omp.h>
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
//#define NUM_COLS 52664
//#define NUM_ROWS 2
//#define NUM_COLS 3
#define WORK_SIZE NUM_ROWS*NUM_COLS


double mean(double *v, int n) {

    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += v[i];

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


    for (int i = 0; i < n; i++) {
        result[i] = v[i] - avg;
    }
}

void prod(double *x, double *y, int n, double *result) {

    for (int i = 0; i < n; i++) {
        result[i] = x[i]*y[i];
    }

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

double* getrowvec(double **m, int row, int nc) {

    double *result = (double*)calloc(BUFFER_SIZE+1,sizeof(double));;

    //int i;
    memcpy(result, m[row], sizeof(double) * nc);
    //for (int i = 0; i < nc; i++) {
    //    result[i] = m[row][i];
    //}
    result[nc] = row;
    return result;

}

double* getcolvec(double **m, int col, int nr) {

    double *result  = (double*)calloc(BUFFER_SIZE+1,sizeof(double));;
    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }
    result[nr] = col;

    return result;

}

void getrowvec(double **m, int row, int nc, double *result) {


    //int i;
    memcpy(result, m[row], sizeof(double) * nc);
    //for (int i = 0; i < nc; i++) {
    //    result[i] = m[row][i];
    //}
    result[nc] = row;

}

void getcolvec(double **m, int col, int nr, double *result) {

    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }
    result[nr] = col;

}



int
main(int argc, char **argv)
{
    double starttime, endtime;
    starttime = MPI_Wtime();
    int ntasks, myrank;
    char name[128];                      
    int namelen;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    MPI_Get_processor_name(name,&namelen);

    /* Get a few OpenMP parameters.                                               */
    int O_P  = omp_get_num_procs();          /* get number of OpenMP processors       */
    int O_T  = omp_get_num_threads();        /* get number of OpenMP threads          */
    int O_ID = omp_get_thread_num();         /* get OpenMP thread ID                  */
    printf("name:%s   M_ID:%d  O_ID:%d  O_P:%d  O_T:%d\n", name,myrank,O_ID,O_P,O_T);

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
    unsigned long x_nr, x_nc, y_nr, y_nc;

    double **X = h5_read("ge_cnv.h5", 1, "/X",         &x_nr, &x_nc);
    double **Y = h5_read("filtered_probes.h5", 1, "/FilteredProbes", &y_nr, &y_nc);

    printf("loaded X, num rows = %d, num cols = %d\n", x_nr, x_nc);
    printf("loaded Y, num rows = %d, num cols = %d\n", y_nr, y_nc);

    unsigned long total_mem = (x_nr * y_nc);
    double **RHO = (double**) malloc(x_nr*sizeof(double*));
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
    int gene, probe, tid, work_completed;
    work_completed = 0;
    double *x = (double*)calloc(BUFFER_SIZE+1,sizeof(double));;
    double *y = (double*)calloc(BUFFER_SIZE+1,sizeof(double));;

    double *xcentered   = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    double *ycentered   = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    double *prod_result = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    int BLOCK_SIZE = NUM_ROWS/ntasks;
    int offset = myrank*BLOCK_SIZE;

    printf("offset = %d, for rank %d, with ntasks = %d, and BLOCK_SIZE = %d\n", 
                offset, myrank, ntasks, BLOCK_SIZE);
    #pragma omp parallel for shared(X, Y, RHO, BLOCK_SIZE, offset, work_completed) private(probe,gene, tid)
    for (gene = offset; gene < BLOCK_SIZE; gene++) {
       for (probe = 0; probe < NUM_COLS; probe++) {
           double *x = getrowvec(X, gene, BUFFER_SIZE);
           double *y = getcolvec(Y, probe, BUFFER_SIZE);
           int x_nc = BUFFER_SIZE;
           int y_nr = BUFFER_SIZE;
           double avgx = mean(x, x_nc);
           sub(x, avgx, x_nc, xcentered);
                    
           double avgy = mean(y, y_nr);
           sub(y, avgy, y_nr, ycentered);
           
           prod(xcentered, ycentered, x_nc, prod_result);
           double sum_prod = sum(prod_result, x_nc);

           double stdX = stddev(x, avgx, x_nc);
           double stdY = stddev(y, avgy, y_nr);
           double rho  = sum_prod/((x_nc-1)*(stdX*stdY));

           //RHO[gene][probe] = rho;
           
           if (work_completed % 10000 == 0) {
               tid = omp_get_thread_num();
               printf("rank = %d, work = %d, result[%d,%d] from %d = %f\n", myrank, work_completed,
                    gene, probe, tid, rho);
           }
           work_completed++;
           free(x);
           free(y);
       }
    }
    printf("********* FINISHED **********\n");
    free(xcentered);
    free(ycentered);
    free(prod_result);
    free(X[0]);
    free(X);
    free(Y[0]);
    free(Y);
    //h5_write(RHO, NUM_ROWS, NUM_COLS, "rho_omp.h5", "/rho");
    free(RHO[0]);
    free(RHO);

  MPI_Finalize();
  return 0;
}
