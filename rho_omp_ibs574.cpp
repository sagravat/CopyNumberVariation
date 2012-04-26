/******************************************************************************
* FILE: rho.cpp
* DESCRIPTION:
*  Hybrid MPI/OpenMP implementation for calculating Pearson Correlation.
*   
*   
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
#define NUM_ROWS 15946
#define NUM_COLS 78995
//#define NUM_ROWS 20004
//#define NUM_COLS 78996
//#define NUM_COLS 52664
//#define NUM_ROWS 10
//#define NUM_COLS 20
#define WORK_SIZE NUM_ROWS*NUM_COLS


typedef struct {
    int probeid;
    char * gene;
    int chr;
} probeinfo;

typedef struct {
    char * gene;
    int chr;
} exprinfo;

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

double *sub(double *v, double avg, int n) {

    double *result = (double*)calloc(BUFFER_SIZE,sizeof(double));;

    for (int i = 0; i < n; i++) {
        result[i] = v[i] - avg;
    }

    return result;
}

void sub(double *v, double avg, int n, double *result) {

    for (int i = 0; i < n; i++) {
        result[i] = v[i] - avg;
    }
}

double *prod(double *x, double *y, int n) {

    double *result = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    for (int i = 0; i < n; i++) {
        result[i] = x[i]*y[i];
    }

    return result;

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

double ztest(double r, int n) {

    return r * sqrt( (n-2)/( 1 - (r*r)));

}

double* getrowvec(double **m, int row, int nc) {

    double *result = (double*)calloc(BUFFER_SIZE,sizeof(double));;

    memcpy(result, m[row], sizeof(double) * nc);
    return result;

}

double* getcolvec(double **m, int col, int nr) {

    double *result  = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }

    return result;

}

void getrowvec(double **m, int row, int nc, double *result) {


    memcpy(result, m[row], sizeof(double) * nc);

}

void getcolvec(double **m, int col, int nr, double *result) {

    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }

}

double **create2dArray(int rows, int cols) {

    unsigned long total_mem = (rows * cols);
    double **array = (double**) calloc(rows, sizeof(double*));
    if (array == NULL) {
        printf("error allcoating array\n");
        exit(EXIT_FAILURE);
    }
    array[0] = (double*)calloc( total_mem, sizeof(double) );
    if (array[0] == NULL) {
        printf("error allcoating array\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned long i=1; i < rows; i++) array[i] = array[0]+i*cols;

    return array;

}

char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

int
main(int argc, char **argv)
{
    double starttime, endtime;
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
    //printf("name:%s   M_ID:%d  O_ID:%d  O_P:%d  O_T:%d\n", name,myrank,O_ID,O_P,O_T);

    FILE *f;
    char line[LINE_SIZE];
    int numlines = 0;

    exprinfo *exprannot = NULL;
    char **glines = NULL;

    f = fopen("gene_list.txt", "r");
    while(fgets(line, LINE_SIZE, f)) {
        glines = (char**)realloc(glines, sizeof(char*)*(numlines+1));
        glines[numlines] = strdup(line);

        char *pch = strtok (line,",");
        char * gene = pch;

        pch = strtok (NULL, ",");
        int chr = atoi(trimwhitespace(pch));

        exprannot = (exprinfo*)realloc(exprannot,sizeof(exprinfo)*(numlines+1));

        exprannot[numlines].gene      = strdup(gene);
        exprannot[numlines].chr       = chr;

        if (!exprannot) printf("not allcoated\n");
        numlines++;
    }
    fclose(f);
    f = fopen("probe_id_mapping.txt", "r");
    numlines = 0;

    free(glines[0]);
    free(glines);

    probeinfo *records = NULL;
    char **lines = NULL;

    while(fgets(line, LINE_SIZE, f)) {              
        lines = (char**)realloc(lines, sizeof(char*)*(numlines+1));
        lines[numlines] = strdup(line);

        char *pch = strtok (line,",");
        int probeid = atoi(pch);

        pch = strtok (NULL, ",");
        char * gene = pch;

        pch = strtok (NULL, ",");
        int chr = atoi(trimwhitespace(pch));

        records = (probeinfo*)realloc(records,sizeof(probeinfo)*(numlines+1));

        records[numlines].probeid   = probeid;
        records[numlines].chr       = chr;
        records[numlines].gene      = strdup(gene);
        if (!records) printf("not allcoated\n");
        numlines++;

    }
    free(lines[0]);
    free(lines);

    fclose(f);


    int NUM_GENES = numlines;
    unsigned long x_nr, x_nc, y_nr, y_nc;

    double **X = h5_read("x.h5", 1, "/X",         &x_nr, &x_nc);
    double **Y = h5_read("filtered_probes.h5", 1, "/FilteredProbes", &y_nr, &y_nc);

    //printf("loaded X, num rows = %d, num cols = %d\n", x_nr, x_nc);
    //printf("loaded Y, num rows = %d, num cols = %d\n", y_nr, y_nc);

    unsigned long total_mem = (x_nr * y_nc);
    double **RHO   = create2dArray(x_nr, y_nc);
    
    int gene, probe, tid, work_completed;
    work_completed = 0;

    
    int BLOCK_SIZE = NUM_ROWS/ntasks;
    int offset = myrank*BLOCK_SIZE;
    int STOP_IDX = offset+BLOCK_SIZE;
    if (NUM_ROWS - STOP_IDX < BLOCK_SIZE)
        STOP_IDX = NUM_ROWS;

   // printf("offset = %d, for rank %d, with ntasks = %d, and BLOCK_SIZE = %d, STOP_IDX = %d\n", 
    //            offset, myrank, ntasks, BLOCK_SIZE, STOP_IDX);

    int num_sig_found = 0; 
    starttime = MPI_Wtime();

    #pragma omp parallel \
     for shared(X, Y, RHO, BLOCK_SIZE, offset, work_completed) \
     private(probe,gene, tid)
    for (gene = offset; gene < STOP_IDX; gene++) {
       for (probe = 0; probe < NUM_COLS; probe++) {
           double *x = getrowvec(X, gene, BUFFER_SIZE);
           double *y = getcolvec(Y, probe, BUFFER_SIZE);

           double avgx = mean(x, BUFFER_SIZE);
           double * xcentered = sub(x, avgx, BUFFER_SIZE);
                    
           double avgy = mean(y, BUFFER_SIZE);
           double * ycentered = sub(y, avgy, BUFFER_SIZE);
           
           double * prod_result = prod(xcentered, ycentered, BUFFER_SIZE);
           double sum_prod = sum(prod_result, BUFFER_SIZE);

           double stdX = stddev(x, avgx, BUFFER_SIZE);
           double stdY = stddev(y, avgy, BUFFER_SIZE);
           double rho  = sum_prod/((BUFFER_SIZE-1)*(stdX*stdY));

           RHO[gene][probe] = rho;
           if (work_completed % 10000 == 0) {
               tid = omp_get_thread_num();
               printf("rank = %d, work = %d, result[%d,%d] from %d = %f\n", myrank, work_completed,
                    gene, probe, tid, rho);
           }
           work_completed++;
           free(x);
           free(y);
           free(xcentered);
           free(ycentered);
           free(prod_result);

       }
    }
    //printf("********* %d FINISHED **********\n", myrank);

    //f = fopen("significant.txt", "a");
    #pragma omp parallel for shared(RHO,exprannot, records)
    for (int i = 0; i < NUM_ROWS; i++) {
        for (int j = 0; j < NUM_COLS; j++) {
            double zscore = ztest(RHO[i][j],BUFFER_SIZE);
            /*
            if (zscore > 5.0) {

                fprintf(f, "%d,%d,%f,%f,%d,%s,%d,%s\n", 
                        i, j, zscore, RHO[i][j], records[j].chr, records[j].gene, exprannot[i].chr,exprannot[i].gene);
                if (i*j%1000 == 0)
                    printf("%d,%d,%f,%f,%d,%s,%d,%s\n", 
                        i, j, zscore, RHO[i][j], records[j].chr, records[j].gene, exprannot[i].chr,exprannot[i].gene);
            }
            */
        }
    }
    //fclose(f);
    endtime   = MPI_Wtime();
    free(X[0]);
    free(X);
    free(Y[0]);
    free(Y);
    printf("rank %d - elapse time - %f\n",myrank, endtime-starttime);
    //for (int i = 0; i < NUM_ROWS; i++) {
        //printf("%s,%d\n", exprannot[i].gene, exprannot[i].chr);
    //}
    //printf("rank %d FINISHED\n",myrank);
    //h5_write(RHO, NUM_ROWS, NUM_COLS, "rho_omp.h5", "/rho");
    free(RHO[0]);
    free(exprannot);
    free(records);
    free(RHO);

  MPI_Finalize();
  return 0;
}
