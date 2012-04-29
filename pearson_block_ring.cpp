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

//#include "h5_read.h"
#include "h5_hyperslab.h"
#include "h5_write.h"

#define WORKTAG 1
#define DIETAG 2

#define BUFFER_SIZE 383

#define LINE_SIZE 32
//#define NUM_ROWS 1594
//#define NUM_COLS 7899
#define NUM_ROWS 15946
#define NUM_COLS 78995
//#define NUM_ROWS 10
//#define NUM_COLS 20
#define WORK_SIZE NUM_ROWS*NUM_COLS


int NUM_PROBES = 0;

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

exprinfo *
get_gene_expr_annot() {

    int numlines = 0;
    char **glines = NULL;
    char line[LINE_SIZE];

    exprinfo *exprannot = NULL;
    FILE *f = fopen("gene_list.txt", "r");

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
    free(glines[0]);
    free(glines);

    return exprannot;

}

probeinfo *
get_probe_info() {

    int numlines = 0;
    char **glines = NULL;
    char line[LINE_SIZE];

    char **lines = NULL;
    probeinfo *records = NULL;

    FILE *f = fopen("probe_id_mapping.txt", "r");

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
    NUM_PROBES = numlines;
    free(lines[0]);
    free(lines);

    fclose(f);

    return records;

}


int
main(int argc, char **argv)
{
    double starttime, endtime;
    int ntasks, myrank;
    char name[128];                      
    int namelen;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

    MPI_Get_processor_name(name,&namelen);

    /* Get a few OpenMP parameters.                                               */
    int O_P  = omp_get_num_procs();          /* get number of OpenMP processors       */
    int O_T  = omp_get_num_threads();        /* get number of OpenMP threads          */
    int O_ID = omp_get_thread_num();         /* get OpenMP thread ID                  */
    //printf("name:%s   M_ID:%d  O_ID:%d  O_P:%d  O_T:%d\n", name,myrank,O_ID,O_P,O_T);
    exprinfo *exprannot = get_gene_expr_annot();
    probeinfo *records = get_probe_info();

    int gene, probe, tid, work_completed;
    work_completed = 0;

    
    int BLOCK_ROW_SIZE = NUM_ROWS/ntasks;
    int BLOCK_COL_SIZE = NUM_COLS/ntasks;
    int y_offset = BLOCK_COL_SIZE*myrank;

    double **X = hyperslabread("x.h5", "/X", BLOCK_ROW_SIZE, BUFFER_SIZE, 0,0);
    double **Y = hyperslabread("filtered_probesT.h5", "/FilteredProbes", 
                                BLOCK_COL_SIZE, BUFFER_SIZE, y_offset,0);
    //double **RHO = create2dArray(BLOCK_ROW_SIZE, BLOCK_COL_SIZE);

    printf("for rank %d, with ntasks = %d, and BLOCK_SIZE = %d, BLOCK_COLS=%d, BUF=%d\n", 
                 myrank, ntasks, BLOCK_ROW_SIZE, BLOCK_COL_SIZE, BUFFER_SIZE);

    int num_sig_found = 0; 

    starttime = MPI_Wtime();

    for (int task = 0; task < ntasks; task++) {

        #pragma omp parallel \
         for shared(X, Y, BLOCK_ROW_SIZE, BLOCK_COL_SIZE,ntasks,task) \
         private(probe,gene, tid )
        for (gene = 0; gene < BLOCK_ROW_SIZE; gene++) {
           for (probe = 0; probe < BLOCK_COL_SIZE; probe++) {
               double *x = getrowvec(X, gene, BUFFER_SIZE);
               double *y = getrowvec(Y, probe, BUFFER_SIZE);

               double avgx = mean(x, BUFFER_SIZE);
               double * xcentered = sub(x, avgx, BUFFER_SIZE);
                        
               double avgy = mean(y, BUFFER_SIZE);
               double * ycentered = sub(y, avgy, BUFFER_SIZE);
               
               double * prod_result = prod(xcentered, ycentered, BUFFER_SIZE);
               double sum_prod = sum(prod_result, BUFFER_SIZE);

               double stdX = stddev(x, avgx, BUFFER_SIZE);
               double stdY = stddev(y, avgy, BUFFER_SIZE);
               double rho  = sum_prod/((BUFFER_SIZE-1)*(stdX*stdY));

               //RHO[gene][probe] = rho;
               double zscore = ztest(rho, BUFFER_SIZE);
               
               //if (work_completed % 100000 == 0) {
                   //tid = omp_get_thread_num();
                   //printf("rank = %d, work = %d, result[%d,%d] from %d = %f\n", myrank, work_completed,
                    //    gene, probe, tid, rho);
               //}
               
               //work_completed++;
               free(x);
               free(y);
               free(xcentered);
               free(ycentered);
               free(prod_result);

           }
        }
        //f = fopen("significant.txt", "a");
        
        //#pragma omp parallel for shared(RHO,exprannot, records)
        //for (int i = 0; i < BLOCK_ROW_SIZE; i++) {
            //for (int j = 0; j < BLOCK_COL_SIZE; j++) {
            //double zscore = ztest(RHO[i][j],BUFFER_SIZE);
                //if (zscore > 5.0) {

                //fprintf(f, "%d,%d,%f,%f,%d,%s,%d,%s\n", 
                 //       i, j, zscore, RHO[i][j], records[j].chr, records[j].gene, exprannot[i].chr,exprannot[i].gene);
                    //printf("%d,%d,%f,%f,%d,%s,%d,%s\n", 
                     //   i, j, zscore, RHO[i][j], records[j].chr, records[j].gene, exprannot[i].chr,exprannot[i].gene);
                //}
            ////}
        //}
        //fclose(f);
    
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Request request;
        MPI_Status status;
        int prev = myrank - 1;
        int next = myrank + 1;
        if (myrank == 0) prev = ntasks - 1;
        if (myrank == ntasks - 1) next = 0;
        printf("rank[%d] sending to %d, and receiving from %d for iteration %d\n", myrank, next, prev,task+1);
        MPI_Isend(&X[0][0], BLOCK_ROW_SIZE*BUFFER_SIZE, MPI_DOUBLE, next, 1, MPI_COMM_WORLD, &request);
        MPI_Irecv(&X[0][0], BLOCK_ROW_SIZE*BUFFER_SIZE, MPI_DOUBLE, prev, 1, MPI_COMM_WORLD, &request);
        MPI_Wait( &request, &status);
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0) {
            printf("completed %d out of %d iterations\n",task+1,ntasks);
        }
    }
    printf("********* %d FINISHED **********\n", myrank);

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
    //free(RHO[0]);
    //free(RHO);
    free(exprannot);
    free(records);

  MPI_Finalize();
  return 0;
}
