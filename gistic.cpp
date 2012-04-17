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
#include <stddef.h>

#include "h5_read.h"
#include "h5_write.h"

#define WORKTAG 1
#define DIETAG 2

#define BUFFER_SIZE 383

#define LINE_SIZE 32
#define NUM_ROWS 383
#define NUM_COLS 78995
//#define NUM_ROWS 20004
//#define NUM_COLS 78996
//#define NUM_COLS 52664
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


double 
mean(double *v, int n) {

    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += v[i];

    }

    return sum/n;
}

double 
stddev(double *v, double avg, int n) {

    double sum = 0;
    for (int i = 0; i < n; i++) 
        sum += (v[i] - avg)*(v[i] - avg);

    return sqrt(sum/(n-1));
}

double *
sub(double *v, double avg, int n) {

    double *result = (double*)calloc(BUFFER_SIZE,sizeof(double));;

    for (int i = 0; i < n; i++) {
        result[i] = v[i] - avg;
    }

    return result;
}

void 
sub(double *v, double avg, int n, double *result) {

    for (int i = 0; i < n; i++) {
        result[i] = v[i] - avg;
    }
}

double *
prod(double *x, double *y, int n) {

    double *result = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    for (int i = 0; i < n; i++) {
        result[i] = x[i]*y[i];
    }

    return result;

}

void 
prod(double *x, double *y, int n, double *result) {

    for (int i = 0; i < n; i++) {
        result[i] = x[i]*y[i];
    }

}

double 
sum(double *v, int n) {

    double result = 0;

    for (int i = 0; i < n; i++) {
        result += v[i];
    }

    return result;

}

double 
ztest(double r, int n) {

    return r * sqrt( (n-2)/( 1 - (r*r)));

}

double* 
getrowvec(double **m, int row, int nc) {

    double *result = (double*)calloc(nc,sizeof(double));;
    memcpy(result, m[row], sizeof(double) * nc);
    return result;

}

double* 
getampthresh(double **m, int sample, int start, int end, double thresh) {

    double *result = (double*)calloc(end-start,sizeof(double));
    int k = 0;
    for (int i = start; i < end; i++) {
        if (m[sample][i] > thresh) 
            result[k++] = 1;
        //printf("m[%d][%d] = %f\n", sample, i, m[sample][i]);
    }
    return result;

}

double* 
getcolvec(double **m, int col, int nr) {

    double *result  = (double*)calloc(BUFFER_SIZE,sizeof(double));;
    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }

    return result;

}

void
normalize(double *v, int n, int samples) {

    int i;
    for (i = 0; i < n; i++) {
        v[i] = v[i]/samples;
    }

}

void 
getrowvec(double **m, int row, int nc, double *result) {


    memcpy(result, m[row], sizeof(double) * nc);

}

void 
getcolvec(double **m, int col, int nr, double *result) {

    int i;
    for (i = 0; i < nr; i++) {
        result[i] = m[i][col];
    }

}

int **
create2dArrayInt(int rows, int cols) {

    unsigned long total_mem = (rows * cols);
    int **array = (int**) calloc(rows, sizeof(int*));
    if (array == NULL) {
        printf("error allcoating array\n");
        exit(EXIT_FAILURE);
    }
    array[0] = (int*)calloc( total_mem, sizeof(int) );
    if (array[0] == NULL) {
        printf("error allcoating array\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned long i=1; i < rows; i++) array[i] = array[0]+i*cols;

    return array;

}

double **
create2dArray(int rows, int cols) {

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

char *
trimwhitespace(char *str)
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

int *
get_chr_probe_start(probeinfo *records) {

    int lastchr = 1;
    int laststart = 0;
    int * result = (int *)calloc(NUM_PROBES,sizeof(int));
    int n = 0;
    for (int i = 0; i < NUM_PROBES; i++) {
        if (records[i].chr != lastchr) {
            result[lastchr-1] = laststart;
            laststart = n;
            lastchr = records[i].chr;
        }
        n++;
    }
    result[lastchr-1] = laststart;

    return result;

}

int *
get_chr_probe_end(probeinfo *records) {

    int lastchr = 1;
    int * result = (int *)calloc(NUM_PROBES,sizeof(int));
    int n = 0;
    for (int i = 0; i < NUM_PROBES; i++) {
        if (records[i].chr != lastchr) {
            result[lastchr-1] = n;
            lastchr = records[i].chr;
        }
        n++;
    }
    result[lastchr-1] += n;

    return result;

}

double **
calc_amp_gscores(double **Y, int *chr_probes_start, int *chr_probes_end) {

    //#pragma omp parallel \
     //for shared(Y,  BLOCK_SIZE, offset, work_completed) \
     //private(probe,sample, tid)
    double **totalamp = create2dArray(22,7685);
    for (int c = 0; c < 22; c++) {
        int len = chr_probes_end[c] - chr_probes_start[c];
        int start = chr_probes_start[c];
        int end   = chr_probes_end[c];
        //printf("len = %d, start = %d, end = %d\n", len, start, end);

        double *amp = (double*)calloc(len,sizeof(double));

        //for (sample = offset; sample < STOP_IDX; sample++) {
        for (int sample = 0; sample < NUM_ROWS; sample++) {
           double *y = getampthresh(Y, sample, start, end, .1);
           for (int i = 0; i < len; i++) {
               amp[i] += y[i];
           }

           free(y);
        }

        normalize(amp,len, BUFFER_SIZE);
        totalamp[c] = &amp[0];
        // TODO: check this above
        //for (int i = 0; i < len; i++) totalamp[c][i] = amp[i];
        free(amp);
        
    }

    return totalamp;
}

double
calculate_max(double **ampscores, int rows, int cols) {

    double max = -1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 100; j++) {
            if (ampscores[i][j] > max) {
                max = ampscores[i][j];
            }
        }

    }

    return max;
}

void
norm_threshold(double *cndata, int len, double threshold, int numsamples) {

    for (int i = 0; i < len; i++) {
        if (cndata[i] >= threshold) {
            cndata[i] = cndata[i]/numsamples;
        } else {
            cndata[i] = 0;
        }
    }

}

double *
get_bins(double binwidth, double max) {

    int numbins = (int)(max/binwidth);
    double *bins = (double*)calloc(numbins, sizeof(double)); 
    for (int i = 0; i < numbins; i++) {
        bins[i] += binwidth;
    }

    return bins;
}

double *
hist(double *data, int datalen, int bins, double binwidth){

    double *hist= (double*)calloc(bins, sizeof(double));

    // create histogram
    for (int i=0; i < datalen; ++i){
        //int bin=int( (data[i]-min)/((max-min)/(bins)) );
        int bin=(int)(data[i]/binwidth);
        hist[bin]++;
    }

    // normalize
    for (int i = 0; i < bins; i++) {
        hist[i] = hist[i]/datalen;
    }
    return hist;
}

int max(int n1, int n2) {
    return (n1 > n2 ? n1 : n2);
}

void 
conv1(const double v1[], size_t n1, const double v2[], size_t n2, double r[])
{
    for (size_t n = 0; n < n1 + n2 - 1; n++)
        for (size_t k = 0; k < max(n1, n2); k++)
            r[n] += (k < n1 ? v1[k] : 0) * (n - k < n2 ? v2[n - k] : 0);
}

void conv(const double *Signal, size_t SignalLen,
          const double *Kernel, size_t KernelLen,
          double *Result)
{
  size_t n;

  for (n = 0; n < SignalLen + KernelLen - 1; n++)
  {
    size_t kmin, kmax, k;

    Result[n] = 0;

    kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
    kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

    for (k = kmin; k <= kmax; k++)
    {
      Result[n] += Signal[k] * Kernel[n - k];
    }
  }
}

int
main(int argc, char **argv)

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
    //printf("name:%s   M_ID:%d  O_ID:%d  O_P:%d  O_T:%d\n", name,myrank,O_ID,O_P,O_T);

    exprinfo *exprannot = get_gene_expr_annot();
    probeinfo *records = get_probe_info();
    int *chr_probes_start = get_chr_probe_start(records);
    int *chr_probes_end  = get_chr_probe_end(records);
    //for (int i = 0; i < 22; i++ ) printf("%d\n", chr_probes_end[i] - chr_probes_start[i]);

    unsigned long y_nr, y_nc;

    double **Y = h5_read("filtered_probes.h5", 1, "/FilteredProbes", &y_nr, &y_nc);
    double binwidth = .01;

    int probe, sample, tid, work_completed;
    work_completed = 0;

    
    int BLOCK_SIZE = NUM_ROWS/ntasks;
    int offset = myrank*BLOCK_SIZE;
    int STOP_IDX = offset+BLOCK_SIZE;
    if (NUM_ROWS - STOP_IDX < BLOCK_SIZE)
        STOP_IDX = NUM_ROWS;

    double ** totalamp  = calc_amp_gscores(Y, chr_probes_start, chr_probes_end);
    double maxamp       = calculate_max(totalamp, 22, 7685);
    printf("max amp = %f\n", maxamp);
    int     numbins     = (int)(maxamp/binwidth);
    double **sample_hist   = create2dArray(NUM_ROWS, numbins);

    for (int i = 0; i < NUM_ROWS; i++) {
        printf("sample: %d\n",i);
        double *samplecnv = getrowvec(Y, i, NUM_COLS);
        norm_threshold(samplecnv, NUM_COLS, .1, NUM_ROWS); 
        double *bins = get_bins(binwidth, maxamp);
        double * h      = hist(samplecnv, NUM_COLS, numbins, binwidth);
        sample_hist[i]   = &h[0];
        for (int k = 0; k < numbins; k++) printf("%f, ", sample_hist[i][k]);
        printf("\n");

        if (i > 3)
            break;



        // GenerateAmplificationNull
        //getampthresh(double **m, int sample, int start, int end, double thresh) {
    }
    printf ("max amp = %f\n", maxamp);

    /*
    for (int c = 0; c < 22; c++) {
        for (int p = 0; p < 20; p++){
            printf("%f,", c,p,totalamp[c][p]);
        }
        printf("\n");
    }
    */

    printf("********* %d FINISHED **********\n", myrank);
    endtime   = MPI_Wtime();
    printf("rank %d - elapse time - %f\n",myrank, endtime-starttime);
    free(Y[0]);
    free(Y);

    /*
    f = fopen("significant.txt", "a");
    //#pragma omp parallel for shared(RHO,exprannot, records)
    for (int i = 0; i < NUM_ROWS; i++) {
        for (int j = 0; j < NUM_COLS; j++) {
            double zscore = ztest(RHO[i][j],BUFFER_SIZE);
            if (zscore > 5.0) {

                fprintf(f, "%d,%d,%f,%f,%d,%s,%d,%s\n", 
                        i, j, zscore, RHO[i][j], records[j].chr, records[j].gene, exprannot[i].chr,exprannot[i].gene);
                if (i*j%1000 == 0)
                    printf("%d,%d,%f,%f,%d,%s,%d,%s\n", 
                        i, j, zscore, RHO[i][j], records[j].chr, records[j].gene, exprannot[i].chr,exprannot[i].gene);
            }
        }
    }
    fclose(f);
    */
    //for (int i = 0; i < NUM_ROWS; i++) {
        //printf("%s,%d\n", exprannot[i].gene, exprannot[i].chr);
    //}
    free(totalamp[0]);
    free(totalamp);
    free(exprannot);
    free(records);
    free(chr_probes_start);
    free(chr_probes_end);

  MPI_Finalize();
  return 0;
}
