/******************************************************************************
* FILE: gistic.cpp
* DESCRIPTION:
*  Implementation for GISTIC algorithm.

@article{Beroukhim11122007,
    author = {Beroukhim, Rameen and Getz, Gad and Nghiemphu, Leia and Barretina, Jordi and Hsueh, Teli and Linhart, David and Vivanco, Igor and Lee, Jeffrey C. and Huang, Julie H. and Alexander, Sethu and Du, Jinyan and Kau, Tweeny and Thomas, Roman K. and Shah, Kinjal and Soto, Horacio and Perner, Sven and Prensner, John and Debiasi, Ralph M. and Demichelis, Francesca and Hatton, Charlie and Rubin, Mark A. and Garraway, Levi A. and Nelson, Stan F. and Liau, Linda and Mischel, Paul S. and Cloughesy, Tim F. and Meyerson, Matthew and Golub, Todd A. and Lander, Eric S. and Mellinghoff, Ingo K. and Sellers, William R.}, 
    title = {Assessing the significance of chromosomal aberrations in cancer: Methodology and application to glioma}, 
    volume = {104}, 
    number = {50}, 
    pages = {20007-20012}, 
    year = {2007}, 
    doi = {10.1073/pnas.0710052104}, 
    URL = {http://www.pnas.org/content/104/50/20007.abstract}, 
    eprint = {http://www.pnas.org/content/104/50/20007.full.pdf+html}, 
    journal = {Proceedings of the National Academy of Sciences} 
}


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

struct mem {
   void *data;
   size_t nbytes;
} ;
typedef struct mem mem;
void *my_malloc(size_t n)
{
   mem *p;
   p = (mem*)malloc( sizeof(*p) + n );
   if(p == NULL) return NULL;
   p->nbytes = n;
   p->data     =   p + 1;        // tricky!
   return p->data;
}
void my_free(void *memp)
{
    mem *p = (mem*) memp;
    if(p == NULL) return ;        // behave just like free()
    p--;                                    // tricky!
    free(p);
}
size_t GetSize(void *memp)
{
    mem *p = (mem*)memp;
    p--;                                 
    return p->nbytes;
}

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
sub(double *v, double val, int n) {

    double *result = (double*)calloc(BUFFER_SIZE,sizeof(double));;

    for (int i = 0; i < n; i++) {
        result[i] = v[i] - val;
    }

    return result;
}

void 
sub2(double *v, double val, int n) {

    for (int i = 0; i < n; i++) {
        v[i] = val - v[i];
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
    double max = 0;
    int chr_max = 0;
    for (int c = 0; c < 22; c++) {
        int len = chr_probes_end[c] - chr_probes_start[c];
        int start = chr_probes_start[c];
        int end   = chr_probes_end[c];
        printf("len = %d, start = %d, end = %d\n", len, start, end);

        double *amp = (double*)calloc(len,sizeof(double));

        //for (int sample = 0; sample < 2; sample++) {
        for (int sample = 0; sample < NUM_ROWS; sample++) {
           // the y elements contain a 1 if the threshold was met
           // otherwise it has a 0
           double *y = getampthresh(Y, sample, start, end, .1);
           // add 1 to the probe element if the sample met the threshold
           for (int i = 0; i < len; i++) {
               amp[i] += y[i];
           }

           free(y);
        }

        normalize(amp,len, BUFFER_SIZE);
        //totalamp[c] = &amp[0];
        // TODO: check this above
        for (int i = 0; i < len; i++) totalamp[c][i] = amp[i];
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

/*
void 
conv1(const double v1[], size_t n1, const double v2[], size_t n2, double r[])
{
    for (size_t n = 0; n < n1 + n2 - 1; n++)
        for (size_t k = 0; k < max(n1, n2); k++)
            r[n] += (k < n1 ? v1[k] : 0) * (n - k < n2 ? v2[n - k] : 0);
}
*/

double *
conv2(double *h, double *x, size_t len1, size_t len2) {

    int i, j;
    double sum = 0;
    double *y = (double*)calloc(len1+len2-1,sizeof(double));

    for (i = 0; i < (len1 + len2) - 1; i++) {
        sum = 0;
        for (j = 0; j <= i; j++) {
            sum += h[j] * x[i-j];
        }
        y[i] = sum;
    }
    free(h);
    return y;
}

double *
conv(double *Signal, size_t SignalLen,
          double *Kernel, size_t KernelLen) {

  size_t n;
  size_t kmin, kmax, k;

  double *Result = (double*)calloc(SignalLen+KernelLen,sizeof(double));
  #pragma omp parallel for shared(SignalLen, KernelLen, Signal, Kernel, Result) \
    private (n,k,kmin,kmax)
  for (n = 0; n < SignalLen + KernelLen - 1; n++)
  {

    Result[n] = 0;
    //Signal[j] = 0;

    kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
    kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

    for (k = kmin; k <= kmax; k++)
    {
      Result[n] += Signal[k] * Kernel[n - k];
      //Signal[j] += Signal[k] * Kernel[n - k];
    }
  }

  free(Signal);
  return Result;
}

double *
truncate(double *v, size_t n) {
    double * newarr = (double*)calloc(n,sizeof(double));

    memcpy(newarr, v, sizeof(double) * n);

    return newarr;
}

void
norm_by_sum(double *v, size_t n) {

    double sum = 0;
    for (int i = 0; i < n; i++) sum += v[i];
    for (int i = 0; i < n; i++) v[i] = v[i]/sum;;
}

void 
cumsum(double *v, size_t n) {

    double sum = 0;
    for (int i = 1; i < n; i++) v[i] = v[i-1] + v[i];
}


void 
setfiltergt(double *v, size_t n, double threshold, double val) {

    for (int i = 0; i < n; i++) 
        if (v[i] > threshold) v[i] = val;
    
}

void
genpvalues(double *v, size_t n) {

    norm_by_sum(v,n);
    cumsum(v,n);
    setfiltergt(v, n, 1.0, 1.0);
    sub2(v,1,n);

}

double **
mapsignificance(double **ampscores, double binwidth, double *pvalues, size_t n) {

    double **sig = create2dArray(23,7685);

    for (int i = 0; i < 22; i++) {
       for (int k = 0; k < 7685; k++) {
            int index = (int)(ampscores[i][k]/binwidth);
            sig[i][k] = pvalues[index];
       } 
    }
    return sig;
}
int comp (const void * elem1, const void * elem2) {
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}


double
find_next_min(double *v, size_t len, size_t idx) {

    double min = 999999;

    #pragma omp parallel for shared(v, min, len, idx) 
    for (size_t i = idx; i < len; i++) {
        if (v[i] < min) {
            min = v[i];
        }
    }

    return min;

}

int 
order_array(double *pvalues, size_t len) {

    int idx;
    double max = -1;
    #pragma omp parallel for shared(pvalues, len, idx, max) 
    for (int i = 0; i < len; i++) {
        if (pvalues[i] > max) {
            max = pvalues[i];
            idx = i;
        }
    }

    double * temp = (double*)calloc(len-1,sizeof(double));
    memcpy(temp,pvalues,(idx-1)*sizeof(double));
    memcpy(temp,&pvalues[idx+1],(len-idx)*sizeof(double));
    pvalues = temp;
    free(temp);
    return idx;
}
double *
benjaminihochberg(double **pvalues2d) {

    size_t len = NUM_ROWS * 7685;
    double *pvalues = (double*)calloc(len,sizeof(double));
    #pragma omp parallel for shared(pvalues,pvalues2d) 
    for (int i = 0; i < NUM_ROWS; i++) {
        for (int j = 0; j < 7685; j++) {
            pvalues[j+i*7685] = pvalues2d[i][j];
        }
    }


    qsort (pvalues, NUM_ROWS*7685, sizeof(*pvalues), comp);
    double *scaled = (double*)calloc(len,sizeof(double));

    #pragma omp parallel for shared(pvalues,scaled,len) 
    for (int i = 0; i < len; i++) {
        scaled[i] = (pvalues[i]*len)/i;
    }

    int *order = (int*)calloc(len,sizeof(int));
    for (int i = 0; i < len; i++) {
        int idx = order_array(pvalues, len);
        len--;
        order[i] = idx;
        printf("%d\n", i);
    }
    free(pvalues); 

    double *qvalues= (double*)calloc(len,sizeof(double));
    double min = 0;

    for (int i = 0; i < len-1; i++) {
        if (scaled[i] == min) {
            qvalues[i] = min;
            min = find_next_min(scaled, len, i+1);
        } else {
            qvalues[i] = min;
        }
    }

    qvalues[len-1] = scaled[len-1];

    return pvalues;
}

void
adjustpvalues(double **sig) {

    double *pvalues = benjaminihochberg(sig);
    free(pvalues);
}
int
main(int argc, char **argv) {

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

    double **Y = h5_read("filtered_probes.h5", 1, -1,-1,"/FilteredProbes", &y_nr, &y_nc);
    double binwidth = .001;

    int probe, sample, tid, work_completed;
    work_completed = 0;

    
    int BLOCK_SIZE = NUM_ROWS/ntasks;
    int offset = myrank*BLOCK_SIZE;
    int STOP_IDX = offset+BLOCK_SIZE;
    if (NUM_ROWS - STOP_IDX < BLOCK_SIZE)
        STOP_IDX = NUM_ROWS;

    double **totalamp = calc_amp_gscores(Y, chr_probes_start, chr_probes_end);
    double maxamp     = calculate_max(totalamp, 22, 7685);
    int    numbins    = (int)(maxamp/binwidth);
    double **amphist  = create2dArray(NUM_ROWS, numbins);
    double *bins      = get_bins(binwidth, maxamp);

    printf("max amp = %f\n", maxamp);

    // ------------ Start Generate Amplification Null -----------------------//
    // generate null distribution for amplification g-scores
    for (int i = 0; i < NUM_ROWS; i++) {

        // geneome wide scores for sample i
        double *samplecnv = getrowvec(Y, i, NUM_COLS);

        // threshold values
        norm_threshold(samplecnv, NUM_COLS, .1, NUM_ROWS); 

        // histogram
        amphist[i]      = hist(samplecnv, NUM_COLS, numbins, binwidth);

        free(samplecnv);

    }
    free(bins);

    printf("num bins = %d\n", numbins);

    double *nulldist = (double*)calloc(numbins,sizeof(double));
    // initialize null distribution
    nulldist = &amphist[0][0];

    printf("start convolution\n");
    // convolve distributions across samples
    for (int i = 1; i < NUM_ROWS; i++) {
           //printf("conv[%d] - %d\n", i, numbins*(i+1));
           nulldist = conv(nulldist, numbins*(i),&amphist[i][0], numbins);
    }
    printf("end convolution\n");

    // truncate
    nulldist = truncate(nulldist, numbins);
    printf("end truncate\n");
    // -------------- End Generate Amplification Null -------------------------//

    printf("start pvals\n");
    genpvalues(nulldist,numbins);
    //for (int i = 0; i < numbins; i++) printf("%f,",nulldist[i]);
    //printf("\n");
    printf("end pvals\n");

    printf("start mapsig\n");
    double **sig = mapsignificance(totalamp, binwidth, nulldist, numbins);
    printf("end mapsig\n");

    /*
    for (int c = 0; c < 22; c++) {
        for (int p = 0; p < 100; p++){
            printf("%10f,", sig[c][p]);
            if (p % 10 == 0) printf("\n");
        }
        printf("\n");
        printf("\n");
    }

    //adjustpvalues(sig);
    for (int c = 0; c < 22; c++) {
        int len = chr_probes_end[c] - chr_probes_start[c];
        int start = chr_probes_start[c];
        int end   = chr_probes_end[c];
        for (int j = 0; j < len; j++) {
             if (sig[c][j] < .000001) {
                 int index = start + j;

                 //printf("%s,%d,%.15f\n", records[index].gene, c+1,sig[c][j]);
                 printf("%s\n", records[index].gene);

             }
        }

    }
    */

    printf("********* %d FINISHED **********\n", myrank);
    endtime   = MPI_Wtime();
    printf("rank %d - elapse time - %f\n",myrank, endtime-starttime);

    free(Y[0]);
    free(Y);

    free(sig[0]);
    free(sig);

    free(nulldist);

    free(totalamp[0]);
    free(totalamp);

    free(exprannot);
    free(records);

    free(chr_probes_start);
    free(chr_probes_end);

  MPI_Finalize();
  return 0;
}
