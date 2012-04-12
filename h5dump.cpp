#include <stdlib.h>
#include <stdio.h>
#include "h5_read.h"



int main(int argc, char* argv[]) 
{

    unsigned long nr, nc;
    if (argc < 2) return 0;
    char * filename = argv[1];
    char * dataset  = argv[2];
    int sampleRate  = atoi(argv[3]);

    double **data= h5_read(filename, sampleRate, dataset,         &nr, &nc);
    printf("loaded data, num rows = %d, num cols = %d, val = %f\n", nr, nc, data[nr-1][nc-1]);
    double total_mem = nr*nc*sizeof(double*);
    double *RHO = (double*) malloc(total_mem);
    if (RHO == NULL) {
        printf("error allcoating array\n");
        exit(EXIT_FAILURE);
    }
     

    return 0;

}

