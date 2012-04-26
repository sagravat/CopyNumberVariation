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

    printf("before load\n");
    double **data= h5_read(filename, sampleRate, -1, -1, dataset,&nr, &nc);
    printf("loaded data, num rows = %d, num cols = %d\n", nr, nc);

    /*
    double total_mem = nr*nc*sizeof(double*);
    double *RHO = (double*) malloc(total_mem);
    if (RHO == NULL) {
        printf("error allcoating array\n");
        exit(EXIT_FAILURE);
    }

    */ 
    free(data[0]);
    free(data);

    return 0;

}

