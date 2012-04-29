/************************************************************
  
  This example shows how to write and read a hyperslab.  It 
  is derived from the h5_read.c and h5_write.c examples in 
  the "Introduction to HDF5".

 ************************************************************/
 
#include "hdf5.h"
#include <stdlib.h>
#include "h5_hyperslab.h"

#define RANK_OUT     2
#define DEBUG 0


double **
hyperslabread(const char * FILE, const char *DATASETNAME,
        int SLAB_ROWS, int SLAB_COLS,
        int offsetrow, int offsetcol)
{

    /* 
     * Data  and output buffer initialization. 
     */
    hid_t       file, dataset;         /* handles */
    hid_t       dataspace;   
    hid_t       memspace; 
    hsize_t     dimsm[2];              /* memory space dimensions */
    hsize_t     dims_out[2];           /* dataset dimensions */      
    herr_t      status;                             

    //double      data_out[SLAB_ROWS][SLAB_COLS]; /* output buffer */
    double ** data_out2;
   
    hsize_t     count[2];              /* size of the hyperslab in the file */
    hsize_t     offset[2];             /* hyperslab offset in the file */
    hsize_t     count_out[2];          /* size of the hyperslab in memory */
    hsize_t     offset_out[2];         /* hyperslab offset in memory */

    int         i, j, k, status_n, rank, rows, cols;

    rows = SLAB_ROWS;
    cols = SLAB_COLS;

/*************************** BEGIN  *******************************/

/* Allocate memory for new integer array[row][col]. First
   allocate the memory for the top-level array (rows).  Make
   sure you use the sizeof a *pointer* to your data type. */

    data_out2 = (double**) malloc(rows*sizeof(double*));
    

/* Allocate a contiguous chunk of memory for the array data values.  
   Use the sizeof the data type. */

    data_out2[0] = (double*)calloc( cols*rows, sizeof(double) );

/* Set the pointers in the top-level (row) array to the
   correct memory locations in the data value chunk. */

    for (i=1; i < rows; i++) data_out2[i] = data_out2[0]+i*cols;

/************************* END *************************************/


/*************************************************************  

  This reads the hyperslab from the sds.h5 file just 
  created, into a 2-dimensional plane of the 3-dimensional 
  array.

 ************************************************************/  

    /*
    for (j = 0; j < SLAB_ROWS; j++) {
        for (i = 0; i < SLAB_COLS; i++) {
                data_out2[j][i] = 0;
        }
    } 
    */
 
    /*
     * Open the file and the dataset.
     */
    file = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (DEBUG)
        printf("before open dataset\n");
    dataset = H5Dopen (file, DATASETNAME, H5P_DEFAULT);

    dataspace = H5Dget_space (dataset);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims (dataspace);
    status_n  = H5Sget_simple_extent_dims (dataspace, dims_out, NULL);
    if (DEBUG)
        printf("\nRank: %d\nDimensions: %lu x %lu \n", rank,
            (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

    /* 
     * Define hyperslab in the dataset. 
     */
    offset[0] = offsetrow;
    offset[1] = offsetcol;
    count[0]  = SLAB_ROWS;
    count[1]  = SLAB_COLS;
    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL, 
                                  count, NULL);


    /*
     * Define the memory dataspace.
     */
    dimsm[0] = SLAB_ROWS;
    dimsm[1] = SLAB_COLS;
    memspace = H5Screate_simple (RANK_OUT, dimsm, NULL);   

    /* 
     * Define memory hyperslab. 
     */
    offset_out[0] = 0;  // which row to start from
    offset_out[1] = 0;  // where column to start from

    count_out[0]  = SLAB_ROWS;  // how many to get from slab
    count_out[1]  = SLAB_COLS;  // how many to get from slab

    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL, 
                                  count_out, NULL);

    /*
     * Read data from hyperslab in the file into the hyperslab in 
     * memory and display.
     */
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                      H5P_DEFAULT, &data_out2[0][0]);



    /*
     * Close and release resources.
     */
    H5Dclose (dataset);
    H5Sclose (dataspace);
    H5Sclose (memspace);
    H5Fclose (file);

    return data_out2;

}     

