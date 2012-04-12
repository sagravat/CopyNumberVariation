/*****************************************************************
  
   h5_write.cpp - Write out dataset using dynamic array

******************************************************************/

#include <hdf5.h>
#include <stdlib.h>
#include "h5_write.h"

#define RANK 2

void
h5_write(double **data, int rows, int cols, const char * filename, const char * dataset_name)
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[2];              /* dataset dimensions */
    herr_t      status;                             
    unsigned long i, j;

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dimsf[0] = rows;
    dimsf[1] = cols;
    dataspace = H5Screate_simple(RANK, dimsf, NULL); 

    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    status = H5Tset_order(datatype, H5T_ORDER_LE);

    dataset = H5Dcreate(file, dataset_name, datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		      H5P_DEFAULT, &data[0][0]);

    //free(data[0]);
    //free(data);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Fclose(file);
 
}     

