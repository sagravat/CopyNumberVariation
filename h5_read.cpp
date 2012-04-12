/******************************************************************
  
 h5_read.cpp  - Read dataset using dynamic arrays

******************************************************************/
 
#include "hdf5.h"
#include <stdlib.h>
#include <omp.h>

#define DEBUG 0

double **
h5_read(const char * filename, int sampleRate, const char * dataset_name, unsigned long *nr, unsigned long *nc)
{
    hid_t       file, dataset;         /* handles */
    hid_t       datatype, dataspace;   
    H5T_class_t clazz;                 /* data type class */
    H5T_order_t order;                 /* data order */
    size_t      size;                  /*
				                        * size of the data element	       
				                        * stored in file
				                        */
    hsize_t     dims_out[2];           /* dataset dimensions */      
    herr_t      status;                             

    double          **data_out; 
    //double          **data_matrix; 
    unsigned long rows, cols, resampledSize;
    unsigned long   i, j, k, status_n, rank;


    /*  Open the file and the dataset  */
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);

    /* Get datatype and dataspace handles and then query
       dataset class, order, size, rank and dimensions  */
    datatype  = H5Dget_type(dataset);     /* datatype handle */ 
    clazz     = H5Tget_class(datatype);
    if (DEBUG)
        if (clazz == H5T_FLOAT) printf("Data set has FLOAT type \n");
    order     = H5Tget_order(datatype);

    if (DEBUG)
        if (order == H5T_ORDER_LE) printf("Little endian order \n");

    size  = H5Tget_size(datatype);
    if (DEBUG)
        printf(" Data size is %d \n", size);

    dataspace = H5Dget_space(dataset);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
   // if (DEBUG)
        printf("rank %d, dimensions %lu rows x %lu cols\n", rank,
	        (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

    rows = dims_out[0];
    cols = dims_out[1];
    //*nc = rows;
    //*nr = cols;
    *nc = cols;
    *nr = rows;
    resampledSize = (cols)/sampleRate;
    //if (DEBUG)
        printf("resampled size = %d\n", resampledSize);

    //m.set_size(resampled_size, cols);

/*************************** BEGIN  *******************************/

/* Allocate memory for new integer array[row][col]. First
   allocate the memory for the top-level array (rows).  Make
   sure you use the sizeof a *pointer* to your data type. */

    data_out = (double**) malloc(rows*sizeof(double*));
    //data_matrix = (double**) malloc(rows*sizeof(double*));
    

/* Allocate a contiguous chunk of memory for the array data values.  
   Use the sizeof the data type. */

    data_out[0] = (double*)malloc( cols*rows*sizeof(double) );
    //data_matrix[0] = (double*)malloc( resampledSize*rows*sizeof(double) );

/* Set the pointers in the top-level (row) array to the
   correct memory locations in the data value chunk. */

    for (i=1; i < rows; i++) data_out[i] = data_out[0]+i*cols;
    //for (i=1; i < rows; i++) data_matrix[i] = data_matrix[0]+i*resampledSize;

/************************* END *************************************/

    /* Initialize data_out array to 0 */
    for (j = 0; j < dims_out[0]; j++) {
        for (i = 0; i < dims_out[1]; i++) {
                data_out[j][i] = 0;
        }
    }
    
        /*
    for (j = 0; j < dims_out[0]; j++) {
        for (i = 0; i < resampledSize; i++) {
                data_matrix[j][i] = 0;
        }
    }
*/
    

    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, &data_out[0][0]);
    
    if (DEBUG) {
        printf("rows = %d\n", rows);
        printf("data = %f\n", data_out[0][0]);
        printf("data = %f\n", data_out[1][0]);
        //printf("m rows = %d\n", m.nr());
    }

    //for (j = 0; j < rows; j++) 
    /*
    k = 0;
    for (j = 0; j < rows; j++) {
	for (i = 0; i < resampledSize; i++) {
            data_matrix[j][i] = data_out[j][k];
        }
        if (cols - k < sampleRate)
            k++;
        else
            k = k + sampleRate;

    }
    *nc = resampledSize;
    */
/*
    for (i = 0; i < cols; i++) {
	    for (j = 0; j < rows; j++) {
            data_matrix[i][j] = data_out[j][i];
        }
    }
*/
    /*
    if (DEBUG) {
        printf("dlib matrix = %f\n", m(0,0));
        printf("dlib matrix = %f\n", m(0,1));
    }
    */

    /* Close/release resources  */
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);

    //free(data_out[0]);
    //free(data_out);

    return data_out;
}     

int
exists(int idx, int *probeids, int len) {

   int i;
   //#pragma omp parallel for shared(probeids, idx, len) private(i)
   for (i = 0; i < len; i++) 
        if (probeids[i] == idx) 
                return i;
   return -1;
}
double **
h5_read2(const char * filename, int * probeids, int len, const char * dataset_name, unsigned long *nr, unsigned long *nc)
{
    hid_t       file, dataset;         /* handles */
    hid_t       datatype, dataspace;   
    H5T_class_t clazz;                 /* data type class */
    H5T_order_t order;                 /* data order */
    size_t      size;                  /*
				                        * size of the data element	       
				                        * stored in file
				                        */
    hsize_t     dims_out[2];           /* dataset dimensions */      
    herr_t      status;                             

    double          **data_out; 
    double          **data_matrix; 
    unsigned long rows, cols, resampledSize;
    unsigned long   i, j, k, status_n, rank;


    /*  Open the file and the dataset  */
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset = H5Dopen(file, dataset_name, H5P_DEFAULT);

    /* Get datatype and dataspace handles and then query
       dataset class, order, size, rank and dimensions  */
    datatype  = H5Dget_type(dataset);     /* datatype handle */ 
    clazz     = H5Tget_class(datatype);
    if (DEBUG)
        if (clazz == H5T_FLOAT) printf("Data set has FLOAT type \n");
    order     = H5Tget_order(datatype);

    if (DEBUG)
        if (order == H5T_ORDER_LE) printf("Little endian order \n");

    size  = H5Tget_size(datatype);
    if (DEBUG)
        printf(" Data size is %d \n", size);

    dataspace = H5Dget_space(dataset);    /* dataspace handle */
    rank      = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
   // if (DEBUG)
        printf("rank %d, dimensions %lu rows x %lu cols\n", rank,
	        (unsigned long)(dims_out[0]), (unsigned long)(dims_out[1]));

    rows = dims_out[0];
    cols = dims_out[1];
    //*nc = rows;
    //*nr = cols;
    *nc = cols;
    *nr = rows;
    resampledSize = (cols)/1;
    //if (DEBUG)
        printf("resampled size = %d\n", resampledSize);

    //m.set_size(resampled_size, cols);

/*************************** BEGIN  *******************************/

/* Allocate memory for new integer array[row][col]. First
   allocate the memory for the top-level array (rows).  Make
   sure you use the sizeof a *pointer* to your data type. */

    data_out = (double**) malloc(rows*sizeof(double*));
    data_matrix = (double**) malloc(rows*sizeof(double*));
    

/* Allocate a contiguous chunk of memory for the array data values.  
   Use the sizeof the data type. */

    data_out[0] = (double*)malloc( cols*rows*sizeof(double) );
    data_matrix[0] = (double*)malloc( len*rows*sizeof(double) );

/* Set the pointers in the top-level (row) array to the
   correct memory locations in the data value chunk. */

    for (i=1; i < rows; i++) data_out[i] = data_out[0]+i*cols;
    for (i=1; i < rows; i++) data_matrix[i] = data_matrix[0]+i*len;

/************************* END *************************************/

    /* Initialize data_out array to 0 */
    for (j = 0; j < dims_out[0]; j++) {
        for (i = 0; i < dims_out[1]; i++) {
                data_out[j][i] = 0;
        }
    }
    
        /*
    for (j = 0; j < dims_out[0]; j++) {
        for (i = 0; i < resampledSize; i++) {
                data_matrix[j][i] = 0;
        }
    }
*/
    

    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		     H5P_DEFAULT, &data_out[0][0]);
    
    if (DEBUG) {
        printf("rows = %d\n", rows);
        printf("data = %f\n", data_out[0][0]);
        printf("data = %f\n", data_out[1][0]);
        //printf("m rows = %d\n", m.nr());
    }

    //for (j = 0; j < rows; j++) 
    k = 0;
    int tid, totalrows = 0;
    #pragma omp parallel for shared(data_matrix, data_out, rows, cols, totalrows) private(j, i, tid)
    for (j = 0; j < rows; j++) {
	    for (i = 0; i < cols; i++) {
            if (totalrows > len*rows) break;
            int probeIdx = exists(i+1, probeids, len);
            if (probeIdx > -1) {
                totalrows++;
                tid = omp_get_thread_num();
                printf("probeIdx = %d, tid = %d\n", probeIdx, tid);
                data_matrix[j][probeIdx] = data_out[j][i];
            } 
        }

    }
    *nc = len;
/*
    for (i = 0; i < cols; i++) {
	    for (j = 0; j < rows; j++) {
            data_matrix[i][j] = data_out[j][i];
        }
    }
*/
    /*
    if (DEBUG) {
        printf("dlib matrix = %f\n", m(0,0));
        printf("dlib matrix = %f\n", m(0,1));
    }
    */

    /* Close/release resources  */
    H5Tclose(datatype);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Fclose(file);

    free(data_out[0]);
    free(data_out);

    return data_matrix;
}     

//int main() {
//  dlib::matrix<double> textonLibrary;
// h5_read("data.h5", 1, "/X");
// 
//}
