#ifndef __H5_HYPERSLAB_H
#define __H5_HYPERSLAB_H

double **
hyperslabread(const char * FILE, const char *DATASETNAME,
                int count_row, int count_col,
                int offsetrow, int offsetcol);


#endif
