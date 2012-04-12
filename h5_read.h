#ifndef __H5_READ_H
#define __H5_READ_H

//#include "dlib/matrix.h"

//using namespace dlib;

//matrix<double>
//h5_read(char * filename, int sampleRate, char * dataset_name);

double **
h5_read(const char * filename, int sampleRate, const char * dataset_name, unsigned long *nr, unsigned long *nc);

double **
h5_read2(const char * filename, int *probeids, int probeidlen, const char * dataset_name, unsigned long *nr, unsigned long *nc);

#endif
