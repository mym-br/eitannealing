#ifndef UTILS_H
#define UTILS_H

#include "settings.h"

#define SIZETHRESH 14

#define LOGAUTO	1
#define LOGCPU	2
#define LOGCUDA	3

#define LOG(a) logMessage(a)
#define LOGN(a,b) logNumber(a, b, LOGAUTO)
#define LOGV(a,b,c) logVector(a, b, c, LOGAUTO)
#define LOGV2(a,b,c,d) logVector(a, b, c, d)
#define LOGM(a,b,c,d) logMatrix(a,b,c,d, LOGAUTO)
#define LOGM2(a,b,c,d,e) logMatrix(a,b,c,d,e)

void logInit(char * filename);

void logClose();

void logMessage(char * msg);

void logNumber(numType * data, char * msg, int mode);

void logVector(numType * data, int size, char * msg, int mode);

void logVector(int * data, int size, char * msg, int mode);

void logMatrix(numType * data, int rows, int cols, char * msg, int mode);

void logMatrix(numType * data, int rows, int * rm, char * msg, int mode);

void logMatrix(int * data, int rows, int cols, char * msg, int mode);

void logMatrix(int * data, int rows, int * rm, char * msg, int mode);

#endif /* UTILS_H */