#ifndef _BINARYUTIL_H_
#define _BINARYUTIL_H_

#include "math.h"
#include "stdlib.h"

int* numToLowMat(int number, int numNodes);

int* edgeToMat(int *edge, int numNodes);

int matToNum(int *mat, int numNodes);


#endif
