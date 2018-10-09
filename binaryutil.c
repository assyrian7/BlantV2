#include "binaryutil.h"

int* numToLowMat(int number, int numNodes){
    int *mat = (int*)malloc(numNodes * numNodes * sizeof(int));
    int bitlength = (numNodes * (numNodes - 1)) / 2;
    int count = 0;
    for(int i = 0; i < numNodes * numNodes; i++){
        mat[i] = 0;
    }
    for(int i = 1; i < numNodes; i++){
        for(int j = 0; j < i; j++){
            if((number & (1 << (bitlength - count - 1))) != 0){
                mat[numNodes * i + j] = 1;
            }
            count++;
        }
    }
    return mat;
}

int* edgeToMat(int *edge, int numNodes){
    int number = 0;
    int bitlength = (numNodes * (numNodes - 1)) / 2;
    int *mat = (int*)malloc(numNodes * numNodes * sizeof(int));
    for(int i = 0; i < numNodes * numNodes; i++){
        mat[i] = 0;
    }
    int i = 0;
    while(i < (sizeof(edge) / sizeof(int) / 2)){
        int i = *edge;
        int j = *(edge + 1);
        if(i == 0 && j == 0) break;
        mat[numNodes * i + j] = 1;
        edge++;
    }
    return mat;
}

int matToNum(int *mat, int numNodes){
    int number = 0;
    int bitlength = (numNodes * (numNodes - 1)) / 2;
    int count = 0;
    for(int i = 1; i < numNodes; i++){
        for(int j = 0; j < i; j++){
            if(mat[numNodes * i + j] == 1) number |= (1 << (bitlength - count - 1));
        }
        count++;
    }
    return number;
}
