#ifndef EXACT_METHOD_H
#define EXACT_METHOD_H
#define ALPHA 0.75;

char* generateLpFile(double** costMatrix, int hubQuantity, int nodeQuantity, char* fileName);
int cplexSolver(char* filePath);

#endif