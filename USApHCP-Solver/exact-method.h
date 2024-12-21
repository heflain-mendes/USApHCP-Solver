#ifndef EXACT_METHOD_H
#define EXACT_METHOD_H

extern int toGenerateLpFile;

int generateLpFile();
int cplexSolver(FILE* outputResult);
int gurobiSolver(FILE* outputResult);

#endif