#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include "shared.h"
#include "exact-method.h"

#define TAM_PATH 100
#define QTD_INSTANCIAS 7

const char* filesNames[] = { "inst10", "inst20", "inst25", "inst40", "inst50", "inst100", "inst200" };
const int hubsQuantity[] = { 2, 2, 3, 4, 4, 10, 10 };

char sourceFile[TAM_PATH];
char targetFile[TAM_PATH];
char gurobiSolFile[TAM_PATH];
char cplexSolFile[TAM_PATH];
int hubQuantity;
double alpha = 0.75;
InstanceEntries instanceEntries;
double** costMatrix;
double** costAggMatrix;

void escreverPath(const char* fileName) {
    //source
    memset(sourceFile, 0, TAM_PATH);
    strcat(sourceFile, "instances/");
    strcat(sourceFile, fileName);
    strcat(sourceFile, ".txt");

    //target(lp)
    memset(targetFile, 0, TAM_PATH);
    strcat(targetFile, "linear_programming/");
    strcat(targetFile, fileName);
    strcat(targetFile, ".lp");

    //gurobi sol
    memset(gurobiSolFile, 0, TAM_PATH);
    strcat(gurobiSolFile, "gurobi_sol/");
    strcat(gurobiSolFile, fileName);
    strcat(gurobiSolFile, ".sol");

    //cplex sol
    memset(cplexSolFile, 0, TAM_PATH);
    strcat(cplexSolFile, "cplex_sol/");
    strcat(cplexSolFile, fileName);
    strcat(cplexSolFile, ".sol");
}


int main()
{
    FILE* outputFile = fopen("output.txt", "w");
    if (outputFile == NULL) {
        perror("Erro ao abrir o arquivo");
        return 1;
    }

    for (int i = 0; i < QTD_INSTANCIAS; i++) {
        escreverPath(filesNames[i]);
        hubQuantity = hubsQuantity[i];
        generateLpFile();
        cplexSolver(outputFile);
        gurobiSolver(outputFile);
    }

    fclose(outputFile);
    return 0;
}
