#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include "shared.h"
#include "exact-method.h"

int toGenerateLpFile = 0;
MethodResoltion methodResoltion = GUROBI;
char sourceFile[] = "instances/inst200.txt";
char targetFile[] = "interactive-optimizers/inst200.lp";
int hubQuantity = 40;
int showLogs = 1;
double alpha = 0.75;
InstanceEntries instanceEntries;
double** costMatrix;
double** costAggMatrix;

int main()
{
    switch (methodResoltion)
    {
    case CPLEX:
        cplexSolver();
        break;
    case GUROBI:
        gurobiSolver();
        break;
    }

    freeCoordinatesInstance();
    freeCostMatriz();
    freeCostAggMatriz();
    return 0;
}


void printCostMatriz() {
    int size = instanceEntries.nodeQuantity;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%8.2lf", costMatrix[i][j]);

            if (j + 1 < size) {
                printf("\t");
            }
            else {
                printf("\n");
            }
        }
    }
}
