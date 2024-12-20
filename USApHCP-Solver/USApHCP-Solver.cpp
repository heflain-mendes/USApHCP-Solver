#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include "shared.h"
#include "exact-method.h"

int toGenerateLpFile = 1;
MethodResoltion methodResoltion = GUROBI;
char sourceFile[] = "instances/inst20.txt";
char targetFile[] = "linear_programming/inst20.lp";
char gurobiSolFile[] = "gurobi_sol/inst20.sol";
char cplexSolFile[] = "cplex_sol/inst20.sol";
int hubQuantity = 10;
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
