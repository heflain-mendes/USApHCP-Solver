#ifndef SHARED_H
#define SHERED_H
#define MAX_PATH_LENGTH 100

typedef struct
{
    double x;
    double y;
} Coordinate;

typedef struct
{
    int nodeQuantity;
    Coordinate* coordinates;
} InstanceEntries;

extern const char* filesNames[];
extern const int hubsQuantity[];

extern char sourceFile[];
extern char targetFile[];
extern char gurobiSolFile[];
extern char cplexSolFile[];
extern int hubQuantity;
extern double alpha;
extern InstanceEntries instanceEntries;
extern double** costMatrix;
extern double** costAggMatrix;

int readInstance();
int freeCoordinatesInstance();
int generateCostMatriz();
int freeCostMatriz();

#endif