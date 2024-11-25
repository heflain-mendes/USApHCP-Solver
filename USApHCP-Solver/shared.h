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
    int hubQuantity;
    Coordinate* coordinates;
} InstanceEntries;

// instances must be in the instances folder
InstanceEntries readInstance(char* instance);
double** generateCostMatrix(InstanceEntries instanceEntries);

#endif