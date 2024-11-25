#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shared.h"

InstanceEntries readInstance(char* instancePath)
{
    InstanceEntries instance;
    // Inicialização padrão
    instance.hubQuantity = 0;
    instance.coordinates = NULL;

    // Abrir o arquivo
    FILE* file = fopen(instancePath, "r");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        printf("Caminho: %s\n", instancePath);
        return instance;
    }

    printf("Arquivo aberto: %s\n", instancePath);

    // Ler a quantidade de hubs
    if (fscanf(file, "%d", &instance.hubQuantity) != 1)
    {
        printf("Erro ao ler a quantidade de hubs no arquivo.\n");
        fclose(file);
        return instance;
    }

    printf("A quantidade de hub: %u\n", instance.hubQuantity);

    // Alocar memória para as coordenadas
    instance.coordinates = (Coordinate*)malloc(instance.hubQuantity * sizeof(Coordinate));
    if (instance.coordinates == NULL)
    {
        printf("Erro ao alocar memória para coordenadas.\n");
        fclose(file);
        instance.hubQuantity = 0;
        return instance;
    }

    printf("Memória alocada: %zu bytes\n", instance.hubQuantity * sizeof(Coordinate));

    // Ler as coordenadas
    for (int i = 0; i < instance.hubQuantity; i++)
    {
        if (fscanf(file, "%lf %lf", &instance.coordinates[i].x, &instance.coordinates[i].y) != 2)
        {
            printf("Erro ao ler as coordenadas do hub %d.\n", i + 1);
            free(instance.coordinates);
            fclose(file);
            return instance;
        }
    }

    printf("Coordenadas lidas!!!\n");

    fclose(file);
    return instance;
}

double** generateCostMatrix(InstanceEntries instanceEntries) {
    int size = instanceEntries.hubQuantity; // Número de coordenadas
    Coordinate* coordinates = instanceEntries.coordinates;

    // Alocar a matriz de custos (matriz dinâmica)
    double** costMatrix = (double**)malloc(size * sizeof(double*));
    if (costMatrix == NULL) {
        printf("Erro ao alocar a matriz de custos.\n");
        return NULL;
    }

    for (int i = 0; i < size; i++) {
        costMatrix[i] = (double*)malloc(size * sizeof(double));
        if (costMatrix[i] == NULL) {
            printf("Erro ao alocar a linha %d da matriz de custos.\n", i);
            // Liberar linhas previamente alocadas antes de retornar
            for (int j = 0; j < i; j++) {
                free(costMatrix[j]);
            }
            free(costMatrix);
            return NULL;
        }
    }

    // Preencher a matriz com os custos (distâncias)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                costMatrix[i][j] = 0.0; // Distância de um ponto para ele mesmo é 0
            }
            else {
                double dx = coordinates[i].x - coordinates[j].x;
                double dy = coordinates[i].y - coordinates[j].y;
                costMatrix[i][j] = sqrt(dx * dx + dy * dy); // Distância Euclidiana
            }
        }
    }

    return costMatrix;
}