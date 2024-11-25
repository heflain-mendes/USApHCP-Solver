#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdlib.h>
#include <stdio.h>
#include <cstring>

#include "shared.h"
#include "exact-method.h"

typedef struct struct_input
{
    char* sourceFile = (char*)calloc(100, sizeof(char));
    char* targetFile = (char*)calloc(100, sizeof(char));
    int nodeQuantity = 0;
} Input;

Input* createInput() {
    Input* input = (Input*)malloc(sizeof(Input));
    if (input == NULL) {
        fprintf(stderr, "Erro ao alocar memória para a struct.\n");
        exit(1);
    }

    input->sourceFile = (char*)calloc(100, sizeof(char));
    if (input->sourceFile == NULL) {
        fprintf(stderr, "Erro ao alocar memória para sourceFile.\n");
        free(input);
        exit(1);
    }

    input->targetFile = (char*)calloc(100, sizeof(char));
    if (input->targetFile == NULL) {
        fprintf(stderr, "Erro ao alocar memória para targetFile.\n");
        free(input->sourceFile);
        free(input);
        exit(1);
    }

    input->nodeQuantity = 0;
    return input;
}

void freeInput(Input* input) {
    if (input != NULL) {
        free(input->sourceFile);
        free(input->targetFile);
        free(input);
    }
}

void menu(Input* input) {
    printf("Preencha os dados abaixo:\n");

    printf("Digite o nome do arquivo de origem: ");
    fgets(input->sourceFile, 100, stdin);
    input->sourceFile[strcspn(input->sourceFile, "\n")] = '\0'; // Remove o '\n' de fgets

    printf("Digite o nome do arquivo de destino: ");
    fgets(input->targetFile, 100, stdin);
    input->targetFile[strcspn(input->targetFile, "\n")] = '\0'; // Remove o '\n' de fgets

    printf("Digite a quantidade de nós: ");
    scanf("%d", &(input->nodeQuantity));
    getchar(); // Limpa o buffer de entrada
}

int main()
{
    Input* input = createInput();
    menu(input);

    InstanceEntries instance = readInstance(input->sourceFile);

    if (instance.coordinates == NULL)
    {
        printf("Falha ao carregar a instância, presente em: %s.\n", input->sourceFile);
        return 1;
    }

    for (int i = 0; i < instance.hubQuantity; i++)
    {
        printf("Hub %d: (%lf, %lf)\n", i + 1, instance.coordinates[i].x, instance.coordinates[i].y);
    }

    double** costMatrix = generateCostMatrix(instance);

    for (int i = 0; i < instance.hubQuantity; i++)
    {
        for (int j = 0; j < instance.hubQuantity; j++)
        {
            printf("%.2lf ", costMatrix[i][j]);
        }
        printf("\n");
    }

    generateLpFile(costMatrix, instance.hubQuantity, input->nodeQuantity, input->targetFile);
    cplexSolver(input->targetFile);

    for (int i = 0; i < instance.hubQuantity; i++)
    {
        free(costMatrix[i]);
    }

    freeInput(input);
    free(costMatrix);
    free(instance.coordinates);
    return 0;
}