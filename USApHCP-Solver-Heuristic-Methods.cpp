#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// Para compilar com g++
#define _CRT_SECURE_NO_WARNINGS

//Constantes
#define TAM_PATH 256  
#define MAX_NODE 200
#define MAX_HUBS 20
#define ALPHA 0.75

// Structs
typedef struct
{
    double x;
    double y;
} Coordinate;

typedef struct
{
    int nodeQuantity;
    Coordinate coordinates[MAX_NODE];
} InstanceEntries;

typedef struct
{
    int node;
    double cost;
} HigherCost;

typedef struct
{
    double rk[MAX_NODE];
    int hubs[MAX_HUBS];
    int allocations[MAX_NODE];
    double FO;
    double timeGetOptimalSolution;
    double totalTime;
} Solution;

// Váriaveis definida pelo usuário
int maxTime = 0;
int hubChangePercentage;
int seed;
double T0;
double Tc;
double alpha;
int SAMax;
char instanceFile[TAM_PATH];
char solutionFile[TAM_PATH];
int amountHubs;

// Variáveis do programa
InstanceEntries instanceEntries;
double costMatrix[MAX_NODE][MAX_NODE];

void printSolution(Solution &solution) {
    FILE* arq = fopen(solutionFile, "w");

    if (arq == NULL) {  // Verifica se o arquivo abriu corretamente
        perror("Erro ao abrir o arquivo");
    }

    fprintf(arq, "FO: %.1lf\n", solution.FO);
    fprintf(arq, "TIME GET OPTMIMAL SOLUTION (seg): %lf\n", solution.timeGetOptimalSolution);
    fprintf(arq, "TOTAL TIME (seg): %.2lf\n", solution.totalTime);

    fprintf(arq, "HUBS: [ ");
    for (int i = 0; i < amountHubs; i++) {
        fprintf(arq, "%d, ", solution.hubs[i] + 1);
    }
    fprintf(arq, "]\n");

    fprintf(arq, "\nALLOCATIONS\n");

    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        fprintf(arq, "%d %d\n", i + 1, solution.allocations[i] + 1);
    }

    fclose(arq);
}

void readParameterIrace(int argc, char* argv[]) {
    char* endptr;

    // SEED
    seed = (int)strtod(argv[1], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Erro ao ler a seed\n");
    }

    // printf("seed: %d\n", seed);

    // INSTANCE
    char* separator = strchr(argv[2], '_');
    if (separator == NULL) {
        fprintf(stderr, "Erro: Formato inválido da instância.\n");
    }

    size_t pathLength = separator - argv[2];
    snprintf(instanceFile, sizeof(instanceFile), "%.*s", (int)pathLength, argv[2]);

    char* hubStr = separator + 1;
    amountHubs = (int)strtod(hubStr, &endptr);

    // printf("instancia: %s\n", instanceFile);
    // printf("qtd hubs: %d\n", amountHubs);
    if (*endptr != '\0') {
        fprintf(stderr, "Erro ao ler o número de hubs\n");
    }

    // Percorre os argumentos opcionais
    for (int i = 3; i < argc - 1; i += 2) {

        char* endptr;
        double value = strtod(argv[i + 1], &endptr);

        if (!strcmp(argv[i], "--HCP")) {
            hubChangePercentage = (int)value;
            // printf("hubChangePercentage: %d\n", hubChangePercentage);
        }
        else if (!strcmp(argv[i], "--T0")) {
            T0 = value;
            // printf("T0: %lf\n", T0);
        }
        else if (!strcmp(argv[i], "--TC")) {
            Tc = value;
            // printf("Tc: %lf\n", Tc);
        }
        else if (!strcmp(argv[i], "--ALPHA")) {
            alpha = value;
            // printf("alpha: %lf\n", alpha);
        }
        else if (!strcmp(argv[i], "--SAMAX")) {
            SAMax = (int)value;
            // printf("SAMax: %d\n", SAMax);
        }
    }
}

void readParameterMetaHeuristic(int argc, char* argv[]) {
    char* endptr;

    // SEED
    seed = (int)strtod(argv[1], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Erro ao ler a seed\n");
    }

    // printf("seed: %d\n", seed);

    // INSTANCE
    strcpy(instanceFile, argv[2]);

    // SOLUTION
    strcpy(solutionFile, argv[3]);

    // QTD HuBS
    amountHubs = (int)strtod(argv[4], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Erro ao ler a quantidade de hubs\n");
    }

    // MAX TIME
    maxTime = (int)strtod(argv[5], &endptr);
    if (*endptr != '\0') {
        fprintf(stderr, "Erro ao ler a seed\n");
    }

    // Percorre os argumentos opcionais
    for (int i = 6; i < argc - 1; i += 2) {

        char* endptr;
        double value = strtod(argv[i + 1], &endptr);

        if (!strcmp(argv[i], "--HCP")) {
            hubChangePercentage = (int)value;
            // printf("hubChangePercentage: %d\n", hubChangePercentage);
        }
        else if (!strcmp(argv[i], "--T0")) {
            T0 = value;
            // printf("T0: %lf\n", T0);
        }
        else if (!strcmp(argv[i], "--TC")) {
            Tc = value;
            // printf("Tc: %lf\n", Tc);
        }
        else if (!strcmp(argv[i], "--ALPHA")) {
            alpha = value;
            // printf("alpha: %lf\n", alpha);
        }
        else if (!strcmp(argv[i], "--SAMAX")) {
            SAMax = (int)value;
            // printf("SAMax: %d\n", SAMax);
        }
    }
}

void readInstance()
{
    // Inicialização padrão
    instanceEntries.nodeQuantity = 0;
    memset(instanceEntries.coordinates, 0, sizeof(instanceEntries.coordinates));

    // Abrir o arquivo
    FILE* file = fopen(instanceFile, "r");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        printf("Caminho: %s\n", instanceFile);
        exit(1);
    }

    // Ler a quantidade de hubs
    if (fscanf(file, "%d", &instanceEntries.nodeQuantity) != 1)
    {
        printf("Erro ao ler a quantidade de hubs no arquivo.\n");
        fclose(file);
    }

    // Ler as coordenadas
    for (int i = 0; i < instanceEntries.nodeQuantity; i++)
    {
        if (fscanf(file, "%lf %lf", &instanceEntries.coordinates[i].x, &instanceEntries.coordinates[i].y) != 2)
        {
            printf("Erro ao ler as coordenadas do hub %d.\n", i + 1);
            fclose(file);
        }
    }

    fclose(file);
}

void calcCostMatriz() {
    int size = instanceEntries.nodeQuantity; // Número de coordenadas
    Coordinate* coordinates = instanceEntries.coordinates;

    memset(costMatrix, 0, sizeof(costMatrix));

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
}

void calcRk(Solution& solution) {
    memset(solution.rk, 0, sizeof(solution.rk)); // zerando o array rk

    for (int k = 0; k < instanceEntries.nodeQuantity; k++) {
        if (costMatrix[k][solution.allocations[k]] > solution.rk[solution.allocations[k]]) {
            solution.rk[solution.allocations[k]] = costMatrix[k][solution.allocations[k]];
        }
    }
}

void calcFO(Solution& solution) {
    calcRk(solution);
    solution.FO = 0;
    double calc;

    for (int m = 0; m < instanceEntries.nodeQuantity; m++) {
        for (int k = 0; k <= m; k++) {
            calc = solution.rk[k] + solution.rk[m] + ALPHA * costMatrix[k][m];

            if (calc > solution.FO) {
                solution.FO = calc;
            }
        }
    }
}

int OrdCoastCres(const void* a, const void* b) {
    // Acessa os custos dos dois elementos
    int custoA = ((HigherCost*)a)->cost;
    int custoB = ((HigherCost*)b)->cost;

    // Ordena de forma crescente pelo custo
    if (custoA < custoB) return -1;
    if (custoA > custoB) return 1;
    return 0;
}

int descarte(int no, int *hubs, HigherCost *higherCost) {
    int hub = higherCost[hubs[0]].node;

    for (int k = 1; k < amountHubs; k++) {
        if (costMatrix[no][higherCost[hubs[k]].node] < costMatrix[no][hub]) {
            hub = higherCost[hubs[k]].node;
        }
    }

    return hub;
}

int getHubMenorCusto(int no, int *hubs){
    int hub = hubs[0];
    double cost = costMatrix[no][hub];

    for (int k = 1; k < amountHubs; k++){
        if(costMatrix[no][hubs[k]] < cost){
            hub = hubs[k];
            cost = costMatrix[no][hubs[k]];
        }
    }

    return hub;
}

void heuristic(Solution& solution) {
    HigherCost higherCost[MAX_NODE];

    // Encontrando os maiores custo para cada nó
    for (int k = 0; k < instanceEntries.nodeQuantity; k++) {
        higherCost[k].node = k;
        higherCost[k].cost = costMatrix[0][k];

        for (int i = 1; i < instanceEntries.nodeQuantity; i++) {
            if (higherCost[k].cost < costMatrix[i][k]) {
                higherCost[k].cost = costMatrix[i][k];
            }
        }
    }

    // Ordenando os maiores custos para cada destino de forma crescente
    qsort(higherCost, instanceEntries.nodeQuantity, sizeof(higherCost[0]), OrdCoastCres);

    // Criando a solução
    //Solution *solution = (Solution*) calloc(1, sizeof(Solution));
    memset(solution.allocations, -1, sizeof(solution.allocations));

    // Definindo os nós que serão hubs, com base no menor custo
    for (int h = 0; h < amountHubs; h++) {
        solution.allocations[higherCost[h].node] = higherCost[h].node;
        solution.hubs[h] = higherCost[h].node;
    }

    // Conectando os nós aos hubs
    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        if (solution.allocations[i] == -1) {
            solution.allocations[i] = getHubMenorCusto(i, solution.hubs);
        }
    }

    calcFO(solution);
}

int isHub(int node, int hubs[]) {
    for (int i = 0; i < amountHubs; i++)
    {
        if (hubs[i] == node) return 1;
    }

    return 0;
}

int changeHub(Solution &solution) {
    int oldHub = rand() % amountHubs;
    int newHub = rand() % instanceEntries.nodeQuantity;

    do
    {
        newHub = rand() % instanceEntries.nodeQuantity;
    } while (isHub(newHub, solution.hubs)); // loop até o newHub não ser um nó que está alocado como hub

   solution.hubs[oldHub] = newHub;

    // Zerando a solução
    memset(solution.allocations, -1, sizeof(solution.allocations));

    // Definindo os nós que serão hubs
    for (int h = 0; h < amountHubs; h++) {
        solution.allocations[solution.hubs[h]] = solution.hubs[h];
    }

    // Conectando os nós aos hubs
    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        if (solution.allocations[i] == -1) {
            solution.allocations[i] = getHubMenorCusto(i, solution.hubs);
        }
    }

    return 0;
}

int changeNode(Solution &solution) {
    int node;
    int newHub;

    do
    {
        node = rand() % instanceEntries.nodeQuantity;
    } while (isHub(node, solution.hubs)); // loop até o node não ser um hub

    do
    {
        newHub = rand() % amountHubs;
    } while (solution.hubs[newHub] == solution.allocations[node]);

    solution.allocations[node] = solution.hubs[newHub];

    return 0;
}

int generateNeighbor(Solution &solution) {
    if (rand() % 100 > hubChangePercentage) {
        changeNode(solution);
    }
    else {
        changeHub(solution);
    }

    calcFO(solution);
    return 0;
}

void simulatedAnnealing(Solution& bestSolution) {
    clock_t hI = clock();

    heuristic(bestSolution);

    Solution solution;
    Solution solucaoVizinha;

    do {
        memcpy(&solution, &bestSolution, sizeof(bestSolution));
        double T = T0;

        while (T > Tc) {
            for (int iterT = 0; iterT < SAMax; iterT++)
            {
                if (maxTime > 0 && ((double)(clock() - hI)) / CLOCKS_PER_SEC >= maxTime) {
                    break;
                }
                memcpy(&solucaoVizinha, &solution, sizeof(solution));
                generateNeighbor(solucaoVizinha);

                if (solucaoVizinha.FO < solution.FO) {
                    memcpy(&solution, &solucaoVizinha, sizeof(solucaoVizinha));
                    if (solution.FO < bestSolution.FO) {
                        memcpy(&bestSolution, &solution, sizeof(solution));
                        bestSolution.timeGetOptimalSolution = ((double)(clock() - hI)) / CLOCKS_PER_SEC;
                    }
                }
                else {
                    double x = rand() % 1001;
                    x = x / 1000.0;

                    if (x < exp(-(solution.FO - solucaoVizinha.FO) / T)) {
                        memcpy(&solution, &solucaoVizinha, sizeof(solucaoVizinha));
                    }
                }
            }

            T = alpha * T;
        }
    } while (maxTime > 0 && ((double)(clock() - hI)) / CLOCKS_PER_SEC < maxTime);

    bestSolution.totalTime = ((double)(clock() - hI)) / CLOCKS_PER_SEC;
}

void irace(int argc, char* argv[]) {
    Solution solution;

    readParameterIrace(argc, argv);

    srand(seed);
    readInstance();
    SAMax = log(instanceEntries.nodeQuantity * amountHubs) * SAMax;
    calcCostMatriz();
    simulatedAnnealing(solution);
    printf("%.1lf", solution.FO);
}

void metaHeuristic(int argc, char* argv[]) {
    readParameterMetaHeuristic(argc, argv);
    Solution solution;
    srand(seed);
    readInstance();
    SAMax = log(instanceEntries.nodeQuantity * amountHubs) * SAMax;
    calcCostMatriz();
    simulatedAnnealing(solution);
    printSolution(solution);
}

int main(int argc, char* argv[])
{
    // irace(argc, argv);
    metaHeuristic(argc, argv);

    return 0;
}