#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "shared.h"
#include "exact-method.h"
#include "gurobi_c++.h"
#include "ilcplex/cplex.h"

int generateLpFile()
{
    int status;
    
    //Lendo instância
    status = readInstance();
    if (status) {
        printf("Erro ao ler a instância\n");
        return 1;
    }

    //Gerando a matriz de custo
    status = generateCostMatriz();
    if (status) {
        printf("Erro ao gerar a matriz de custo\n");
        freeCoordinatesInstance();
        freeCostMatriz();
        return 1;
    }

    //Abrindo(criando) o arquivo .lp
    FILE* file = fopen(targetFile, "w");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        printf("Caminho: %s\n", targetFile);
        freeCoordinatesInstance();
        freeCostMatriz();
        return 1;
    }

    //Escrevendo a FO
    fprintf(file, "Minimize\n");
    fprintf(file, "obj: z\n");

    //Escrevendo o subject
    fprintf(file, "\nSubject To\n");

    //Escrevendo a restrição: cada nó é atendido por apenas um hub
    for (int i = 1; i <= instanceEntries.nodeQuantity; i++)
    {
        for (int k = 1; k <= instanceEntries.nodeQuantity; k++)
        {
            fprintf(file, "1.0 X_%d_%d", i, k);

            if (k + 1 <= instanceEntries.nodeQuantity)
            {
                fprintf(file, " + ");
            }
        }

        fprintf(file, " = 1\n");
    }

    fprintf(file, "\n");

    //Escrevendo a restrição: um nó só pode ser hub se tiver um outro nó alocado a ele
    for (int i = 1; i <= instanceEntries.nodeQuantity; i++)
    {
        for (int k = 1; k <= instanceEntries.nodeQuantity; k++)
        {
            fprintf(file, "1.0 X_%d_%d - 1.0 X_%d_%d <= 0\n", i, k, k, k);
        }
        fprintf(file, "\n");
    }

    //Escrevendo a restrição: a solução não pode ter mais de p nó alocados como hub
    for (int k = 1; k <= instanceEntries.nodeQuantity; k++)
    {
        fprintf(file, "1.0 X_%d_%d ", k, k);

        if (k + 1 <= instanceEntries.nodeQuantity)
        {
            fprintf(file, "+ ");
        }
    }
    fprintf(file, "= %d\n\n", hubQuantity);

    //Escrevendo a restrição: raio
    for (int i = 0; i < instanceEntries.nodeQuantity; i++) {
        for (int k = 0; k < instanceEntries.nodeQuantity; k++) {
            fprintf(file, "r_%d - %lf X_%d_%d >= 0\n", k + 1, costMatrix[i][k], i + 1, k + 1);
        }
    }
    fprintf(file, "\n");
    
    for (int m = 0; m < instanceEntries.nodeQuantity; m++) {
        for (int k = 0; k <= m; k++) {
            fprintf(file, "z - r_%d - r_%d >= %lf\n", k + 1, m + 1, alpha * costMatrix[k][m]);
        }
    }
    fprintf(file, "\n");

    //Escrevendo os bounds
    fprintf(file, "BOUNDS\n");
    for (int k = 0; k < instanceEntries.nodeQuantity; k++) {
        fprintf(file, "0 <= r_%d <= inf\n", k);
    }

    //Escrevendo os bin
    fprintf(file, "\nBinary\n");

    for (int i = 1; i <= instanceEntries.nodeQuantity; i++)
    {
        for (int j = 1; j <= instanceEntries.nodeQuantity; j++)
        {
            fprintf(file, "X_%d_%d\n", i, j);
        }
    }

    //Escrevendo o fim do arquivo .lp
    fprintf(file, "\nEnd");


    //Fechando o arquivo
    fclose(file);

    //Liberando memória
    freeCoordinatesInstance();
    freeCostMatriz();

    //retornando sucesso
    return 0;
}

int cplexSolver(FILE* outputResult){

    //Caso não seja passado um arquivo a saída será impressa na tela
    if (outputResult == NULL) {
        outputResult = stdout;
    }

    // Declaração de variáveis
    int status;     // Variável para armazenar o status dos métodos
    CPXENVptr env = NULL;     // Ambiente do CPLEX
    CPXLPptr lp = NULL;       // Modelo do problema
    double objval;            // Valor da função objetivo
    int num_vars;             // Número de variáveis no modelo
    clock_t start_time, end_time; // Variáveis para medir o tempo
    double elapsed_time; // Tempo total
    int colspace = 1024; // Espaço inicial do buffer de nomes
    double lower_bound, upper_bound, obj_val, mip_gap; // saídas do sistema


    //Iniciando ambiente CPLEX
    env = CPXopenCPLEX(&status);
    if (env == NULL) {
        printf("Erro ao iniciar o ambiente do CPLEX.\n");
        return status;
    }

    //Criando o problema
    lp = CPXcreateprob(env, &status, "problema");
    if (lp == NULL) {
        printf("Erro ao criar o problema.\n");
        CPXcloseCPLEX(&env);
        return status;
    }

    //Lendo instância do arquivo .lp
    status = CPXreadcopyprob(env, lp, targetFile, NULL);
    if (status) {
        fprintf(stderr, "Erro ao carregar o arquivo LP.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return status;
    }

    //Definindo o tempo limite de execução
    status = CPXsetdblparam(env, CPXPARAM_TimeLimit, 3600.0);
    if (status) {
        printf("Erro ao configurar o limite de tempo.\n");
        CPXfreeprob(env, &lp);
        return status;
    }

    //Exibindo o log do CPLEX
    status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    if (status) {
        printf("Erro ao ativar a visualização de execução do CPLEX");
    }

    //Definindo o tipo de problema
    CPXchgprobtype(env, lp, CPXPROB_MILP);


    //Obtendo o tempo de inicio do processamento
    start_time = clock();

    //Otimizando o modelo
    status = CPXmipopt(env, lp);
    if (status) {
        fprintf(stderr, "Erro ao otimizar o modelo.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return status;
    }

    //Obtendo o tempo final
    end_time = clock();

    //Calculando o tempo total da otimização
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    //imprimindo cabeçalho
    fprintf(outputResult, "CPLEX - node %d hub %d -------\n", instanceEntries.nodeQuantity, hubQuantity);

    //obtendo o status da função objetivo
    status = CPXgetstat(env, lp);
    fprintf(outputResult, "Status da solução: %d\n", status);

    //Obtendo valor da função objetivo
    status = CPXgetobjval(env, lp, &objval);
    if (status) {
        fprintf(outputResult, "Erro ao obter o valor da função objetivo.\n");
    }
    else {
        fprintf(outputResult, "Valor da função objetivo: %.1f\n", objval);
    }

    //Obtendo o limitante inferior
    status = CPXgetbestobjval(env, lp, &lower_bound);
    if (status) {
        fprintf(outputResult, "Erro ao obter o limitante inferior.\n");
    }
    else {
        fprintf(outputResult, "Limitante inferior: %.1f\n", lower_bound);
    }

    //Obtendo o gap
    status = CPXgetmiprelgap(env, lp, &mip_gap);
    if (status) {
        fprintf(outputResult, "Erro ao obter o gap relativo.\n");
    }
    else {
        fprintf(outputResult, "Gap relativo: %.2f%%\n", mip_gap * 100);
    }

    //Obtendo o tempo total para encontrar a solução
    fprintf(outputResult, "Tempo total de execução: %.2f segundos\n", elapsed_time);

    //rodapé
    fprintf(outputResult, "-----------------------------------------\n\n");

    // Salvar a solução em um arquivo
    status = CPXsolwrite(env, lp, cplexSolFile);

    //limpando memória
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    //Retornando sucesso
    return 0;
}

int gurobiSolver(FILE* outputResult) {

    //Caso não seja passado um arquivo a saída será impressa na tela
    if (outputResult == NULL) {
        outputResult = stdout;
    }

    // Declaração de variáveis
    int status;     // Variável para armazenar o status dos métodos
    GRBenv* env = NULL;     // Ambiente
    GRBmodel* model = NULL; // Modelo
    int error;              // Código de erro
    double objval;          // Valor da função objetivo
    int numvars;            // Número de variáveis no modelo
    double* sol;            // Valores das variáveis
    char** varnames;        // Nomes das variáveis
    double lower_bound = 0;     // Limitante inferior ou superior
    double mipgap;          // Gap da solução
    clock_t start_time, end_time; //Calculo do tempo
    double elapsed_time;

    // Inicializar o ambiente
    error = GRBloadenv(&env, "gurobi.log");
    if (error) {
        printf("Error ao inicializar o ambiente");
        GRBfreeenv(env);
        return error;
    }

    // Carregar o modelo a partir do arquivo .lp
    error = GRBreadmodel(env, targetFile, &model);
    if (error) {
        printf("Error ao carregar o arquivo lp");
        GRBfreemodel(model);
        GRBfreeenv(env);
        return error;
    }

    // Definir o limite de tempo (1 hora = 3600 segundos)
    error = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, 3600);
    if (error) {
        printf("Error ao definir o tempo limite de execução");
        GRBfreemodel(model);
        GRBfreeenv(env);
        return error;
    }

    // Configurar o gap
    error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MIPGAP, 0.00001);
    if (error) {
        printf("Erro ao configurar o gap: %s\n", GRBgeterrormsg(env));
    }

    // Registrar o tempo inicial
    start_time = clock();

    // Otimizar o modelo
    error = GRBoptimize(model);
    if (error) {
        printf("Error ao otimizar o modelo");
        GRBfreemodel(model);
        GRBfreeenv(env);
        return error;
    }

    // Registrar o tempo final
    end_time = clock();
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    //imprimindo cabeçalho
    fprintf(outputResult, "GUROBI - node %d hub %d -------\n", instanceEntries.nodeQuantity, hubQuantity);

    // Obter o status da solução
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
    if (error) {
        fprintf(outputResult, "Error ao obter o status da solução");
    }
    else {
        fprintf(outputResult, "Status da solução: %d\n", status);
    }

    // Obter o valor da função objetivo
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) {
        fprintf(outputResult, "Erro ao obter o valor da função objetivo.\n");
    }
    else {
        fprintf(outputResult, "Valor da função objetivo: %.1f\n", objval);
    }

    // Obter o limitante inferior global
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJBOUND, &lower_bound);

    //Obtendo o gap
    error = GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &mipgap);
    if (error) {
        fprintf(outputResult, "Erro ao obter o gap relativo.\n");
    }
    else {
        fprintf(outputResult, "Gap relativo: %.2f%%\n", mipgap * 100);
    }

    fprintf(outputResult, "Tempo total de execução: %.2f segundos\n", elapsed_time);

    //rodapé
    fprintf(outputResult, "-----------------------------------------\n\n");

    error = GRBwrite(model, gurobiSolFile);
    if (error) {
        printf("Erro ao salvar a solução: %s\n", GRBgeterrormsg(env));
    }

    //limpando memória
    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}