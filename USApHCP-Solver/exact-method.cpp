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
    
    status = readInstance();
    if (status) {
        printf("Erro ao ler a instância\n");
        return 1;
    }

    status = generateCostMatriz();
    if (status) {
        printf("Erro ao gerar a matriz de custo\n");
        return 1;
    }

    status = generateCostAggMatriz();
    if (status) {
        printf("Erro ao gerar a matriz agregada de custo\n");
        return 1;
    }

    // Abrir o arquivo
    if (showLogs) {
        printf("Abrindo arquivo da instância\n");
    }

    FILE* file = fopen(targetFile, "w");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        printf("Caminho: %s\n", targetFile);
        return 1;
    }

    // escrevendo a FO
    if (showLogs) {
        printf("Escrevendo arquivo lp\n");
    }

    fprintf(file, "Minimize\n");
    fprintf(file, "obj: z\n");

    // subject
    fprintf(file, "\nSubject To\n");
    
    /*fprintf(file, "- 1.0 z + ");

    double r;

    for (int i = 0; i < instanceEntries.nodeQuantity; i++)
    {
        for (int j = 0; j < instanceEntries.nodeQuantity; j++)
        {
            fprintf(file, "%lf X_%d_%d", costAggMatrix[i][j], i + 1, j + 1);

            if (i + 1 < instanceEntries.nodeQuantity || j + 1 < instanceEntries.nodeQuantity) {
                fprintf(file, " + ");
            }
        }
    }

    fprintf(file, " <= 0\n\n");*/

    // Cada nó é atendido por apenas um hub
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

    // Um nó só pode ser hub se tiver um outro nó alocado a ele
    for (int i = 1; i <= instanceEntries.nodeQuantity; i++)
    {
        for (int k = 1; k <= instanceEntries.nodeQuantity; k++)
        {
            fprintf(file, "1.0 X_%d_%d - 1.0 X_%d_%d <= 0\n", i, k, k, k);
        }
        fprintf(file, "\n");
    }

    // Não pode ter mais de p nó alocados como hub
    for (int k = 1; k <= instanceEntries.nodeQuantity; k++)
    {
        fprintf(file, "1.0 X_%d_%d ", k, k);

        if (k + 1 <= instanceEntries.nodeQuantity)
        {
            fprintf(file, "+ ");
        }
    }
    fprintf(file, "= %d\n\n", hubQuantity);

    //raio
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

    //bounds
    fprintf(file, "BOUNDS\n");
    for (int k = 0; k < instanceEntries.nodeQuantity; k++) {
        fprintf(file, "0 <= r_%d <= inf\n", k);
    }

    //bin
    fprintf(file, "\nBinary\n");

    for (int i = 1; i <= instanceEntries.nodeQuantity; i++)
    {
        for (int j = 1; j <= instanceEntries.nodeQuantity; j++)
        {
            fprintf(file, "X_%d_%d\n", i, j);
        }
    }

    fprintf(file, "\nEnd");

    fclose(file);
    
    if (showLogs) {
        printf("Arquivo lp escrito com sucesso!!!\n");
    }

    return 0;
}

int cplexSolver(){
    int status, error;

    if (toGenerateLpFile) {
        status = generateLpFile();

        if (status) {
            printf("Erro ao gerar arquivo lp\n");
            return status;
        }
    }

    // Declaração de variáveis
    CPXENVptr env = NULL;     // Ambiente do CPLEX
    CPXLPptr lp = NULL;       // Modelo do problema
    double objval;            // Valor da função objetivo
    int num_vars;             // Número de variáveis no modelo
    double* x = NULL;         // Solução das variáveis
    clock_t start_time, end_time; // Variáveis para medir o tempo
    double elapsed_time;
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
    error = CPXreadcopyprob(env, lp, targetFile, NULL);
    if (error) {
        fprintf(stderr, "Erro ao carregar o arquivo LP.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return status;
    }

    //Definindo o tempo limite de execução
    error = CPXsetdblparam(env, CPXPARAM_TimeLimit, 3600.0);
    if (error) {
        printf("Erro ao configurar o limite de tempo.\n");
    }

    //Definindo se o log do CPLEX será ou não mostrado
    if (showLogs) {
        error = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
        if (error) {
            printf("Erro ao ativar a visualização de execução do CPLEX");
        }
    }
    else {
        error = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
        if (error) {
            printf("Erro ao desativar a visualização de execução do CPLEX");
        }
    }

    //Definindo o tipo de problema
    CPXchgprobtype(env, lp, CPXPROB_MILP);


    //Obtendo o tempo de inicio do processamento
    //CPXgettime(env, &start_time);
    start_time = clock();

    //Otimizando o modelo
    error = CPXmipopt(env, lp);
    if (error) {
        fprintf(stderr, "Erro ao otimizar o modelo.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return error;
    }

    //Obtendo o tempo final
    //CPXgettime(env, &end_time);
    end_time = clock();
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    //obtendo o status da função objetivo
    status = CPXgetstat(env, lp);
    printf("Status da solução: %d\n", status);

    //Obtendo valor da função objetivo
    error = CPXgetobjval(env, lp, &objval);
    if (error) {
        fprintf(stderr, "Erro ao obter o valor da função objetivo.\n");
    }
    else {
        printf("Valor da função objetivo: %.1f\n", objval);
    }

    //Obtendo o limitante inferior
    error = CPXgetbestobjval(env, lp, &lower_bound);
    if (error) {
        fprintf(stderr, "Erro ao obter o limitante inferior.\n");
    }
    else {
        printf("Limitante inferior: %.1f\n", lower_bound);
    }

    //Obtendo o limitante superior
    error = CPXgetobjval(env, lp, &upper_bound);
    if (error) {
        fprintf(stderr, "Erro ao obter o limitante superior.\n");
    }
    else {
        printf("Limitante superior: %.1f\n", upper_bound);
    }


    //Obtendo o gap
    error = CPXgetmiprelgap(env, lp, &mip_gap);
    if (error) {
        fprintf(stderr, "Erro ao obter o gap relativo.\n");
    }
    else {
        printf("Gap relativo: %.2f%%\n", mip_gap * 100);
    }

    //Obtendo o tempo total para encontrar a solução
    printf("Tempo total de execução: %.2f segundos\n", elapsed_time);

    // Salvar a solução em um arquivo
    status = CPXsolwrite(env, lp, cplexSolFile);
    if (showLogs) {
        if (status) {
            printf("Erro ao salvar a solução: %d\n", status);
        }
        else {
            printf("Solução salva com sucesso em 'solution.sol'\n");
        }
    }

    //limpando memória
    free(x);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

int gurobiSolver() {
    int status;

    if (toGenerateLpFile) {
        status = generateLpFile();

        if (status) {
            printf("Erro ao gerar arquivo lp\n");
            return status;
        }
    }

    GRBenv* env = NULL;     // Ambiente
    GRBmodel* model = NULL; // Modelo
    int error;              // Código de erro
    double objval;          // Valor da função objetivo
    int numvars;            // Número de variáveis no modelo
    double* sol;            // Valores das variáveis
    char** varnames;        // Nomes das variáveis
    double lower_bound =0 , upper_bound=0;        // Limitante inferior ou superior
    double mipgap;          // Gap da solução
    clock_t start_time, end_time; //Calculo do tempo
    double elapsed_time;

    // Inicializar o ambiente
    error = GRBloadenv(&env, "gurobi.log");
    if (error) {
        printf("Error ao inicializar o ambiente");
        return 1;
    }

    // Carregar o modelo a partir do arquivo .lp
    error = GRBreadmodel(env, targetFile, &model);
    if (error) {
        printf("Error ao carregar o arquivo lp");
        return 1;
    }

    // Definir o limite de tempo (1 hora = 3600 segundos)
    error = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, 3600);
    if (error) {
        printf("Error ao definir o tempo limite de execução");
        return 1;
    }

    // Configurar o gap
    error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MIPGAP, 0.0001);
    if (error) {
        printf("Erro ao configurar o gap: %s\n", GRBgeterrormsg(env));
        return 1;
    }

    // Registrar o tempo inicial
    start_time = clock();

    // Otimizar o modelo
    error = GRBoptimize(model);
    if (error) {
        printf("Error ao otimizar o modelo");
        return 1;
    }

    // Registrar o tempo final
    end_time = clock();
    elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Obter o status da solução
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
    if (error) {
        printf("Error ao obter o status da solução");
        return 1;
    }
    printf("Status da solução: %d\n", status);

    // Obter o valor da função objetivo
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) {
        fprintf(stderr, "Erro ao obter o valor da função objetivo.\n");
    }
    else {
        printf("Valor da função objetivo: %.1f\n", objval);
    }

    // Obter o limitante inferior global
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJBOUND, &lower_bound);
    if (error) {
        fprintf(stderr, "Erro ao obter o limitante inferior global.\n");
    }
    else {
        printf("Limitante inferior global: %.1f\n", lower_bound);
    }

    // Obter o limite superior (igual ao valor da solução viável)
    upper_bound = objval; // Só é válido se a solução viável foi encontrada
    printf("Limitante superior global: %.1f\n", upper_bound);

    //Obtendo o gap
    error = GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &mipgap);
    if (error) {
        fprintf(stderr, "Erro ao obter o gap relativo.\n");
    }
    else {
        printf("Gap relativo: %.2f%%\n", mipgap * 100);
    }

    printf("Tempo total de execução: %.2f segundos\n", elapsed_time);

    error = GRBwrite(model, gurobiSolFile);
    if (error) {
        printf("Erro ao salvar a solução: %s\n", GRBgeterrormsg(env));
        return 1;
    }
    if (showLogs) {
        printf("Solução salva com sucesso em 'solution.sol'\n");
    }
    //limpando memória
    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}