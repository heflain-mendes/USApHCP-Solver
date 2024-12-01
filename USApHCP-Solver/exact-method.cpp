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
    fprintf(file, "- 1.0 z + ");

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

    fprintf(file, " <= 0\n\n");

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


    fprintf(file, "Binary\n");

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
    int status = generateLpFile();

    if (status) {
        printf("Erro ao gerar arquivo lp\n");
        return status;
    }

    // Declaração de variáveis
    CPXENVptr env = NULL;     // Ambiente do CPLEX
    CPXLPptr lp = NULL;       // Modelo do problema
    double objval;            // Valor da função objetivo
    int num_vars;             // Número de variáveis no modelo
    double* x = NULL;         // Solução das variáveis
    double start_time, end_time; // Variáveis para medir o tempo
    char** colname;
    char* colnamebuf = NULL;
    int surplus;
    int colspace = 1024; // Espaço inicial do buffer de nomes

    env = CPXopenCPLEX(&status);
    if (env == NULL) {
        printf("Erro ao iniciar o ambiente do CPLEX.\n");
        return status;
    }

    lp = CPXcreateprob(env, &status, "problema");
    if (lp == NULL) {
        printf("Erro ao criar o problema.\n");
        CPXcloseCPLEX(&env);
        return status;
    }

    status = CPXreadcopyprob(env, lp, targetFile, NULL);
    if (status) {
        fprintf(stderr, "Erro ao carregar o arquivo LP.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return status;
    }

    status = CPXsetdblparam(env, CPXPARAM_TimeLimit, 3600.0);
    if (status) {
        printf("Erro ao configurar o limite de tempo.\n");
    }

    if (showLogs) {
        status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
        if (status) {
            printf("Erro ao ativar a visualização de execução do CPLEX");
        }
    }
    else {
        status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
        if (status) {
            printf("Erro ao desativar a visualização de execução do CPLEX");
        }
    }

    CPXchgprobtype(env, lp, CPXPROB_LP);

    CPXgettime(env, &start_time);

    status = CPXprimopt(env, lp);
    if (status) {
        fprintf(stderr, "Erro ao otimizar o modelo.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return status;
    }

    CPXgettime(env, &end_time);

    status = CPXgetobjval(env, lp, &objval);
    if (status) {
        fprintf(stderr, "Erro ao obter o valor da função objetivo.\n");
    }
    else {
        printf("Valor da função objetivo: %f\n", objval);
    }

    // Recuperar os nomes das variáveis
    //num_vars = CPXgetnumcols(env, lp);
    //colname = (char**)malloc(num_vars * sizeof(char*));
    //colnamebuf = (char*)malloc(colspace * sizeof(char));

    //if (colname == NULL || colnamebuf == NULL) {
    //    fprintf(stderr, "Erro de alocação de memória.\n");
    //    free(colname);
    //    free(colnamebuf);
    //    CPXfreeprob(env, &lp);
    //    CPXcloseCPLEX(&env);
    //    return 1;
    //}

    //status = CPXgetcolname(env, lp, colname, colnamebuf, colspace, &surplus, 0, num_vars - 1);
    //if (status == CPXERR_NEGATIVE_SURPLUS) {
    //    // Realoca buffer caso necessário
    //    colspace = -surplus;
    //    free(colnamebuf);
    //    colnamebuf = (char*)malloc(colspace * sizeof(char));
    //    status = CPXgetcolname(env, lp, colname, colnamebuf, colspace, &surplus, 0, num_vars - 1);
    //}

    //x = (double*)malloc(num_vars * sizeof(double));
    //if (x == NULL) {
    //    fprintf(stderr, "Erro ao alocar memória.\n");
    //    CPXfreeprob(env, &lp);
    //    CPXcloseCPLEX(&env);
    //    return -1;
    //}

    //status = CPXgetx(env, lp, x, 0, num_vars - 1);
    //if (status) {
    //    fprintf(stderr, "Erro ao obter os valores das variáveis.\n");
    //}
    //else {
    //    printf("Solução:\n");
    //    for (int i = 0; i < num_vars; i++) {
    //        printf("%s = %f\n", colname[i], x[i]);
    //    }
    //}

    double lower_bound, upper_bound, obj_val, mip_gap;

    status = CPXgetbestobjval(env, lp, &lower_bound);
    if (status) {
        fprintf(stderr, "Erro ao obter o limitante inferior.\n");
    }
    else {
        printf("Limitante inferior: %f\n", lower_bound);
    }

    status = CPXgetmipobjval(env, lp, &upper_bound);
    if (status) {
        fprintf(stderr, "Erro ao obter o limitante superior.\n");
    }
    else {
        printf("Limitante superior: %f\n", upper_bound);
    }

    status = CPXgetmiprelgap(env, lp, &mip_gap);
    if (status) {
        fprintf(stderr, "Erro ao obter o gap relativo.\n");
    }
    else {
        printf("Gap relativo: %f%%\n", mip_gap * 100);
    }

    printf("Tempo total de execução: %f segundos\n", end_time - start_time);

    free(x);
    //free(colname);
    free(colnamebuf);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

int gurobiSolver() {
    int status = generateLpFile();

    if (status) {
        printf("Erro ao gerar arquivo lp\n");
        return status;
    }

    GRBenv* env = NULL;     // Ambiente
    GRBmodel* model = NULL; // Modelo
    int error;              // Código de erro
    double objval;          // Valor da função objetivo
    int numvars;            // Número de variáveis no modelo
    double* sol;            // Valores das variáveis
    char** varnames;        // Nomes das variáveis
    double objbound;        // Limitante inferior ou superior
    double mipgap;          // Gap da solução

    // Inicializar o ambiente
    error = GRBloadenv(&env, "gurobi.log");
    if (error) {
        printf("Error: %s\n", GRBgeterrormsg(env));
        return 1;
    }

    // Carregar o modelo a partir do arquivo .lp
    error = GRBreadmodel(env, targetFile, &model);
    if (error) {
        printf("Error: %s\n", GRBgeterrormsg(env));
        return 1;
    }

    // Definir o limite de tempo (1 hora = 3600 segundos)
    error = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, 3);
    if (error) {
        printf("Error: %s\n", GRBgeterrormsg(env));
        return 1;
    }

    // Registrar o tempo inicial
    clock_t start_time = clock();

    // Otimizar o modelo
    error = GRBoptimize(model);
    if (error) {
        printf("Error: %s\n", GRBgeterrormsg(env));
        return 1;
    }

    // Registrar o tempo final
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Obter o status da solução
    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &status);
    if (error) {
        printf("Error: %s\n", GRBgeterrormsg(env));
        return 1;
    }

    printf("Time to solve: %.2f seconds\n", elapsed_time);

    if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
        // Obter o valor da função objetivo
        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
        if (!error) {
            printf("Optimal objective value: %f\n", objval);
        }

        // Obter os limitantes e o gap
        GRBgetdblattr(model, GRB_DBL_ATTR_OBJBOUND, &objbound);
        GRBgetdblattr(model, GRB_DBL_ATTR_MIPGAP, &mipgap);
        printf("Objective bound: %f\n", objbound);
        printf("MIP gap: %f\n", mipgap);

        // Obter o número de variáveis no modelo
        error = GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &numvars);
        if (error) {
            printf("Error: %s\n", GRBgeterrormsg(env));
            return 1;
        }

        // Alocar memória para valores e nomes das variáveis
        sol = (double*)malloc(numvars * sizeof(double));
        varnames = (char**)malloc(numvars * sizeof(char*));
        for (int i = 0; i < numvars; i++) {
            varnames[i] = (char*)malloc(100 * sizeof(char)); // Buffer para o nome
        }

        if (sol == NULL || varnames == NULL) {
            printf("Memory allocation error.\n");
            return 1;
        }

        // Obter os valores das variáveis
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, numvars, sol);
        if (error) {
            printf("Error: %s\n", GRBgeterrormsg(env));
            free(sol);
            for (int i = 0; i < numvars; i++) free(varnames[i]);
            free(varnames);
            return 1;
        }

        // Obter os nomes das variáveis
        error = GRBgetstrattrarray(model, GRB_STR_ATTR_VARNAME, 0, numvars, varnames);
        if (error) {
            printf("Error: %s\n", GRBgeterrormsg(env));
            free(sol);
            for (int i = 0; i < numvars; i++) free(varnames[i]);
            free(varnames);
            return 1;
        }

        // Imprimir os valores das variáveis
        /*printf("Variable values:\n");
        for (int i = 0; i < numvars; i++) {
            printf("%s = %f\n", varnames[i], sol[i]);
        }*/

        if (varnames) {
            //for (int i = 0; i < numvars; i++) {
            //    if (varnames[i]) free(varnames[i]); // Libera cada ponteiro individual
            //}
            free(varnames); // Libera o array de ponteiros
        }

        free(sol);
    }
    else {
        printf("No optimal solution found. Status: %d\n", status);
    }

    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}
