#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "shared.h"
#include "exact-method.h"
#include "gurobi_c++.h"
#include "ilcplex/cplex.h"

char* generateLpFile(double** costMatrix, int hubQuantity, int nodeQuantity, char* filePath)
{
    // Abrir o arquivo
    FILE* file = fopen(filePath, "w");
    if (file == NULL)
    {
        perror("Erro ao abrir o arquivo");
        printf("Caminho: %s\n", filePath);
        return NULL;
    }

    // escrevendo a FO
    fprintf(file, "Minimize\n");
    fprintf(file, "obj: z\n");

    double alpha = ALPHA;

    // subject
    fprintf(file, "\nSubject To\n");
    fprintf(file, "fo:  - 1.0 z + ");

    double r;

    for (int i = 0; i < hubQuantity; i++)
    {
        for (int j = 0; j < hubQuantity; j++)
        {
            for (int m = 0; m < hubQuantity; m++)
            {
                for (int k = 0; k < hubQuantity; k++)
                {
                    r = costMatrix[i][k] + alpha * costMatrix[k][m];
                    
                    fprintf(file, "%lf X_%d_%d + %lf X_%d_%d",
                        r, i + 1, k + 1, costMatrix[m][j], j + 1, m + 1);

                    if (i + 1 < hubQuantity || j + 1 < hubQuantity || k + 1 < hubQuantity || m + 1 < hubQuantity)
                    {
                        fprintf(file, " + \\\n");
                    }
                }
            }
        }
    }

    fprintf(file, " <= 0\n\n");

    // Cada nó é atendido por apenas um hub
    for (int i = 1; i <= hubQuantity; i++)
    {
        fprintf(file, "customer_service_%d: ", i);
        for (int k = 1; k <= hubQuantity; k++)
        {
            fprintf(file, "1.0 X_%d_%d", i, k);

            if (k + 1 <= hubQuantity)
            {
                fprintf(file, " + ");
            }
        }

        fprintf(file, " = 1\n");
    }

    fprintf(file, "\n");

    // Um nó só pode ser hub se tiver um outro nó alocado a ele
    for (int i = 1; i <= hubQuantity; i++)
    {
        for (int k = 1; k <= hubQuantity; k++)
        {
            fprintf(file, "facility_restriction_%d_%d: 1.0 X_%d_%d - 1.0 X_%d_%d <= 0\n", i, k, i, k, k, k);
        }
        fprintf(file, "\n");
    }

    // Não pode ter mais de p nó alocados como hub
    fprintf(file, "facility_limit: ");
    for (int k = 1; k <= hubQuantity; k++)
    {
        fprintf(file, "1.0 X_%d_%d ", k, k);

        if (k + 1 <= hubQuantity)
        {
            fprintf(file, "+ ");
        }
    }
    fprintf(file, "= %d\n\n", nodeQuantity);


    fprintf(file, "Binary\n");

    for (int i = 1; i <= hubQuantity; i++)
    {
        for (int j = 1; j <= hubQuantity; j++)
        {
            fprintf(file, "X_%d_%d\n", i, j);
        }
    }

    fprintf(file, "\nEnd");
}

int cplexSolver(char *filePath){
    // Declaração de variáveis
    CPXENVptr env = NULL;     // Ambiente do CPLEX
    CPXLPptr lp = NULL;       // Modelo do problema
    int status;               // Para verificar erros
    double objval;            // Valor da função objetivo
    int num_vars;             // Número de variáveis no modelo
    double* x = NULL;         // Solução das variáveis
    double start_time, end_time; // Variáveis para medir o tempo

    env = CPXopenCPLEX(&status);
    if (env == NULL) {
        printf("Erro ao iniciar o ambiente do CPLEX.\n");
        return status;
    }

    CPXgettime(env, &start_time);

    lp = CPXcreateprob(env, &status, "problema");
    if (lp == NULL) {
        printf("Erro ao criar o problema.\n");
        CPXcloseCPLEX(&env);
        return status;
    }

    status = CPXreadcopyprob(env, lp, filePath, NULL);
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

    status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
    if (status) {
        printf("Erro ao ativar a visualização de execução do CPLEX");
    }

    CPXchgprobtype(env, lp, CPXPROB_LP);

    status = CPXdualopt(env, lp);
    if (status) {
        fprintf(stderr, "Erro ao otimizar o modelo.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return status;
    }

    status = CPXgetobjval(env, lp, &objval);
    if (status) {
        fprintf(stderr, "Erro ao obter o valor da função objetivo.\n");
    }
    else {
        printf("Valor da função objetivo: %f\n", objval);
    }

    num_vars = CPXgetnumcols(env, lp);
    x = (double*)malloc(num_vars * sizeof(double));
    if (x == NULL) {
        fprintf(stderr, "Erro ao alocar memória.\n");
        CPXfreeprob(env, &lp);
        CPXcloseCPLEX(&env);
        return -1;
    }

    status = CPXgetx(env, lp, x, 0, num_vars - 1);
    if (status) {
        fprintf(stderr, "Erro ao obter os valores das variáveis.\n");
    }
    else {
        printf("Solução:\n");
        for (int i = 0; i < num_vars; i++) {
            printf("x[%d] = %f\n", i, x[i]);
        }
    }

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

    status = CPXgetobjval(env, lp, &obj_val);
    if (status) {
        fprintf(stderr, "Erro ao obter o valor da solução.\n");
    }
    else {
        printf("Valor ótimo da função objetivo: %f\n", obj_val);
    }

    CPXgettime(env, &end_time);
    printf("Tempo total de execução: %f segundos\n", end_time - start_time);

    free(x);
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    return 0;
}

//int gurobiSolver(char* filePath) {
//    GRBenv* env = NULL;    // Ambiente
//    GRBmodel* model = NULL; // Modelo
//    int error;             // Código de erro
//    int optimstatus;       // Status da otimização
//    double objval;         // Valor da função objetivo
//
//    // Inicializar o ambiente
//    error = GRBloadenv(&env, "gurobi.log");
//    if (error) {
//        printf("Error: %s\n", GRBgeterrormsg(env));
//        return 1;
//    }
//
//    // Carregar o modelo a partir do arquivo .lp
//    error = GRBreadmodel(env, filePath, &model);
//    if (error) {
//        printf("Error: %s\n", GRBgeterrormsg(env));
//        return 1;
//    }
//
//    // Definir o limite de tempo (1 hora = 3600 segundos)
//    error = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, 3600.0);
//    if (error) {
//        printf("Error: %s\n", GRBgeterrormsg(env));
//        return 1;
//    }
//
//    // Otimizar o modelo
//    error = GRBoptimize(model);
//    if (error) {
//        printf("Error: %s\n", GRBgeterrormsg(env));
//        return 1;
//    }
//
//    // Obter o status da solução
//    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
//    if (error) {
//        printf("Error: %s\n", GRBgeterrormsg(env));
//        return 1;
//    }
//
//    // Verificar se a solução foi encontrada
//    if (optimstatus == GRB_OPTIMAL) {
//        // Obter o valor da função objetivo
//        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
//        if (error) {
//            printf("Optimal objective value: %f\n", objval);
//        }
//        else {
//            printf("Error: %s\n", GRBgeterrormsg(env));
//
//        }
//    }
//    else if (optimstatus == GRB_TIME_LIMIT) {
//        printf("Time limit reached. Partial solution may be available.\n");
//    }
//    else {
//        printf("No optimal solution found. Status: %d\n", optimstatus);
//    }
//
//    GRBfreemodel(model);
//    GRBfreeenv(env);
//
//    return 0;
//}