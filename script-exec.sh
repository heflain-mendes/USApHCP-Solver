#!/bin/bash

# Definindo as duas listas
instancias=("inst10" "inst20" "inst25" "inst40" "inst50" "inst100" "inst200")
hubs=("2" "2" "3" "4" "4" "10" "10")
path_instances="Instances"
path_solutions="Solutions"
tempo=300

# Verificando se as listas têm o mesmo tamanho
if [ ${#instancias[@]} -ne ${#hubs[@]} ]; then
    echo "Erro: As listas têm tamanhos diferentes!"
    exit 1
fi

# Percorrendo as listas usando índices
for i in "${!instancias[@]}"
do
    # Gerando uma seed aleatória entre 1 e 1000
    seed=$RANDOM

    # Definindo o caminho do arquivo de instância
    instance_file="${path_instances}/${instancias[$i]}.txt"

    # Chamando o programa externo e passando os parâmetros
    for ((j = 1; j <= 3; j++)); do
        # Definindo o nome do arquivo de solução com base na iteração
        solution_file="${path_solutions}/${instancias[$i]}_${j}.txt"

        # Executando o programa
        ./USApHCP-Solver "$seed" "$instance_file" "$solution_file" "${hubs[$i]}" "$tempo" --HCP 47 --ALPHA 0.9945 --SAMAX 726 --T0 3916.4429 --TC 0.0099
    done
done