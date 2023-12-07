# Programação Concorrente e Distribuida - Rainbow Game of Life em OpenMPI

Este repositório apresenta uma implementação do [Rainbow Game of Life](https://people.kth.se/~gunnarj/LIFE/WLIF/wlcframes.html).

Primeiramente, cada processo é associado a um grupo de linhas do tabuleiro completo.

A cada iteração, processos vizinhos trocam informações sobre suas linhas na fronteira por meio de comunicação MPI. Após isso, eles atualizam o valor de cada uma de suas linhas, agora que tem conhecimento daquelas na fronteira.

Ao final da execução, executa MPI_Reduce centrado no processo 0 para calcular o número de células vivas em todo o tabuleiro.
