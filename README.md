# PowerFlow_Solver

## Descrição do Projeto
**PowerFlow_Solver** é um programa desenvolvido em Python com o objetivo de resolver problemas de fluxo de potência em sistemas de energia elétrica. É implementado dois algoritmos clássicos de análise de fluxo de potência: Gauss-Seidel e Newton-Raphson. O programa tem como finalidade calcular as perdas de um sistemas IEEE com 14 barras. Esse programa foi desenvolvido como atividade avaliativa da componente curricular de Sistemas de Potência do curso de Engenharia Elétrica da Universidade Federal do Oeste da Bahia.

## Equipe
- João Pedro Soares Raccolto
- Matheus Pereira Alves

## Funcionalidades

- Resolução de fluxo de potência:
  - Cálculo das tensões nodais e ângulos de fase.
  - Determinação de potências ativa e reativa em barras específicas.
  - Cálculo das perdas durante a transmissão.
- Algoritmos disponíveis:
  - Gauss-Seidel:
    - Método iterativo simples e eficiente para sistemas pequenos.
    - Requer menos memória, mas pode convergir mais lentamente.
  - Newton-Raphson:
    - Método iterativo rápido e robusto para sistemas maiores.
    - Utiliza a matriz Jacobiana para convergência mais eficiente.