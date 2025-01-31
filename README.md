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
   
## Metódo de utilização
  # Gauss-Seidel
    
    Para que seja possível utilizar e testar o código de Gauss-Seidel é necessário fazer o download dos arquivos "Matriz Admitância", "Barras" e "impedância", logo após o download, busque no explorador de arquivos o link do diretório em que cada uma das pastas do excel se encontra e substitua o link do seu computador no código, para o link da tabela da matriz admitância substitua na variavel "admitancia" no código, para o link das barras substitua na variável "barras", para o link da admitância substitua na variável "impedâncias".        
 ![image](https://github.com/user-attachments/assets/1f1960e1-3035-4d85-a48e-249b65472559)
 ![Captura de Tela (62)](https://github.com/user-attachments/assets/01c60a0f-fd32-430b-84bf-204fc05e9ade)

