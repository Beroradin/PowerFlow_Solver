#Programa para cálculo do fluxo de potência em um sistema elétrico de potência.

#Passos para execução do programa:
#01 - Especificar os parâmetros das barras e das linhas de transmissão e suas condições iniciais.
#02 - Determinar as injeções líquidas de potência e a matriz admitância.
#03 - Detereminar os resíduos de potências das barras PQ.
#04 - Determinar o resíduo de potência ativa das barra PV.
#05 - Determinar os elementos de J1, J2, J3 e J4 da matriz jacobiana.
#06 - Resolver a equação linear da matriz jacobiana (determinando o fluxo de potência).
#07 - Indicar novas tensões e ângulos.
#08 - Repetir os passos 3 a 7 até que os resíduos sejam menores que a tolerância estabelecida.
#09 - Determinar a injeção líquida de potência na barra slack.
#10 - Determinar o fluxo de potência nas linhas e as perdas.

#bibliotecas
import cmath as cmt
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# Classe para o cálculo do fluxo de potência
class NR:
    def __init__(self): # Dunder
        self.__Sbase = 100e6       # 100 MVA Base
        self.__dados = dict()      # dicionário para armazenar os dados das barras

        self.__Sesp = dict()
        self.__Sbarras = dict()
        self.__ligacoes = dict()   # Muitos dicionários kkk
        self.__ybus = list()       # Finalmente uma lista
        self.__V = dict()
        self.__I = dict()
        self.__fluxS = dict()
        self.__Perdas = 0

        self.__tensaoPlot = dict()
        self.__angPlot = dict()
        self.__list_tensao = list()
        self.__list_ang = list()

        self.nPV = int()
        self.nPQ = int()

        #Atributos das submatrizes jacobianas e da jacobiana
        self.__J1 = list()
        self.__J2 = list()
        self.__J3 = list()
        self.__J4 = list()
        self.__Jacobiana = list()

    def setBarras(self, barra, code, tensao, ang, carga, geracao):
        # code 1 -> slack / 2 -> PQ / 3 -> PV
        self.__dados[barra] = {'code': code, 'tensao': tensao,
                                'ang': mt.radians(ang), 'carga': (carga/self.__Sbase),'geracao': (geracao/self.__Sbase)} #transformando o ângulo em radianos
        self.__tensaoPlot[barra] = [tensao]
        self.__angPlot[barra] = [ang]
    
    def printBarras(self):
        # printando as barras
        print('\n===================================DADOS DAS BARRAS===================================')
        print('Sbase = ', self.__Sbase, 'VA')
        for i in self.__dados: print(self.__dados[i])
        print('======================================================================')

    def setSesp(self):
        # Método utilizado para calcular a potência espeficificada em cada barra.
        for i in self.__dados:
            if self.__dados[i]['code'] == 2:   # Barra PQ
                self.__Sesp[i] = {'Pesp' : np.real(self.__dados.get(i)['geracao'] - self.__dados.get(i)['carga']),
                                  'Qesp' : float(np.imag(self.__dados.get(i)['geracao'] - self.__dados.get(i)['carga']))}
            elif self.__dados[i]['code'] == 3: # Barra PV
                self.__Sesp[i] = {'Pesp' : np.real(self.__dados.get(i)['geracao']),
                                  'Qesp' : float(np.imag(self.__dados.get(i)['geracao']))}
                
        print('\n======================================POTÊNCIAS ESPECIFICADAS======================================')
        print(self.__Sesp, 'pu')
        print('====================================================================================================')
                
    def ligacoes(self, b1, b2, impedancia=None, admitancia=None):
        # Método utilizado para especificar as ligações entre as barras.
        if impedancia is not None:
            admitancia = 1 / impedancia
        
        if admitancia is not None:
            impedancia= 1 / admitancia

        self.__ligacoes[(b1, b2)] = { 
            'impedancia' : impedancia,
            'admitancia' : admitancia
        }
        
    def printLigacoes(self):
        print('\n===================================LIGAÇÕES===================================') 
        for i in self.__ligacoes: print('Ligação:', i, '\t', self.__ligacoes[i]) # Utilizando a tupla
        print('===========================================================================')

    def __printYbus(self):
        print('\n===================================MATRIZ ADMITÂNCIA===================================')
        for i in self.__ybus: print(i)
        print('=========================================================================================')


    def __printYbus(self):
        print('\n===================================MATRIZ ADMITÂNCIA===================================')
        for i in self.__ybus: print(i)
        print('=========================================================================================')

    def Ybus(self):
        # Método utilizado para calcular a matriz admitância do sistema.
        self.__ybus = np.zeros((len(self.__dados), len(self.__dados)), dtype=complex)

        for i in range(len(self.__ybus)):
            lin = []
            for j in range(len(self.__ybus)):
                if i == j:
                    lin.append(0)
                else:
                    if self.__ligacoes.__contains__(tuple([i + 1, j + 1])):
                        lin.append(-self.__ligacoes.get(tuple([i + 1, j + 1]))['admitancia'])
                    elif self.__ligacoes.__contains__(tuple([j + 1, i + 1])):
                        lin.append(-self.__ligacoes.get(tuple([j + 1, i + 1]))['admitancia'])
                    else:
                        lin.append(0)
            for j in range(len(self.__ybus)):
                if i == j:
                    lin[j] = -1 * sum(lin)
            self.__ybus[i] = lin
        self.__printYbus()

        for i in self.__dados:
            if self.__dados.get(i)['code'] == 2:
                self.nPQ += 1
            elif self.__dados.get(i)['code'] == 3:
                self.nPV += 1

    def Sinjetada(self):

        self.__Sinjetada = dict()
        self.__deltaPeQ = []
        self.__ResiduoP = []
        self.__ResiduoQ = []

        # Método utilizado para calcular as potências injetadas nas barras.
        for i in self.__dados:
            soma1 = []
            soma2 = []
            if self.__dados[i]['code'] != 1: # Não pode ser barra slack
                for j in self.__dados:
                    soma1.append(
                        abs(self.__ybus[i - 1][j - 1]) * abs(self.__dados.get(i)['tensao']) * 
                        abs(self.__dados.get(j)['tensao']) * 
                        mt.cos(np.angle(self.__ybus[i - 1][j - 1]) - self.__dados.get(i)['ang'] + self.__dados.get(j)['ang'])
                    )
                    

                    soma2.append(
                        -abs(self.__ybus[i - 1][j - 1]) * abs(self.__dados.get(i)['tensao']) * 
                        abs(self.__dados.get(j)['tensao']) * 
                        mt.sin(np.angle(self.__ybus[i - 1][j - 1]) - self.__dados.get(i)['ang'] + self.__dados.get(j)['ang']) * 1j
                    )
            
                self.__ResiduoP.append(
                    np.real(
                        self.__Sesp.get(i)['Pesp'] - sum(soma1)
                    )
                )
                if self.__dados[i]['code'] == 2:
                    self.__ResiduoQ.append(
                        np.imag(
                            self.__Sesp.get(i)['Qesp'] * 1j - sum(soma2)
                        )
                    )
        
        for i in range(len(self.__ResiduoP)):
            self.__deltaPeQ.append(self.__ResiduoP[i])
        for i in range(len(self.__ResiduoQ)):
            self.__deltaPeQ.append(self.__ResiduoQ[i])
        print('===================================RESÍDUOS===================================')
        for i in self.__deltaPeQ: 
            print(i)

    def __calcularJ1(self, list_ang, nPQ, nPV):
        # Método utilizado para calcular a submatriz J1.
        self.__J1 = np.ones((nPQ + nPV, nPQ + nPV))

        mainDiagonal = []
        outDiagonal = []

        for i in list_ang:
            soma = []
            for j in range(1, len(self.__dados) + 1, 1):
                if i != j:
                    soma.append(
                        abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(i)['tensao'])
                        * abs(self.__dados.get(j)['tensao'])
                        * cmt.sin(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )
                                 )
                    )
            mainDiagonal.append(sum(soma)) # Atribuindo o valor da soma na diagonal principal

        for i in list_ang:
            for j in list_ang:
                if i != j:
                    outDiagonal.append(
                        -abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(i)['tensao'])
                        * abs(self.__dados.get(j)['tensao'])
                        * cmt.sin(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )    
                        )
                    )
        m = 0
        for i in range(len(list_ang)):
            for j in range(len(list_ang)):
                if i == j:
                    self.__J1[i][j] = np.real(mainDiagonal[j])
                else:
                    self.__J1[i][j] = np.real(outDiagonal[m])
                    m += 1
        
        #print('\nJ1 = \n', self.__J1)

        return self.__J1
    
    def __calcularJ2(self, list_tensao, list_ang, nPQ, nPV):
        # Método utilizado para calcular a submatriz J2.
        self.__J2 = np.ones((nPQ + nPV, nPQ))

        mainDiagonal = []
        outDiagonal = []
     
        for i in list_ang:
            soma = []
            a = 0; # Auxiliar para percorrer a lista de tensões
            for j in range(1, len(self.__dados) + 1, 1):
                if i != j:
                    soma.append(
                        abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(j)['tensao'])
                        * cmt.cos(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )
                                 )
                    )
            a = (2 * abs(self.__dados.get(i)['tensao']) * abs(self.__ybus[i - 1][i - 1])
                 * cmt.cos(cmt.phase(self.__ybus[i - 1][i - 1])))
            mainDiagonal.append(a + sum(soma)) # Atribuindo o valor da soma na diagonal principal

        for i in list_ang:
            for j in list_tensao:
                if i != j:
                    outDiagonal.append(
                        abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(i)['tensao'])
                        * cmt.cos(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )    
                        )
                    )
        
        m = 0
        for i in range(nPQ + nPV):
            k = nPV
            for j in range(nPQ):
                if i < nPV:
                    self.__J2[i][j] = np.real(outDiagonal[m])
                elif i >= nPV:
                    if i - nPV == j:
                        self.__J2[i][j] = np.real(mainDiagonal[j + nPV])
                        k += 1
                    else:
                        self.__J2[i][j] = np.real(outDiagonal[m])
                        m += 1
        #print('\nk = ', k)
        #print('\nJ2 = \n', self.__J2)

        return self.__J2
    
    def __calcularJ3(self, list_tensao, list_ang, nPQ, nPV):
        # Método utilizado para calcular a submatriz J3.
        self.__J3 = np.ones((nPQ,nPQ + nPV))

        mainDiagonal = []
        outDiagonal = []
     
        for i in list_ang:
            soma = []
            for j in range(1, len(self.__dados) + 1, 1):
                if i != j:
                    soma.append(
                        abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(i)['tensao'])
                        * abs(self.__dados.get(j)['tensao'])
                        * cmt.cos(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )
                                 )
                    )
            mainDiagonal.append(sum(soma)) # Atribuindo o valor da soma na diagonal principal

        for i in list_ang:
            for j in list_tensao:
                if i != j:
                    outDiagonal.append(
                        -abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(i)['tensao'])
                        * abs(self.__dados.get(j)['tensao'])
                        * cmt.cos(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )    
                        )
                    )
        
        m = 0
        for i in range(nPQ):
            for j in range(nPQ + nPV):
                if j < nPV:
                    self.__J3[i][j] = np.real(outDiagonal[m])
                    m += 1
                elif j >= nPV:
                    if j - nPV == i:
                        self.__J3[i][j] = np.real(mainDiagonal[i + nPV])
                    else:
                        self.__J3[i][j] = np.real(outDiagonal[m])
                        m += 1
        #print('\nJ3 = \n', self.__J3)

        return self.__J3
    
    def __calcularJ4(self, list_tensao, list_ang, nPQ, nPV):
        # Método utilizado para calcular a submatriz J4.
        self.__J4 = np.ones((nPQ, nPQ))

        mainDiagonal = []
        outDiagonal = []
     
        for i in list_ang:
            soma = []
            a = 0; 
            for j in range(1, len(self.__dados) + 1, 1):
                if i != j:
                    soma.append(
                        abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(j)['tensao'])
                        * cmt.sin(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )
                                 )
                    )
            a = (2 * abs(self.__dados.get(i)['tensao']) * abs(self.__ybus[i - 1][i - 1])
                 * cmt.sin(cmt.phase(self.__ybus[i - 1][i - 1])))
            mainDiagonal.append(-a - sum(soma)) # Atribuindo o valor da soma na diagonal principal

        for i in list_ang:
            for j in list_tensao:
                if i != j:
                    outDiagonal.append(
                        -abs(self.__ybus[i - 1][j - 1])
                        * abs(self.__dados.get(i)['tensao'])
                        * cmt.sin(cmt.phase(self.__ybus[i - 1][j - 1]- 
                                            self.__dados.get(i)['ang'] + 
                                            self.__dados.get(j)['ang']
                                           )    
                        )
                    )
        
        m = 0
        for i in range(nPQ):
            for j in range(nPQ):
                if i == j:
                    self.__J4[i][j] = np.real(mainDiagonal[j + nPV])	
                else:
                    self.__J4[i][j] = np.real(outDiagonal[m])
                    m += 1
        #print('\nJ4 = \n', self.__J4)

        return self.__J4
    
    def setJacobiana(self, list_tensao, list_ang):

        self.__Jacobiana = []
        self.__list_tensao = list_tensao
        self.__list_ang = list_ang
        nXn = len(list_tensao) + len(list_ang)

        J1 = self.__calcularJ1(list_ang, self.nPQ, self.nPV)
        J2 = self.__calcularJ2(list_tensao, list_ang, self.nPQ, self.nPV)
        J3 = self.__calcularJ3(list_tensao, list_ang, self.nPQ, self.nPV)
        J4 = self.__calcularJ4(list_tensao, list_ang, self.nPQ, self.nPV)

        self.__Jacobiana = np.zeros((nXn, nXn)) # Matriz Jacobiana é quadrada

        self.__Jacobiana = np.block([
            [J1, J2],
            [J3, J4]
        ])
        
        print('===================================JACOBIANA===================================')
        for i in self.__Jacobiana: print(i)
        print('=========================================================================================')
        

    def linearSolver(self):
        # Método utilizado para calcular os resultados do sistema linear da equação de potência.

        self.__x = []

        # Converter a Jacobiana para matriz esparsa antes de resolver
        J_sparse = sp.csc_matrix(self.__Jacobiana)

        # Resolver J * x = -F usando solver para matrizes esparsas
        self.__x = spla.spsolve(J_sparse, self.__deltaPeQ)

        # Verificar se a solução está correta
        correto = np.allclose(J_sparse @ self.__x, self.__deltaPeQ)
        print('\nO Sistema linear foi calculado corretamente?', correto)

        ang = []
        tensao = []
        for i in range(len(self.__x)):
            if i < (self.nPQ + self.nPV):
                ang.append(self.__x[i])
            else:
                tensao.append(self.__x[i])     # Primeiro salva os ângulos e depois as tensões

        m = 0
        sorted_barras = sorted(self.__dados.keys())     # Ordenando as barras evito pegar valores nulos com get(i) no índice zero
        for barra_num in sorted_barras:
            if self.__dados[barra_num]['code'] != 1:    # Barras PQ e PV precisam atualizar o ang
                self.__dados[barra_num]['ang'] += float(np.real(ang[m]))
                self.__angPlot[barra_num].append(self.__dados[barra_num]['ang'])  
                m += 1
        m = 0
        for barra_num in sorted_barras:
            if self.__dados[barra_num]['code'] == 2:    # Barras PQ precisam atualizar a tensão
                self.__dados[barra_num]['tensao'] += float(np.real(tensao[m]))
                self.__tensaoPlot[barra_num].append(self.__dados[barra_num]['tensao'])
                m += 1
    
    def novaInjecao(self):
        # Método utilizado para calcular o novo valor de injeção de potência aparente nas barras PV e slack.

        self.__Sbarras = dict()

        # Método utilizado para calcular as potências injetadas nas barras.
        for i in self.__dados:
            soma1 = []
            soma2 = []
            if self.__dados[i]['code'] != 2: # Não pode ser barra PQ
                for j in self.__dados:
                    soma1.append(
                        abs(self.__ybus[i - 1][j - 1]) * abs(self.__dados.get(i)['tensao']) * 
                        abs(self.__dados.get(j)['tensao']) * 
                        mt.cos(np.angle(self.__ybus[i - 1][j - 1]) - self.__dados.get(i)['ang'] + self.__dados.get(j)['ang'])
                    )

                    soma2.append(
                        -abs(self.__ybus[i - 1][j - 1]) * abs(self.__dados.get(i)['tensao']) * 
                        abs(self.__dados.get(j)['tensao']) * 
                        mt.sin(np.angle(self.__ybus[i - 1][j - 1]) - self.__dados.get(i)['ang'] + self.__dados.get(j)['ang']) * 1j
                    )
            if self.__dados[i]['code'] == 1:
                self.__Sbarras[i] = {'P' : np.real(sum(soma1)), 'Q' : np.imag(sum(soma2))} 
            elif self.__dados[i]['code'] == 3:
                self.__Sbarras[i] = {'Q' : np.imag(sum(soma2))}  # Somente a Reativa para PV
            
        for i in self.__dados:
            if self.__dados[i]['code'] == 1:
                self.__dados[i]['geracao'] = self.__Sbarras.get(i)['P'] + self.__Sbarras.get(i)['Q'] * 1j
            elif self.__dados[i]['code'] == 3:
                self.__dados[i]['geracao'] = np.real(self.__Sbarras.get(i)['geracao'] + self.__Sbarras.get(i)['Q'] * 1j)

    
    def solveCircuito(self, erro=None, iteracoes=None, list_tensao=None, list_ang=None):
        # Método utilizado para resolver o circuito, calculando o fluxo de potência.
        self.__list_tensao = list_tensao
        self.__list_ang = list_ang
        self.count = 1

        self.Ybus()
        self.Sinjetada()
        self.setJacobiana(list_tensao = self.__list_tensao, list_ang = self.__list_ang)
        self.linearSolver()

        if iteracoes is None and erro is not None:
            pEq = list(map(abs, self.__deltaPeQ))
            teste = list(map(lambda m: True if (m < erro) else False, pEq))
            stop = teste.count(False)
            while True:
                self.Sinjetada()
                self.setJacobiana(list_tensao = self.__list_tensao, list_ang = self.__list_ang)
                self.linearSolver()
                self.count += 1
                pEq = list(map(abs, self.__deltaPeQ))
                teste = list(map(lambda m: True if (m < erro) else False, pEq))
                stop = teste.count(False)
                if stop == 0:
                    break
        elif iteracoes is not None and erro is None:
            while self.count < iteracoes:
                self.Sinjetada()
                self.setJacobiana(list_tensao = self.__list_tensao, list_ang = self.__list_ang)
                self.linearSolver()
                self.count += 1

        if iteracoes is not None:
            print('\n==================================NÚMERO DE ITERAÇÕES==================================')
            print('Número de iterações = ', self.count)
        elif erro is not None:
            print('\nConvergiu para um erro de', erro, '.')
            print('Convergiu em', self.count, 'iterações.')

    def __printTensao(self):
        print('\n===================================TENSÕES===================================')
        for i in self.__V: print('Barra: \t', i, '\tTensão = \t', self.__V.get(i), ' p.u.')
        print('\n=============================================================================')


    def Tensoes(self, print=None):
        # Método utilizado para calcular as tensões em cada barra em p.u.
        self.__V = dict()
        for i in self.__dados:
            self.__V[i] = cmt.rect(self.__dados.get(i)['tensao'], self.__dados.get(i)['ang'])
        if print:
            self.__printTensao()

    def __printCorrentes(self):
        print('\n===================================CORRENTES===================================')
        for i in self.__I: print('Corrente: \t', i, '\tValor = \t', self.__I.get(i), ' p.u.')
        print('\n=============================================================================')
    
    def Correntes(self, print=None):
        # Método utilizado para calcular os valores das correntes em cada linha considerando o ângulo das tensôes, seus valores são calculados como o somatório de todas as correntes da barra sob análise.

        self.__I = dict()
        self.Tensoes(print=None)
        
        for i in self.__dados:
            for j in self.__dados:
                if i == j:
                    continue  # Ignora a diagonal (será calculada depois)
                if (i, j) in self.__ligacoes or (j, i) in self.__ligacoes:
                    # Calcula I_ij = (V_i - V_j) * Y_ij
                    Y_ij = self.__ybus[i-1][j-1]
                    V_i = self.__V.get(i)
                    V_j = self.__V.get(j)
                    self.__I[(i, j)] = (V_i - V_j) * Y_ij
                    
        # Corrente total injetada na barra (I_ii = Σ I_ij para j ≠ i)
        for i in self.__dados:
            I_total = 0
            for j in self.__dados:
                if i != j and ((i, j) in self.__I or (j, i) in self.__I):
                    I_total += self.__I.get((i, j), 0)
            self.__I[(i, i)] = I_total
        
        if print:
            self.__printCorrentes()

    def fluxoS(self, printTensao=None, printCorrentes=None):
        # Método utilizado para calcular o fluxo de potências nas linhas e barras do sistema.

        self.__fluxS = dict()
        self.Tensoes(print=printTensao)
        self.Correntes(print=printCorrentes)

        # Calcula fluxo de potência nas linhas
        for (i, j) in self.__I:
            if i == j:
                continue  # Ignora correntes injetadas
            V_i = self.__V.get(i)
            I_ij = self.__I[(i, j)]
            self.__fluxS[(i, j)] = V_i * np.conjugate(I_ij)

        print('\n===================================FLUXO DE POTÊNCIA===================================')
        for key in self.__fluxS:
            print(f'Ligação: {key}\tValor: {self.__fluxS[key]:.6f} p.u.')

        # Atualiza geração nas barras Slack e PV (se necessário)
        for i in self.__dados:
            if self.__dados[i]['code'] == 1:  # Barra Slack
                self.__dados[i]['geracao'] = self.__fluxS.get((i, i), 0)

    def Perdas(self):
        # Método utilizado para calcular as perdas do sistema.

        self.__Perdas = 0
        processed_pairs = set()  # Para evitar dupla contagem

        for (i, j) in self.__fluxS:
            if i == j or (j, i) in processed_pairs:
                continue
            S_ij = self.__fluxS.get((i, j), 0)
            S_ji = self.__fluxS.get((j, i), 0)
            perda_linha = S_ij + S_ji
            self.__Perdas += perda_linha
            processed_pairs.add((i, j))

        print('\n===================================PERDAS===================================')
        print(f'Perdas Totais: {self.__Perdas:.6f} p.u.')

    def __plotTensao(self):
        # Método utilizado para plotar a tensão das barras após o algoritmo de Newton-Raphson.

        x = self.count
        barras = []
        y = []
        for i in self.__dados:
            if self.__dados.get(i)['code'] == 2:
                barras.append(i)
        for i in barras:
            y.append(self.__tensaoPlot.get(i))
        for i in range(len(barras)):
            plt.subplot(len(barras), 1, i + 1)
            plt.plot(range(x + 1), y[i])
            plt.title('Variação da tensão na barra' + str(barras[i]) + ' X número de iterações')
            plt.xlabel('Iterações')
            plt.ylabel('Tensão na barra' + str(barras[i]) + ' p.u.')
            plt.grid(True)
        plt.tight_layout()
        plt.show()

    def __plotAng(self):
        # Método utilizado para plotar a tensão das barras após o algoritmo de Newton-Raphson.

        x = self.count
        barras = []
        y = []
        for i in self.__dados:
            if self.__dados.get(i)['code'] != 1:
                barras.append(i)
        for i in barras:
            y.append(self.__angPlot.get(i))
        for i in range(len(barras)):
            plt.subplot(len(barras), 1, i + 1)
            plt.plot(range(x + 1), y[i])
            plt.title('Variação do ângulo na barra' + str(barras[i]) + ' X número de iterações')
            plt.xlabel('Iterações')
            plt.ylabel('Ângulo na barra' + str(barras[i]) + ' p.u.')
            plt.grid(True)
        plt.tight_layout()
        plt.show()
    
    def plotDados(self, tensao=None, ang=None):
        # Método utilizado para plotar os dados das barras após o algoritmo de Newton-Raphson.

        if tensao:
            self.__plotTensao()
        if ang:
            self.__plotAng()



