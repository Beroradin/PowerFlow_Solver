import pandas as pd
import numpy as np
import time 

comeco = time.time()
# Configurações iniciais
erro_max = 0.000001  # Critério de convergência
K_max = 1000  # Número máximo de iterações
contador = 0
erro = 1

# Carregamento dos arquivos
admitancia = r"C:\Users\João Pedro\Desktop\Simulação_SEP\matriz admitância.xlsx"
barras = r"C:\Users\João Pedro\Desktop\Simulação_SEP\Barras.xlsx"
impedancias = r"C:\Users\João Pedro\Desktop\SEP\impedância.xlsx"
# Carregando a matriz de admitância
matriz_admt = pd.read_excel(admitancia, index_col=1, header=0).iloc[:, 1:].copy()
for i in range(len(matriz_admt)):
    for j in range(len(matriz_admt)):
        valor = str(matriz_admt.iloc[i, j]).replace(",", ".").replace("i", "j")
        matriz_admt.iloc[i, j] = complex(valor)

# Carregando os dados das barras
tipo_barras = pd.read_excel(barras, index_col=0, header=0)
impedancias = pd.read_excel(impedancias, index_col=0, header=0)
vetor_tensao = [complex(str(valor)) for valor in tipo_barras["VOLTAGE MAGNITUDE"]]
vetor_pot_ativa = [
    (float(tipo_barras.iloc[i]["GENERATOR(MW)"]) - float(tipo_barras.iloc[i]["LOAD(MW)"])) / 100
    for i in range(len(tipo_barras))
]
vetor_pot_reativa = [
    (float(tipo_barras.iloc[i]["GENERATOR(MVAR)"]) - float(tipo_barras.iloc[i]["LOAD(MVAR)"])) / 100
    for i in range(len(tipo_barras))
]
carga_reativa = [float(tipo_barras.iloc[i]["LOAD(MVAR)"]) / 100 for i in range(len(tipo_barras))]  # QL(k)

vetor_resistencia = [
    float(impedancias.iloc[i]["RESISTÊNCIA"]) for i in range(len(impedancias))
]
vetor_reatancia = [
    float(impedancias.iloc[i]["REATÂNCIA"]) for i in range(len(impedancias))
]

# Iterações
while (erro > erro_max) and (contador < K_max):
    vetor_tensao_antiga = vetor_tensao.copy()
    contador += 1
    erro = 0

    for k in range(len(tipo_barras)):
        if tipo_barras.index[k] == 1:  # Barra Slack
            continue  # A tensão da barra slack não muda

        # Calcula a soma de YV
        YV = sum(matriz_admt.iloc[k, n] * vetor_tensao[n] for n in range(len(tipo_barras)) if k != n)

        if tipo_barras.index[k] == 0:  # Barra PQ
            vetor_tensao[k] = (1 / matriz_admt.iloc[k, k]) * ((vetor_pot_ativa[k] + 1j * vetor_pot_reativa[k]) / vetor_tensao[k].conjugate() - YV)

        elif tipo_barras.index[k] == 2:  # Barra PV
            # Calcula a potência reativa
            Q_calc = -np.imag(vetor_tensao[k].conjugate() * (YV + matriz_admt.iloc[k, k] * vetor_tensao[k]))
            # Calcula Q_liq como a diferença entre Q_calc e QL
            Q_liq = Q_calc - carga_reativa[k]

            # Atualiza a tensão
            vetor_tensao[k] = (1 / matriz_admt.iloc[k, k]) * ((vetor_pot_ativa[k] + 1j * Q_liq) / vetor_tensao[k].conjugate() - YV)
            vetor_tensao[k] = abs(vetor_tensao_antiga[k]) * (vetor_tensao[k] / abs(vetor_tensao[k]))

        # Atualiza o erro
        erro = max(erro, abs(vetor_tensao[k] - vetor_tensao_antiga[k]))

print(f"Tempo de execução: {time.time() - comeco:.2f} segundos\n")
# Resultados
print(f"Convergiu após {contador} iterações com erro: {erro:.8f}")
for i, tensao in enumerate(vetor_tensao):
    modulo = abs(tensao)  # Calcula o módulo da tensão
    angulo = np.angle(tensao)  # Calcula o ângulo da tensão (em radianos)
    angulo = np.degrees(angulo)
    print("Tensão na barra {:d}: {:.3f} ∠ {:.3f}°".format(i + 1, modulo, angulo))

P_gerada_nova = np.zeros(len(tipo_barras))
Q_gerada_nova = np.zeros(len(tipo_barras))
P_consumida = [
    float(tipo_barras.iloc[i]["LOAD(MW)"]) for i in range(len(tipo_barras))
]
Q_consumida = [
    float(tipo_barras.iloc[i]["LOAD(MVAR)"]) for i in range(len(tipo_barras))
]

for i in range(len(tipo_barras)):
    I_injetada = np.dot(matriz_admt.iloc[i, :], vetor_tensao)
    S_injetada = vetor_tensao[i] * np.conj(I_injetada)  
    P_consumida[i] = P_consumida[i] / 100
    Q_consumida[i] = Q_consumida[i] / 100
    P_gerada_nova[i] = np.real(S_injetada) + P_consumida[i]
    Q_gerada_nova[i] = -np.imag(S_injetada) + Q_consumida[i]

print("\nPotências geradas após a convergência:")

for i in range(len(tipo_barras)):
    print(f"Barra {i + 1}: |P = {P_gerada_nova[i] * 100:.2f} MW| Q = {Q_gerada_nova[i] * 100:.2f} MVar|")

potencia_ativa_linhas = np.zeros(impedancias.shape[0])
potencia_reativa_linhas = np.zeros(impedancias.shape[0])    
#print(impedancias.columns)
vetor_de = [
    int(impedancias.iloc[i]["DE"]) - 1 for i in range(len(impedancias))
]
vetor_para = [
    int(impedancias.iloc[i]["PARA"]) - 1  for i in range(len(impedancias))   
]
#print(vetor_de, vetor_para)
for i in range(impedancias.shape[0]):
    de_Barra = vetor_de[i]
    para_Barra = vetor_para[i]
    R = vetor_resistencia[i]
    X = vetor_reatancia[i]
    Z = R + 1j * X
    V_de = vetor_tensao[de_Barra]
    V_para = vetor_tensao[para_Barra]
    I = (V_de - V_para) / Z
    potencia_ativa_linhas[i] = np.real(V_de * np.conj(I))
    potencia_reativa_linhas[i] = np.imag(V_de * np.conj(I))

print("\nPotências ativas e reativas das linhas:")
for i in range(impedancias.shape[0]):
    print(f"Barra {int(impedancias.iloc[i,0])}-{int(impedancias.iloc[i,1])}: |P = {potencia_ativa_linhas[i]*100:.2f} MW| Q = {potencia_reativa_linhas[i]*100:.2f} MVar|")

perdas_ativas_linhas = np.zeros(impedancias.shape[0])
perdas_reativas_linhas = np.zeros(impedancias.shape[0])
perdas_totais_ativas = 0
perdas_totais_reativa = 0
for i in range(impedancias.shape[0]):
    de_Barra = vetor_de[i]
    para_Barra = vetor_para[i]
    R = vetor_resistencia[i]
    X = vetor_reatancia[i]
    Z = R + 1j * X
    V_de = vetor_tensao[de_Barra]
    V_para = vetor_tensao[para_Barra]
    corrente = (V_de - V_para) / Z
    perdas_ativas_linhas[i] = np.real(corrente) ** 2 * R + np.imag(corrente) ** 2 * R
    perdas_reativas_linhas[i] = np.real(corrente) ** 2 * X + np.imag(corrente) ** 2 * X
    perdas_totais_ativas += perdas_ativas_linhas[i]
    perdas_totais_reativa += perdas_reativas_linhas[i]

print("\nPerdas de potência ativa (em MW) e reativa (em MVar) em cada linha:")
for i in range(impedancias.shape[0]):
    print(f"Barra {int(impedancias.iloc[i, 0])}-{int(impedancias.iloc[i, 1])}: |P = {perdas_ativas_linhas[i] * 100:.4f} MW| Q = {perdas_reativas_linhas[i] * 100:.4f} MVar|")
print("\nPerdas totais de potência ativa (em MW) e reativa (em MVar):")
print(f"P = {perdas_totais_ativas * 100:.4f} MW | Q = {perdas_totais_reativa * 100:.4f} MVar")