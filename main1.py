from powerflow1 import NR 


Grid = NR() # Objeto Grid

Grid.setBarras(1, 1, 1.06, 0.0, 30.38e6 + 17.78e6 * 1j, 40e6 + -40e6 * 1j)     # Slack
Grid.setBarras(2, 3, 1.045, 0.0, 0 + 0 * 1j, 232e6 + 0 * 1j)                   # PV
Grid.setBarras(3, 3, 1.01, 0.0, 131.88e6 + 26.6e6 * 1j, 0 + 0 * 1j)            # PV
Grid.setBarras(4, 2, 1.0, 0.0, 66.92e6 + 10e6 * 1j, 0 + 0 * 1j)                # PQ
Grid.setBarras(5, 2, 1.0, 0.0, 10.64e6 + 2.24e6 * 1j, 0 + 0 * 1j)              # PQ
Grid.setBarras(6, 3, 1.07, 0.0, 15.68e6 + 10.5e6 * 1j, 0 + 0 * 1j)             # PV
Grid.setBarras(7, 2, 1.0, 0.0, 0 + 0 * 1j, 0 + 0 * 1j)                         # PQ
Grid.setBarras(8, 3, 1.09, 0.0, 0 + 0 * 1j, 0 + 0 * 1j)                        # PV
Grid.setBarras(9, 2, 1.0, 0.0, 41.3e6 + 23.24e6 * 1j, 0 + 0 * 1j)              # PQ
Grid.setBarras(10, 2, 1.0, 0.0, 12.6e6 + 8.12e6 * 1j, 0 + 0 * 1j)              # PQ
Grid.setBarras(11, 2, 1.0, 0.0, 4.9e6 + 2.52e6 * 1j, 0 + 0 * 1j)               # PQ
Grid.setBarras(12, 2, 1.0, 0.0, 8.54e6 + 2.24e6 * 1j, 0 + 0 * 1j)              # PQ
Grid.setBarras(13, 2, 1.0, 0.0, 18.9e6 + 8.12e6 * 1j, 0 + 0 * 1j)              # PQ
Grid.setBarras(14, 2, 1.0, 0.0, 20.86e6 + 7.0e6 * 1j, 0 + 0 * 1j)              # PQ

Grid.printBarras()
Grid.setSesp()

Grid.ligacoes(1, 2, impedancia = 0.01938 + 0.05914j)
Grid.ligacoes(2, 3, impedancia = 0.04699 + 0.19797j)
Grid.ligacoes(2, 4, impedancia = 0.05811 + 0.17632j)
Grid.ligacoes(1, 5, impedancia = 0.05403 + 0.22304j)
Grid.ligacoes(2, 5, impedancia = 0.05695 + 0.17388j)
Grid.ligacoes(3, 4, impedancia = 0.06701 + 0.17103j)
Grid.ligacoes(4, 5, impedancia = 0.01335 + 0.04211j)
Grid.ligacoes(5, 6, impedancia = 0.00000 + 0.25202j)
Grid.ligacoes(4, 7, impedancia = 0.00000 + 0.20912j)
Grid.ligacoes(7, 8, impedancia = 0.00000 + 0.17615j)
Grid.ligacoes(4, 9, impedancia = 0.00000 + 0.55618j)
Grid.ligacoes(9, 10, impedancia = 0.03181 + 0.08450j)
Grid.ligacoes(6, 11, impedancia = 0.09498 + 0.19890j)
Grid.ligacoes(6, 12, impedancia = 0.12291 + 0.25581j)
Grid.ligacoes(6, 13, impedancia = 0.06615 + 0.13027j)
Grid.ligacoes(9, 14, impedancia = 0.12711 + 0.27038j)
Grid.ligacoes(10, 11, impedancia = 0.08205 + 0.19207j)
Grid.ligacoes(12, 13, impedancia = 0.22092 + 0.19988j)

Grid.printLigacoes()

#Grid.Ybus()

#Grid.Sinjetada()

#Grid.setJacobiana([2], [2,3])

Grid.solveCircuito(erro = 1e-5, list_tensao=[4, 5, 7, 9, 10, 11, 12, 13, 14], list_ang=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
#Grid.Tensoes(print=True)
#Grid.Correntes(print=True)
Grid.fluxoS(printTensao = True, printCorrentes = True)
Grid.Perdas()
Grid.plotDados(tensao=True, ang=True)


