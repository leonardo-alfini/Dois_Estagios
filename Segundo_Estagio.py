import sys #argv
import math
import cmath
import numpy as np
from numpy import *
from numpy.linalg import inv

informacoes = []
with open(sys.argv[1], 'r') as file:
	for line in file:
		informacoes.append(line.strip())
	file.close()

aux = []
curto = []

for i in informacoes:
	aux.append(i.replace(";", " ").split())

for j in aux:
	temp = []
	for k in j:
		temp.append(float(k))
	curto.append(temp)

informacoes = []
with open(sys.argv[2], 'r') as file:
	for line in file:
		informacoes.append(line.strip())
	file.close()

aux = []
dados = []

for i in informacoes:
	aux.append(i.replace(";", " ").split())

for j in aux:
	temp = []
	for k in j:
		temp.append(int(k))
	dados.append(temp)

ND = dados[0][0]  # Número de disjuntores da SE
NN = dados[0][1]  # Número de nós da SE
barra_curto = dados[0][2] # Barra em que ocorreu o curto
erro = 0.001     # Erro admissível

D = np.zeros((ND,2), dtype = int)
i = 0
while(i < ND):
        j = 0
        while(j < 2):
                D[i][j] = dados[i + 1][j]
                j += 1
        i += 1

shunt = np.zeros(len(dados[0])-3, dtype = int) # Vetor com os shunts
print(shunt) 
i = 3
while(i < len(dados[0])):
       shunt[i-3] = dados[0][i]
       i += 1
       
tc = 0
iteracao = 1  
while(iteracao != 0):  # Laço principal, retorna caso seja necessário alterar o erro

        A = np.zeros((NN,ND), dtype = int)   # Matriz de incidência
        
        # Formação da matriz de incidência
        i = 0
        while(i < NN):
                j = 0
                while(j < ND):
                        if(D[j][0] == i + 1):
                                A[i][j] = -1    # Convenção das correntes entrando no nó: se o disjuntor estiver em "de" --> -1                
                        if(D[j][1] == i + 1):
                                A[i][j] =  1    # Convenção das correntes entrando no nó: se o disjuntor estiver em "para" --> +1  
                        j += 1
                i += 1

        if(tc == 0):
                print("Matriz de incidência: ")
                print(A)        # Imprime matriz de incidência
        
        # Formação do vetor b para solução do sistema Ax + b = 0 (para cada fase). Convenção das correntes de injeção saindo do nó
                        
        b_A = np.zeros((NN,1),dtype=np.complex_)
        b_B = np.zeros((NN,1),dtype=np.complex_)
        b_C = np.zeros((NN,1),dtype=np.complex_)

        # Linha de b correspondente à barra de falta recebe o valor do curto-circuito

        b_A[barra_curto-1][0] += complex(curto[0][10]*math.cos(curto[0][11]*math.pi/180),curto[0][10]*math.sin(curto[0][11]*math.pi/180))
        b_B[barra_curto-1][0] += complex(curto[0][12]*math.cos(curto[0][13]*math.pi/180),curto[0][12]*math.sin(curto[0][13]*math.pi/180))
        b_C[barra_curto-1][0] += complex(curto[0][14]*math.cos(curto[0][15]*math.pi/180),curto[0][14]*math.sin(curto[0][15]*math.pi/180))

        # Injeção de corrente pelos elementos shunt
        
        i = 0
        j = 0
        while(i < len(curto)):
                if(curto[i][0] == 0):
                        b_A[shunt[j]-1][0] -= complex(curto[i][10]*math.cos(curto[i][11]*math.pi/180),curto[i][10]*math.sin(curto[i][11]*math.pi/180))
                        b_B[shunt[j]-1][0] -= complex(curto[i][12]*math.cos(curto[i][13]*math.pi/180),curto[i][12]*math.sin(curto[i][13]*math.pi/180))
                        b_C[shunt[j]-1][0] -= complex(curto[i][14]*math.cos(curto[i][15]*math.pi/180),curto[i][14]*math.sin(curto[i][15]*math.pi/180))
                        j += 1
                i += 1

        i = 2
        j = 1
        while(i < len(curto)):    # Linhas subsequentes de b recebem os valores das contribuições na mesma ordem do arquivo do ANAFAS
               if(curto[i][0] != 0):
                       b_A[j][0] -= complex(curto[i][10]*math.cos(curto[i][11]*math.pi/180),curto[i][10]*math.sin(curto[i][11]*math.pi/180))
                       b_B[j][0] -= complex(curto[i][12]*math.cos(curto[i][13]*math.pi/180),curto[i][12]*math.sin(curto[i][13]*math.pi/180))
                       b_C[j][0] -= complex(curto[i][14]*math.cos(curto[i][15]*math.pi/180),curto[i][14]*math.sin(curto[i][15]*math.pi/180))
                       j += 1
               i += 1

        linhas_usadas = np.zeros(len(A))  # Vetor para definir as linhas já usadas no escalonamento
        j = 0
        while(j < len(A[0])):
                i = 0
                while(i < len(A)):
                        if(A[i][j] != 0 and linhas_usadas[i] == 0):  # Escolhe elemento não-nulo para ser o pivô, cuja linha não tenha sido usada ainda
                                pivo = A[i][j]  # Define pivô
                                linha = i   # Marca linha do pivô
                                linhas_usadas[i] = 1   # Registra linha que será usada
                                i = len(A)   # Força saída do laço quando define pivô
                        else:
                                i += 1
                k = 0
                while(k < len(A[0])):  # Laço para tornar o pivô unitário
                      A[linha][k] /= pivo
                      k += 1
                b_A[linha][0] /= pivo
                b_B[linha][0] /= pivo
                b_C[linha][0] /= pivo

                i = 0
                while(i < len(A)):  # Operação de escalonamento entre linhas
                      if(i != linha and A[i][j] != 0):  # Condição para não haver operação da linha do pivô com ela mesma, ou com linha em que o elemento da coluna do pivô já seja nulo 
                              k = 0
                              fator = A[i][j]
                              while(k < len(A[0])):
                                      A[i][k] -= A[linha][k]*fator
                                      k += 1
                              b_A[i][0] -= b_A[linha][0]*fator
                              b_B[i][0] -= b_B[linha][0]*fator
                              b_C[i][0] -= b_C[linha][0]*fator
                      i += 1
                j += 1

        # Vetores para armazenar as correntes dos disjuntores (módulo e ângulo)
        ans_A = np.zeros((NN,2))
        ans_B = np.zeros((NN,2))
        ans_C = np.zeros((NN,2))
        
        i = 0
        while(i < NN):
                ans_A[i][1] = cmath.phase(b_A[i][0])*(180/math.pi)
                ans_A[i][0] = abs(b_A[i][0])
                ans_B[i][1] = cmath.phase(b_B[i][0])*(180/math.pi)
                ans_B[i][0] = abs(b_B[i][0])
                ans_C[i][1] = cmath.phase(b_C[i][0])*(180/math.pi)
                ans_C[i][0] = abs(b_C[i][0])
                i += 1

        i = 0
        while(i < NN):
                j = 0
                soma = 0
                while(j < ND):
                        soma += A[i][j]
                        j += 1
                if(soma == 0):          # Verifica se há linhas de zeros em A
                        if(ans_A[i][0] < erro and ans_B[i][0] < erro and ans_C[i][0] < erro): # Caso o b referente à linha de zeros seja menor que o erro, a linha é considerada redundante e é então eliminada
                                A = np.delete(A, i, 0)
                                ans_A = np.delete(ans_A, i, 0)
                                ans_B = np.delete(ans_B, i, 0)
                                ans_C = np.delete(ans_C, i, 0)
                i += 1

        i = 0
        while(i < len(A)):
                IB = 100000/(curto[i][3]*3**(1/2))  # Corrente base
                ans_A[i][0] *= IB
                ans_B[i][0] *= IB
                ans_C[i][0] *= IB
                i += 1
                        
        if(len(A) == ND):       # Se o número de linhas for igual ao de colunas, o sistema possui solução única e imprime resultados

                dj = np.zeros((ND,1), dtype = int)

                i = 0
                while(i < len(A)):
                        j = 0
                        while(j < ND):
                                if(A[i][j] == 1):
                                        dj[i][0] = j + 1
                                j += 1
                        i += 1
                
                print("Fase A:")
                i = 0
                while(i < ND):
                        print("Id",dj[i][0],"=",round(ans_A[i][0],3),"<",round(ans_A[i][1],3),"° [A]")
                        i += 1
                print("Fase B:")
                i = 0
                while(i < ND):
                        print("Id",dj[i][0],"=",round(ans_B[i][0],3),"<",round(ans_B[i][1],3),"° [A]")
                        i += 1
                print("Fase C:")
                i = 0
                while(i < ND):
                        print("Id",dj[i][0],"=",round(ans_C[i][0],3),"<",round(ans_C[i][1],3),"° [A]")
                        i += 1
                        
                with open("resultados", "w") as file:
                        file.write("Corrente nos disjuntores")
                        
                        file.write("\n" + "\n" + "Fase A:")
                        i = 0
                        while(i < ND):
                                teste = "\n" + "ID" + str(dj[i][0]) + " = " + str(round(ans_A[i][0],3))
                                file.write(teste)
                                i += 1
                                
                        file.write("\n" + "\n" + "Fase B:")
                        i = 0
                        while(i < ND):
                                teste = "\n" + "ID" + str(dj[i][0]) + " = " + str(round(ans_B[i][0],3))
                                file.write(teste)
                                i += 1
                                
                        file.write("\n" + "\n" + "Fase C:")
                        i = 0
                        while(i < ND):
                                teste = "\n" + "ID" + str(dj[i][0]) + " = " + str(round(ans_C[i][0],3))
                                file.write(teste)
                                i += 1                        
                enc = 1
        
        else:
                enc = 2
                tc += 1   # t corresponde ao número de tentativas de correção caso o sistema não tenha solução (nº linhas diferente do nº de colunas)
        
        if(enc == 2):
                if(tc == 1):
                        erro *= 10  # Tenta aumentar o erro admissível em 10x em relação ao erro original
                if(tc == 2):
                        erro *= 10  # Tenta aumentar o erro admissível em 100x em relação ao erro original
                if(tc == 3):
                        enc = 1  # Admite que o sistema não tem solução e encerra o programa
                        print("O sistema não tem solução.")
        
        if(enc == 1):
                iteracao = 0    # Sai do while principal e o programa é finalizado
