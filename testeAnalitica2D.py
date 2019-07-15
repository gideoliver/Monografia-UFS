# -*- coding: utf-8 -*-
"""
Created on Sat Jun 01 19:11:16 2019

@author: gideo
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import minhasFuncoes as mf
from matplotlib import cm

#### #################        DADOS DO PRPBLEMA                              ##############

k = 1.0e-3               # coeficiente de transferencia de calpor
r = 1.0               # massa especifica
C = 1.0               # calor especifico
s = k/(r*C)           # coeficiente de dufusao

#################################################################
########################           DOMINIO                         ####################

M = 12                           # mumeos de pontos sem a fronteira
L = 1.0                                    # comprimento do cubo 
h = L/(M+1)                              #espacamento dos pontos
t = 0.1                                   #espaçamento do tempo
x = np.linspace(0,L,M+2)
y = np.linspace(0,L,M +2)
z = np.linspace(0,L,M+2)

X,Y,Z = np.meshgrid(x,y,z)
U,V = np.meshgrid(x,z)

#################################################################
#####################   FUNCOES AUXILIARES   ############################

########                                 funcao da velocidade                                          #########
n,m,p = np.shape(X)

#u = -1.0e-10*(Y-1)**500 - (a-1.0e-10)* (Z-1)**500 +a 
u= np.zeros((n,m,p))
v =  np.zeros((n,m,p))
w = np.zeros((n,m,p))

auxnu=np.zeros((m,1))
auxnv=np.zeros((m,1))
auxnw=np.zeros((m,1))
for i in range(m):
    auxnu[i] = max([valor for linha in u[i] for valor in linha])
    auxnv[i] = max([valor for linha in v[i] for valor in linha])
    auxnw[i] = max([valor for linha in w[i] for valor in linha])
nu = max([valor for linha in auxnu for valor in linha])
nv= max([valor for linha in auxnv for valor in linha])
nw= max([valor for linha in auxnw for valor in linha])


###########                             Funcao da fonte                                    #################
 #(fazer o loop i.........5000)
b = t/h**3
q= 0

############                                Fronteira                                                  #############

f = 0

#################################################################
############                 COEFICIENTES DE CONTRIBUICAO           ###################

He = h*s+(0.5*h**2)*max(-nu,0)                           # leste
Hw = h*s+(0.5*h**2)*max(nu,0)                           # oeste
Hn = 0                           # norte     
Hs = 0                            # sul
Ht = h*s+(0.5*h**2)*max(-nw,0)                          # em cima  
Hb = h*s+(0.5*h**2)*max(nw,0)                           # em baixo
Sh = He+Hw+Hn+Hs+Ht+Hb                        
Hp = Sh -8*h*s                                              # ponto

Lpe = b*Hp +1                                     # ponto no metodo explicito
Lpi = b*Hp -1                                       # ponto no metos implicito
Le  = b*He                                            # contribuicao a leste
Lw = b*Hw                                          # contribuicao a oeste
Ln  = b*Hn                                            # contribuicao a norte
Ls  = b*Hs                                           # contribuicao a sul
Lt  = b*Ht                                           # contribuicao a cima
Lb  = b*Hb                                          # contribuicao a baixo


#################################################################

#                                                      MATRIZES DE CONTRIBUIÇÃO 

#############                                   Matriz A  (vizinhos)                # #################

Ae = mf.matrizA(m,Lpe,Ln,Ls,Le,Lw,Lb,Lt)   ## metodo explicito
Ai = mf.matrizA(m,Lpi,Ln,Ls,Le,Lw,Lb,Lt)    ## metodo implicito

#print Lpe
############                                    Matriz F  (fronteira)                             ##############

F = mf.matrizF(m,0,0,f,f,f,f) 

###########                                     Matriz G (fonte)                              ################

p =mf.centrobase(m)                         # posição no centro da base da malha com incluido a fronteira
G =mf.matrizG(m,p[0],p[1],p[2],q)    # matriz da fonte incluindo a fronteira



##########         IMPLEMETACAO  2D           ####################

B=10000 ##Numero de iteracoes

Pc = mf.centroPlaca(M)
#print(Pc)

Me = np.zeros((B,8))  ## coleta os dados para analise do erro

Taux = np.zeros((m,m))
T0e = np.zeros((m**3,1))
T0i = np.zeros((m**3,1))

d= np.pi/L

for i in range (m):
    for j in range(m):
        Taux[i,j] = np.sin(d*i*h)*np.sin(d*j*h) + f 

Taux1 = Taux.reshape((m**2,1))

for i in range (m**2):
    T0e[i] = Taux1[i]
    T0i[i] = Taux1[i]

######### METODO EXPLICITO
fig1 = plt.figure()

ax1 = fig1.add_subplot(131, projection = '3d')

teste =1 # numeros de ciclos a cada plot

Te= np.zeros((m**3,1))
for i in range(B):
    q= 0
    G =mf.matrizG(m,p[0],p[1],p[2],q)
    if (i==0):
        Te = T0e
    else:
        Te = np.dot(Ae,T0e) + F + G
    Qe = mf.plot2(Te,f,h) ### matriz que acumula as coordenadas com  respectivos valores dos pontos (explicito) 
    T0e =Te
    Z0 = Qe[:,3].reshape((m,m))
    Me[i,0]= i
    Me[i,1]= Qe[Pc,3]
    #print Me
#    if(i == teste):
#    #if (i%10==0):
#       ax1.plot_trisurf(Qe[:,0],Qe[:,2],Qe[:,3], cmap='jet', linewidth=0.1)
#       ax1.set_title('EXPLICITO', color = 'black')
#      # ax.scatter(Qe[:,0],Qe[:,2],Qe[:,3])
       
##### METODO IMPLICITO
ax2 =fig1.add_subplot(132, projection = '3d')
A = np.linalg.inv(Ai)
Ti= np.zeros((m**3,1))
for i in range(B):
    q= 0
    G =mf.matrizG(m,p[0],p[1],p[2],q)
    if(i==0):
        Ti = T0i
    else:
        T = -T0i - F - G
        Ti = np.dot(A,T)
        
    Qi = mf.plot2(Ti,f,h) ### matriz que acumula as coordenadas com  respectivos valores dos pontos (implicito) 
    T0i = Ti
    Z0 = Qi[:,3].reshape((m,m))
    Me[i,2] = Qi[Pc,3]
    #if (i%10==0):
#    if(i == teste):    
#       ax2.plot_trisurf(Qi[:,0],Qi[:,2],Qi[:,3], cmap='jet', linewidth=0.1)
#       #print '*****',np.shape(Qi[:,0]),np.shape(Qi[:,2])
#       ax2.set_title('IMPLICITO', color = 'black')
#       
#       ax.scatter(Qe[:,0],Qe[:,2],Qe[:,3])

########## METODO ANALITICO
#fig2 = plt.figure()
ax3 = fig1.add_subplot(133, projection = '3d')
d = np.pi/L
W= np.zeros((m**2,1))
s1= np.sqrt(s)
for i in range (B):
    W = np.sin(d*U)*np.sin(d*V)*np.exp((-2*(s1*d)**2)*t*i)
#    if (i%10==0):
#    if(i == teste):
#        u0 = U.reshape((m**2,1))
#        v0 = V.reshape((m**2,1))        
#        w0 = W.reshape((m**2,1))        
#        
#        #print np.shape(u0), np.shape(v0),np.shape(w0)
#        u0 = u0.flatten()
#        v0 = v0.flatten()
#        w0 = w0.flatten()
#        
#        ax3.plot_trisurf(u0, v0, w0,  cmap='jet', linewidth=0.1)
#        ax3.set_title('ANALITICO', color = 'black')
    W=W.reshape((m**2,1))
    Me[i,3] = W[Pc]


##### Plotando a matriz de erro #######

for i in range(B):
    Me[i,4]=( (Me[i,3]-Me[i,1])/Me[i,3]) *100  ## erro relativo explicito
    Me[i,5]=( (Me[i,3]-Me[i,2])/Me[i,3]) *100 ## erro relativo implicito
    Me[i,6]=np.abs(Me[i,3]-Me[i,1])    ## erro absoluto explicito
    Me[i,7]=np.abs(Me[i,3]-Me[i,2])    ## erro absoluto implicito



#
#print 'Erro médio relativo Explicito: ', np.mean(np.abs(Me[:,4])), '%'
#print 'Erro maximo relativo Explicito: ', np.max(np.abs(Me[:,4])), '%' 
#print
#print 'Erro médio relativo Implicito: ', np.mean(np.abs(Me[:,5])), '%'
#print 'Erro maximo relativo Implicito: ',  np.max(np.abs(Me[:,5])), '%' 
#
#fig3 = plt.figure()
#ax5 = fig3.add_subplot(1,1,1)
#ax5.plot(Me[:,0]*t,Me[:,1], color = 'r')  ## curva do metodo explicito
#ax5.plot(Me[:,0]*t,Me[:,2], color = 'b')  ## curva do metodo implicito
#ax5.plot(Me[:,0]*t,Me[:,3], color = 'g')  ## Curva solução analitica


fig4 = plt.figure()
ax4 = fig4.add_subplot(1,1,1)
plt.plot(Me[:,0]*t, Me[:,6], color ='r')   # propagacao erro explicito
plt.plot(Me[:,0]*t, Me[:,7], color ='b')  # propagacao erro implicito 
ax4.text(0.05, 0.95, "Erro 2D", transform=ax4.transAxes)

#print Me

#ax3.set_zlim(0,0.3)
 
plt.show
    
 

###############################  FIM     #############################
