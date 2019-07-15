# -*- coding: utf-8 -*-
"""
Created on Thu May 30 21:41:44 2019

@author: gideo
"""

##############  TESET 1 DIMENSÃO ################
# barra finita, sem fonte 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import minhasFuncoes as mf

#### #################        DADOS DO PRPBLEMA                              ##############

k =1.0e-3       # coeficiente de transferencia de calpor
r = 1.0               # massa especifica
C = 1.0               # calor especifico
s = k/(r*C)           # coeficiente de dufusao

#################################################################
########################           DOMINIO                         ####################

M = 12                         # mumeos de pontos sem a fronteira
L = 1.0                                    # comprimento do cubo 
h = L/(M+1)                              #espacamento dos pontos
t = 0.1                                  #espaçamento do tempo
x = np.linspace(0,L,M+2)
y = np.linspace(0,L,M +2)
z = np.linspace(0,L,M+2)

X,Y,Z = np.meshgrid(x,z,y)
#
##################################################################
######################   FUNCOES AUXILIARES   ############################
#
#########                                 funcao da velocidade                                          #########
n,m,p = np.shape(X)
a=0.1  ## controla a norma da velocidade

#u = -1.0e-10*(Y-1)**500 - (a-1.0e-10)* (Z-1)**500 +a 
#u = a*(Y**2) * (Y-L)*(Z**2)*(Z-L)

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

c = mf.CFL(nu,t,h)   #  condicao CFL eixo x
#print (c)
#print (u,nu,nv,nw)

############                             Funcao da fonte                                    #################
 #(fazer o loop i.........5000)
b = t/h**3
#q= b*2.0e+4* np.exp ( -1.0e-6*(i-2400.0)**2) +25.0
q = 0

#############                                Fronteira                                                  #############

f = 0

##################################################################
#############                 COEFICIENTES DE CONTRIBUICAO           ###################

He = h*s+(0.5*h**2)*max(-nu,0)                           # leste
Hw = h*s+(0.5*h**2)*max(nu,0)                           # oeste
Hn = 0                                                             # norte     
Hs = 0                             # sul
Ht = 0                           # em cima  
Hb = 0                          # em baixo
Sh = He+Hw+Hn+Hs+Ht+Hb                        
Hp = Sh -4*h*s                                              # ponto

Lpe = b*Hp +1                                     # ponto no metodo explicito
Lpi = b*Hp -1                                       # ponto no metos implicito
Le = b*He                                            # contribuicao a leste
Lw = b*Hw                                          # contribuicao a oeste
Ln =b*Hn                                            # contribuicao a norte
Ls =b*Hs                                           # contribuicao a sul
Lt = b*Ht                                           # contribuicao a cima
Lb = b*Hb                                          # contribuicao a baixo


####,##############################################################
#
##                                                      MATRIZES DE CONTRIBUIÇÃO 
#
##############                                   Matriz A  (vizinhos)                # #################

Ae = mf.matrizA(m,Lpe,Ln,Ls,Le,Lw,Lb,Lt)   ## metodo explicito
Ai = mf.matrizA(m,Lpi,Ln,Ls,Le,Lw,Lb,Lt)    ## metodo implicito


#############                                    Matriz F  (fronteira)                             ##############
#
F = mf.matrizF(m,0,0,f,f,0,0) 
#
############                                     Matriz G (fonte)                              ################
#
p =mf.centrobase(m)                         # posição no centro da base da malha com incluido a fronteira
G =mf.matrizG(m,p[0],p[1],p[2],q)    # matriz da fonte incluindo a fronteira


 
##########         IMPLEMETACAO  MODELO 1D         ####################
#
B=10000  ##Numero de iteracoes

xf=(M+1)//2  ## ponto para analise do erro
Me = np.zeros((B,8))  ## coleta os dados para analise do erro

T0e = np.zeros((m**3,1))
T0i = np.zeros((m**3,1))
d = np.pi/L

for i in range (m**3):
    if(i<m):
        T0e[i]= np.sin(d*i*h) +f
        T0i[i]= np.sin(d*i*h) +f


Te= np.zeros((m**3,1))
#
#fig1 = plt.figure() 
#ax1 = fig1.add_subplot(1,3,1) 


#########         IMPLEMETACAO  E PLOT DO METODO EXPLICITO           ####################

U= np.zeros((m,1))
erroe=np.zeros((B,1))

for i in range(B):
   q = 0
   G =mf.matrizG(m,p[0],p[1],p[2],q)
   
   if (i ==0):
        Te = T0e
   else:
        Te = np.dot(Ae,T0e) + F + G     
      
   Qe = mf.plot1(Te,f,h) ### matriz que acumula as coordenadas com  respectivos valores dos pontos (explicito)     
   T0e=Te 
   Me[i,0]=i
   Me[i,1]=Qe[xf,3]  ## valor expicito
   #print Me
#   if (i%20==0): 
#        #print Qe[:,3]
#        ax1.plot(Qe[:,0],Qe[:,3],c='r')  
#        ax1.set_title(' METODO EXPLICITO', color = 'black')  
#        #print Qe[:,3]         
#        plt.show()
     
#### #########                   METODO IMPLICITO            ########################

Ti= np.zeros((m**3,1))

A = np.linalg.inv(Ai)

#fig3=plt.figure()
#ax3 = fig1.add_subplot(1,3,2) 
for i in range(B):
        q = 0        
        G =mf.matrizG(m,p[0],p[1],p[2],q)
        if(i==0):
            Ti=T0i
        else:
            T = -T0i - F - G
            Ti = np.dot(A,T)
           
        Qi= mf.plot1(Ti,f,h)  ### matri que acumula as coordenadas com  respectivos valores dos pontos (Implicito)
        T0i = Ti 
        Me[i,2]=Qi[xf,3] ## valor implicito
        
#        if (i%20==0): 
#            #print Qe[:,3]
#            ax3.plot(Qi[:,0],Qi[:,3],c='b')  
#            ax3.set_title(' METODO IMPLICITO', color = 'black')
#            
#            plt.show()
#print"Qi = ", Qi      
#print "Qe = ", Qe
#
#erroQ = np.abs(Qi[:,3]-Qe[:,3])
#
#print 'erro medio = ', np.mean(erroQ)
#print 'erro max = ', np.max(erroQ)
#print 'erro min = ', np.min(erroQ)



### SOLUCAO ANALITICA ################
s1 =np.sqrt(s)
#fig2 = plt.figure()
#ax2 = fig1.add_subplot(1,3,3)

Ua = np.zeros((B,len(x)))
for i in range(B):
    for j in range(len(x)):
        Ua[i,j] = np.sin(d*x[j]) * np.exp(-i*t*(d*s1)**2)
        Me[i,3]=Ua[i,xf]  ## valor analitico


#for i in range(B):
#   U = np.sin(d*x) * np.exp(-i*t*(d*s1)**2) 
#   if(i%20==0):         
#         ax2.set_title(' ANALITICA', color = 'black')
#         ax2.plot(x,U,c='g')


#
#fig2 = plt.figure()
#ax4 = fig2.add_subplot(1,1,1)
#
#ax4.plot(x,Qi[:,3], c='b')
#ax4.plot(x,Qe[:,3], c='r')
#ax4.plot(x,U, c='g')

#### plotando o erro de um ponto Xf

for i in range(B):
    Me[i,4]=( (Me[i,3]-Me[i,1])/Me[i,3]) *100  ## erro relativo explicito
    Me[i,5]=( (Me[i,3]-Me[i,2])/Me[i,3]) *100 ## erro relativo implicito
    Me[i,6]=np.abs(Me[i,3]-Me[i,1])    ## erro absoluto explicito
    Me[i,7]=np.abs(Me[i,3]-Me[i,2])    ## erro absoluto implicito



erroMedE = np.mean( np.abs(Me[:,3]-Me[:,1]))
erroMedi = np.mean( np.abs(Me[:,3]-Me[:,2])) 

#print 'Erro médio relativo Explicito: ', np.mean(np.abs(Me[:,4])), '%'
#print 'Erro maximo relativo Explicito: ', np.max(np.abs(Me[:,4])), '%' 
#print
#print 'Erro médio relativo Implicito: ', np.mean(np.abs(Me[:,5])), '%'
#print 'Erro maximo relativo Implicito: ',  np.max(np.abs(Me[:,5])), '%'

#fig3 = plt.figure()
#ax3 = fig3.add_subplot(1,1,1)
#ax3.plot(Me[:,0]*t,Me[:,1], color = 'r')  ## curva do netodo explicito
#ax3.plot(Me[:,0]*t,Me[:,2], color = 'b') ## curva do metodo implicito
#ax3.plot(Me[:,0]*t,Me[:,3], color = 'g') ## curva da soulução anlitica

fig4 = plt.figure()
ax4 = fig4.add_subplot(1,1,1)
#plt.loglog(Me[:,0]*t, Me[:,6], color ='r')   # propagacao erro explicito
#plt.loglog(Me[:,0]*t, Me[:,7], color ='b')  # propagacao erro implicito 
#plt.semilogy(Me[:,0]*t, Me[:,6], color ='r')   # propagacao erro explicito
#plt.semilogy(Me[:,0]*t, Me[:,7], color ='b')  # propagacao erro implicito 
#plt.plot(np.log(Me[:,0]*t), np.log(Me[:,6]), color ='r')   # propagacao erro explicito
#plt.plot(np.log(Me[:,0]*t), np.log(Me[:,7]), color ='b')  # propagacao erro implicito 


plt.plot(Me[:,0]*t, Me[:,6], color ='r')   # propagacao erro explicito
plt.plot(Me[:,0]*t, Me[:,7], color ='b')  # propagacao erro implicito
ax4.text(0.05, 0.95, "Erro 1D", transform=ax4.transAxes)


#print Me

plt.show












#############

