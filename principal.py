# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 08:13:53 2019

@author: gideo
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import minhasFuncoes as mf

#### #################        DADOS DO PRPBLEMA                              ##############

k = 1.0e-3              # coeficiente de transferencia de calpor
r = 1.0               # massa especifica
C = 1.0             # calor especifico
s = k/(r*C)           # coeficiente de dufusao

#################################################################
########################           DOMINIO                         ####################

M = 13                          # mumeos de pontos sem a fronteira
L = 1.0                                    # comprimento do cubo 
h = L/(M+1)                              #espacamento dos pontos
t = 0.1                                  #espaçamento do tempo
x = np.linspace(0,L,M+2)
y = np.linspace(0,L,M +2)
z = np.linspace(0,L,M+2)

X,Y,Z = np.meshgrid(x,z,y)

#################################################################
#####################   FUNCOES AUXILIARES   ############################

########                                 funcao da velocidade                                          #########
n,m,p = np.shape(X)
a=0.1  ## controla a norma da velocidade


#u = -1.0e-10*(Y-1)**500 - (a-1.0e-10)* (Z-1)**500 +a 
u = a*(Y**2) * (Y-L)*(Z**2)*(Z-L)
#u= np.zeros((n,m,p))
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

print 'Condição CFL: ', c

###########                             Funcao da fonte                                    #################
 #(fazer o loop i.........5000)
b = t/h**3
#q= b*2.0e+4* np.exp ( -1.0e-6*(i-2400.0)**2) +25.0
#q = 5000

#################################################################
############                 COEFICIENTES DE CONTRIBUICAO           ###################

He = h*s+(0.5*h**2)*max(-nu,0)                           # leste
Hw = h*s+(0.5*h**2)*max(nu,0)                           # oeste
Hn = h*s+(0.5*h**2)*max(-nv,0)                           # norte     
Hs = h*s+(0.5*h**2)*max(nv,0)                            # sul
Ht = h*s+(0.5*h**2)*max(-nw,0)                          # em cima  
Hb = h*s+(0.5*h**2)*max(nw,0)                           # em baixo
Sh = He+Hw+Hn+Hs+Ht+Hb                        
Hp = Sh -12*h*s                                              # ponto

Lpe = b*Hp +1                                     # ponto no metodo explicito
Lpi = b*Hp -1                                       # ponto no metos implicito
Le = b*He                                            # contribuicao a leste
Lw = b*Hw                                          # contribuicao a oeste
Ln = b*Hn                                            # contribuicao a norte
Ls = b*Hs                                           # contribuicao a sul
Lt = b*Ht                                           # contribuicao a cima
Lb = b*Hb                                          # contribuicao a baixo

#print c, nu, Lpe
#################################################################


############                             fator de contribuição da  Fronteira                #############


f = 25

fn = f*Ln
fs = f*Ls
fe = f*Le
fw = f*Lw
fb = f*Lb
ft = f*Lt

fel = np.zeros(M**2)
fwl = np.zeros(M**2)

F = mf.matrizF(M,fn,fs,fe,fw,fb,ft) # Matris de contribuicao da fronteira

#print fn,fs,fe,fw,fb,ft
#                                                      MATRIZES DE CONTRIBUIÇÃO 

#############                                   Matriz A  (vizinhos)                # #################

Ae = mf.matrizA(M,Lpe,Ln,Ls,Le,Lw,Lb,Lt)   ## metodo explicito
Ai = mf.matrizA(M,Lpi,Ln,Ls,Le,Lw,Lb,Lt)    ## metodo implicito



###########                                     Matriz G (fonte)                              ################

p =mf.centrobase(M)                         # posição no centro da base da malha com incluido a fronteira
G =np.zeros((M**3,1))   # matriz da fonte incluindo a fronteira
#print Ae
#print
#print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHh'
#print Ai 
###########################################################################
########## ######     PLOTAGEM DO CAMPO DE VELOCIDADE   ############################

#fig = plt.figure()
#ax = fig.add_subplot(111, projection = '3d')
#ax.quiver(X, Y, Z, u, v ,w ,  length=0.1)
#plt.show

#####$$$$$$$$$$$$$$     Dados  do problema

B=2001  ##Numero de iteracoes
r = 5   # Posicao do corte

############################ FONTE DE CALOR ###############

q = np.zeros(B)
for i in range(B):
    q[i] = 200 #b*2.0e+4* np.exp ( -1.0e-6*(i*t-2400.0)**2)

#print q

##########################
##########    IMPLEMETACAO FRONTEIRA TEMPERATURA CONSTANTE (METODO EXPLICITO) ####################

T0e = 25*np.ones((M**3,1))
Te= np.zeros((M**3,1))


#######
#for i in range(B):
#    if (i==0):        
#        Te = T0e +  mf.matrizG(M,p[0],p[1],p[2],q[i])
#    else:
#        G =mf.matrizG(M,p[0],p[1],p[2],q[i])
#        Te = np.dot(Ae,T0e) + F + G
#    Qe = mf.plot(Te,f,h) ### matriz que acumula as coordenadas com  respectivos valores dos pontos (explicito) 
#    T0e=Te 
#
#    if (i%500==0): 
#        fig1 = plt.figure()  
#        escalas =mpl.colors.Normalize(vmin =25, vmax =1000)

#print Qe

#####################  IMPLEMENTACAO FRONTEIRA LIVRE (METODO EXPLICITO)########################################

for i in range (B):
    if (i==0):
        Te = T0e +  mf.matrizG(M,p[0],p[1],p[2],q[i])
    else:
        G =mf.matrizG(M,p[0],p[1],p[2],q[i])
        contador=0
        j=0   
        k=0
        while contador < M**3:
            if contador % M ==0:
                fel[j] = Te[contador]*Le
                j = j+1
            if(contador+1)%M==0:
                fwl[k] = Te[contador]*Lw
                k=k+1
            contador = contador +1
    Fl = mf.matrizFl(M,fn,fs,fel,fwl,fb,ft)   # matriz de contricuicao da fronteira
    T0e = np.dot(Ae,Te) + Fl + G
    Te=T0e
    Qe = mf.plott(Te,f,fel,Le,fwl,Lw,h)  # matriz auxiliar de plotagem

    if (i%500==0): 
        fig1 = plt.figure()  
        escalas =mpl.colors.Normalize(vmin =25, vmax =1200) 
#print fe  

############################################################################################3
############ corte vertical XY
                
#        XYe = mf.cortXY(Qe,m,r*h)  
#        ax1 = fig1.add_subplot(1,1,1)
#        im1 = ax1.contourf(XYe[0], XYe[1], XYe[2], 20, alpha=.75, cmap='jet', norm=escalas)
#        C1 = ax1.contour(XYe[0], XYe[1], XYe[2], 20, colors='black', linewidth=.5)
#        ax1.clabel(C1,inline=True,fmt='%1.1f',fontsize=10)
#        fig1.colorbar(im1)
#        ax1.set_title(' METODO EXPLICITO', color = 'black')       
#        plt.show()
#print ' eixo XY em z = ', r*h        
        
############ corte vertical YZ
##        
#        XYe = mf.cortZY(Qe,m,r*h)  
#        ax1 = fig1.add_subplot(1,1,1)
#        im1 = ax1.contourf(XYe[0], XYe[1], XYe[2], 20, alpha=.75, cmap='jet', norm=escalas)
#        C1 = ax1.contour(XYe[0], XYe[1], XYe[2], 20, colors='black', linewidth=.5)
#        ax1.clabel(C1,inline=True,fmt='%1.1f',fontsize=10)
#        fig1.colorbar(im1)
#        ax1.set_title(' METODO EXPLICITO', color = 'black')       
#        plt.show()
#print ' eixo ZY em x = ', r*h              

     
############## corte horizontal
#        
        XZe = mf.cortXZ(Qe,m,r*h)        
        ax1 = fig1.add_subplot(1,1,1)
        im1 = ax1.contourf(XZe[0], XZe[1], XZe[2], 20, alpha=0.75, cmap='jet', norm=escalas)
        C1 = ax1.contour(XZe[0], XZe[1], XZe[2], 20, colors='black', linewidth=.5)
        ax1.clabel(C1,inline=True,fmt='%1.1f',fontsize=10)
        fig1.colorbar(im1)
        ax1.set_title(' METODO EXPLICITO', color = 'black')
        plt.show()
print ' eixo XZ em y = ', r*h              

############# IMPLEMENTACAO FRONTEIRA TEMPERATURA CONSTANTE (METODO IMPLICITO) ##################
    
#T0i = 25*np.ones((M**3,1))
#Ti= np.zeros((M**3,1))
#A = np.linalg.inv(Ai)
#
#
#for i in range(B):
#    if (i==0):
#        Ti = T0i + mf.matrizG(M,p[0],p[1],p[2],q[i])
#    else:
#        G =mf.matrizG(M,p[0],p[1],p[2],q[i])
#        T = -T0i - F - G
#        Ti = np.dot(A,T)
#    Qi= mf.plot(Ti,f,h)  ### matri que acumula as coordenadas com  respectivos valores dos pontos (Implicito)
#    T0i = Ti   
#     
#    if (i%500 ==0):
#        fig1 = plt.figure()  
#        escalas =mpl.colors.Normalize(vmin =25, vmax =1000)    

########### corte verftical XY
#            
#        XYi = mf.cortXY(Qi,m,r*h)        
#        ax2 = fig1.add_subplot(1,1,1)
#        im2 = ax2.contourf(XYi[0], XYi[1], XYi[2], 20, alpha=.75, cmap='jet', norm=escalas)
#        C2 = ax2.contour(XYi[0], XYi[1], XYi[2], 20, colors='black', linewidth=.5)
#        ax2.clabel(C2,inline=True,fmt='%1.1f',fontsize=10)
#        fig1.colorbar(im2)
#        ax2.set_title(' METODO IMPLICITO', color = 'black')
#        plt.show  
#print ' eixo XY em z = ', r*h      
    
############  Corte vertical ZY    
#    
#          
#        XYi = mf.cortZY(Qi,m,r*h)        
#        ax2 = fig1.add_subplot(1,1,1)
#        im2 = ax2.contourf(XYi[0], XYi[1], XYi[2], 20, alpha=.75, cmap='jet', norm=escalas)
#        C2 = ax2.contour(XYi[0], XYi[1], XYi[2], 20, colors='black', linewidth=.5)
#        ax2.clabel(C2,inline=True,fmt='%1.1f',fontsize=10)
#        fig1.colorbar(im2)
#        ax2.set_title(' METODO IMPLICITO', color = 'black')
#        plt.show
#print ' eixo ZY em x = ', r*h         

#            
############# corte horizontal
#           
#        XZi = mf.cortXZ(Qi,m,r*h)
#        ax2 = fig1.add_subplot(1,1,1)
#        im2 = ax2.contourf(XZi[0], XZi[1], XZi[2], 20, alpha=.75, cmap='jet', norm=escalas)
#        C2 = ax2.contour(XZi[0], XZi[1], XZi[2], 20, colors='black', linewidth=.5)
#        ax2.clabel(C2,inline=True,fmt='%1.1f',fontsize=10)
#        fig1.colorbar(im2)
#        ax2.set_title(' METODO IMPLICITO', color = 'black')
#        plt.show  
#print ' eixo XZ em y= ', r*h                

###############################  FIM     #############################
