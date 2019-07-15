# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 08:12:56 2019

@author: gideo
"""
import numpy as np

###################### FUNCAO TESTE DA CONDICAO DE CFL ####################################
def CFL(v, dt, dx):
    a = v * np.divide(dt,dx)
    if a <= 1:
        return a
    else:
        return "CFL nao atendido condicao deve ser menor ou igual a 1, ela esta:  ", a 
# v é norma da velocidade
# dt passo do tempo
# dx passo do espaco
############################################################################################

######################### FUNCAO QUE ACHA AS COORDENADAS DADA A POSICAO DO PONTO ###########
def indices(m,n):
    if (n >= m**3):
        return " estouro "
    else:
        i = n%m 
        j = int(n/m**2) 
        a = (n -(i-1))/m
        k = a%m    
        return int(i) ,int(j), int(k)
# n = a posicao do ponto na malha ordenada Po, P1, ....Pn
# m = a quantidade de pontos distribuido ao longo de um eixo do cubo, sem contar as bordas 
#############################################################################################


#################### FUNCAO QUE INFORMA: (p,F,F,F,F,F,F,P), sendo: p = Posição no sistema de coordenadas, 
##  F = fronteira que pertence e P = posisção do ponto. 
##  DADAO AS COORDEADAS DO PONTO ##################
def indicaP(i,j,k,m):
    mf = m+2
    if(0 <= i< mf and 0<=j< mf and 0<=k< mf):
        p = (mf*(k) +i) + (j)*mf**2         
        if(p% (mf**2) < mf): 
            a = 'BOTTON' 
        else:
            a = 'null'            
        if(mf**2 -mf <= p%(mf**2) <= mf**2):
            b = 'TOP'
        else:
            b = 'null'            
        if(p%mf == 0):
            c = 'WEST'
        else:
            c = 'null'            
        if(p%mf == mf-1):
            d = 'EAST'
        else:
            d = 'null'        
        if (p < mf**2):
            e = 'SOUTH'
        else:
            e = 'null'         
        if (mf**2*(mf-1) <= p < mf**3):
            f = 'NORTH'
        else:
            f = 'null'
        if(0 < i < mf-1 and 0<j< mf-1 and 0<k<= mf-1):
            P = (j-1)*(mf-2)**2 + (k-1)*(mf-2) +(i-1)   
        else:
            P = 'Fronteira'                
    return p,f,e,d,c,a,b,P
# mf = o numeros de pontos ao logo do eixo no dominio incluido a fronteira
#############################################################################################################

################################### MONTA A MATRIZ DE CONTRIBUICAO #################################

def matrizA(m,lp,ln,ls,le,lw,lb,lt):
    if m < 3:
        A = "Nao é possivel construir o cubo do dominio"   
    else:
        A = np.zeros((m**3,m**3))
        for i in range(m**3):
            A[i,i] = lp #Ponto
            # Direcao Y: Norte e sul
            if (i < (m**3-m**2)):
                A[i,i+ m**2 ] = ln #Norte
                A[i+ m**2,i] = ls #Sul
            
             # Direcao Z: Topo e boton
            if(i%(m**2) < m*(m-1)):                
                A[i,i+m] = lb  #Boton
                A[i+m,i] = lt  #Topo
    
            # Direcao X: Leste e Oeste
            if ((i+1)% m != 0):
                A[i,i+1] = le #Leste
                A[i+1,i] = lw #Oeste 
    return A
    
# m = numeros de ponto no eixo sem as fronteiras. l(p,n,s,e,w,b,t) sao os respectivos coeficientes  

############### MONTA A MATRIZ DE CONTRIBUIÇÃO DA FRONTEIRA ######################## 

def matrizFl(m,fn,fs,fe,fw,fb,ft):
    if m< 3:
        F="Não é possivel contruir o cubo do dominio"
    else:
        F= np.zeros((m**3,1))
        k=0
        j=0
        for i in range(m**3):
            if (i%m**2 < m):
                F[i] = F[i]  + fb
            if(m**2 -m <= i%(m**2) < m**2):
                F[i] = F[i] + ft 
            if(i%m ==0):
                F[i] = F[i] +fw[k]
                k=k+1
            if(i%m ==m-1):
                F[i] = F[i] + fe[j]
                j=j+1
            if(i< m**2):
                F[i] = F[i]+fs
            if(m**2*(m-1)<= i<m**3):
                F[i] = F[i]+fn
    return F 




def matrizF(m,fn,fs,fe,fw,fb,ft):
    if m< 3:
        F="Não é possivel contruir o cubo do dominio"
    else:
        F= np.zeros((m**3,1))
        for i in range(m**3):
            if (i%m**2 < m):
                F[i] = F[i]  + fb
            if(m**2 -m <= i%(m**2) < m**2):
                F[i] = F[i] + ft 
            if(i%m ==0):
                F[i] = F[i] +fe
            if(i%m ==m-1):
                F[i] = F[i] + fe
            if(i< m**2):
                F[i] = F[i]+fs
            if(m**2*(m-1)<= i<m**3):
                F[i] = F[i]+fn
    return F 

# m= muneros de pontos no eixo sem as fronteiras. f(n,s,e,w,b,t) são os fatores de contribuicao da fronteiras pelas respectivas faces

#######################  FUNCAO MONTA A MATRIZ QUE POSICIONA A  FONTE  #####################

def matrizG(m,x,y,z,q): 
    G = np.zeros((m**3,1))
    p = indicaP(x,y,z,m)[7]
    G[p] =  q  
    return G
    
#m = numeros de pontos em cada eixo, sem as fronteiras
# x,y,z sao as coordenadas do local da fonte
# q = t* a função da fonte

#######################   FUNCAO QUE LOCALIZA O PONTO CENTRAL ############
def centrobase(m):
    if m%2 != 0:
        i = (m-1)/2 +1
        k = (m-1)/2 +1 
        j = 1
        p=[i,j,k]
    
    if m%2 == 0:
        i = (m/2) +1
        k = (m/2)+1
        j = 1
        p=[i,j,k]      
    return p

def centroPlaca(m):
    M =m+2
    if m%2 != 0:
        for i in range(M**2):
            if(i%(M)== 0.5*(m+1) and i//(M) == 0.5*(m+1)):
                p=i
    if m%2 ==0:
        for i in range(m**2):
            if (i%m ==0.5*m and i/(m+2) ==0.5*m):
                p=i
    return p
######################## FUNCAO QUE  AUXILIAR PARA PLOTAR ###################
# acumula as coordenaa dao ponto + valor do ponto + fronteira. Recebe um vetor ordem m e sai uma matiz m x 4
def plot(X,f,dx):
    M = len(X)
    n = (M**(1.0/3.0))
    n1 = int(n)    
    if ((n - n1) > 0.5):
        n = n1 + 1
    else:
        n = n1 
    n =int(n+2)
    A = np.zeros((n**3,4))
    j=0  
    for i in range(n**3):
        if(indices(n,i)[0] ==0.0 or indices(n,i)[1] ==0.0 or indices(n,i)[2] ==0.0 or indices(n,i)[0] ==n-1 or indices(n,i)[1] ==n-1 or indices(n,i)[2] ==n-1):
            A[:,0][i] = indices(n,i)[0]*dx
            A[:,1][i] = indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx
            A[:,3][i] = f      
        else:
            A[:,0][i] =  indices(n,i)[0]*dx
            A[:,1][i] =  indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx
            A[:,3][i] = X[j]
            j =j+1   
    return A


def plott(X,f,fe,Le,fw,Lw,dx):
    M = len(X)
    n = (M**(1.0/3.0))
    n1 = int(n)    
    if ((n - n1) > 0.5):
        n = n1 + 1
    else:
        n = n1 
    n =int(n+2)
    A = np.zeros((n**3,4))
    j=0
    k=0
    l=0
    m=0
    N = n**3
    for i in range(N):      
            A[:,0][i] = indices(n,i)[0]*dx
            A[:,1][i] = indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx                     
            if( indices(n,i)[1] ==0.0 or indices(n,i)[2] ==0.0 or indices(n,i)[1] ==n-1 or indices(n,i)[2] ==n-1):
                A[:,3][i] = f
                m=m+1
            if(indices(n,i)[0] ==0.0 and indices(n,i)[1] !=0.0 and indices(n,i)[2] !=0.0 and indices(n,i)[0] !=n-1 and indices(n,i)[1] !=n-1 and indices(n,i)[2] !=n-1):
                A[:,3][i] = fe[k]*(Le**(-1))
                k = k+1
            if(indices(n,i)[0] ==n-1 and indices(n,i)[1] !=0.0 and indices(n,i)[2] !=0.0 and indices(n,i)[1] !=n-1 and indices(n,i)[2] !=n-1):
                A[:,3][i] = fw[l]*(Lw**(-1))
                l=l+1
            if(indices(n,i)[0] !=0.0 and indices(n,i)[1] !=0.0 and indices(n,i)[2] !=0.0 and indices(n,i)[0] !=n-1 and indices(n,i)[1] !=n-1 and indices(n,i)[2] !=n-1):
                A[:,3][i] = X[j]
                j =j+1   
    return A
 



def plot1(X,f,dx):
    M = len(X)
    n = (M**(1.0/3.0))
    n1 = int(n)    
    if ((n - n1) > 0.5):
        n = n1 + 1
    else:
        n = n1 
    A = np.zeros((n,4))
    for i in range(n):
        if(indices(n,i)[0] ==0.0 or indices(n,i)[0] ==n-1.0 ):
            A[:,0][i] = indices(n,i)[0]*dx
            A[:,1][i] = indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx
            A[:,3][i] = f
        else:
            A[:,0][i] =  indices(n,i)[0]*dx
            A[:,1][i] =  indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx
            A[:,3][i] = X[i]         
    return A

def plot2(X,f,dx):
    M = len(X)
    n = (M**(1.0/3.0))
    n1 = int(n)    
    if ((n - n1) > 0.5):
        n = n1 + 1
    else:
        n = n1 
    A = np.zeros((n**2,4))
    for i in range(n**2):
        if(indices(n,i)[0] ==0.0 or indices(n,i)[2] ==0.0 or indices(n,i)[0] ==n-1.0 or indices(n,i)[2] ==n-1.0 ):
            A[:,0][i] = indices(n,i)[0]*dx
            A[:,1][i] = indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx
            A[:,3][i] = f
        else:
            A[:,0][i] =  indices(n,i)[0]*dx
            A[:,1][i] =  indices(n,i)[1]*dx
            A[:,2][i] = indices(n,i)[2]*dx
            A[:,3][i] = X[i]
              
    return A

##############################################################################

################ #########                    PLOTS                      #######################

# funcao pega corta dos planos para imprimir 
#recebe  uma matriz m3 x 4 e sai uma matriz de ordem m2  para os eixo dois eixos e o outro é o valor do ponto

def cortXY (Q,m,z):    ## apresenta p plano xy
    n =m**2
    F= np. zeros((n,3))
    i = 0
    for j in range(m**3):        
       if (Q[j][2]==z):         
           F[i][0] = Q[j][0] 
           F[i][1] = Q[j][1]
           F[i][2] = Q[j][3]
           i = i+1
    X = F[:,0].reshape((m,m))
    Y = F[:,1].reshape((m,m))
    Z = F[:,2].reshape((m,m))
    return X,Y, Z

def cortZY (Q,m,x):    ## apresenta p plano xy
    n =m**2
    F= np. zeros((n,3))
    i = 0
    for j in range(m**3):        
       if (Q[j][0]==x):         
           F[i][0] = Q[j][2] 
           F[i][1] = Q[j][1]
           F[i][2] = Q[j][3]
           i = i+1
    X = F[:,0].reshape((m,m))
    Y = F[:,1].reshape((m,m))
    Z = F[:,2].reshape((m,m))
    return X,Y, Z


def cortXZ (Q,m,y):   ## apresenta o plnano XZ
    n =m**2
    F= np. zeros((n,3))
    i = 0
    for j in range(m**3):        
       if (Q[j][1]==y):         
           F[i][0] = Q[j][0] 
           F[i][1] = Q[j][2]
           F[i][2] = Q[j][3]
           i = i+1
    X = F[:,0].reshape((m,m))
    Y = F[:,1].reshape((m,m))
    Z = F[:,2].reshape((m,m))
    return X,Y,Z

###########################   FIM   ###########################################


# resolver o proble da matriz X que entra com n**3 elemento e deve ser colocada no centro do cubo com m**3 elementos  