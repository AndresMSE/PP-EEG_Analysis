#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''Funciones generadas en el notebook "Singular Spectrum Analysis - Fundamentos'''

import scipy as sc  #Paq. p-análisis numérico 
import numpy as np  #Paq. p-cálculo numérico
from scipy import linalg
import math # math module


def ScreeD(signal,L): # La función toma una señal y un parámetro L válido
    '''PREPARACIÓN'''
    N=len(signal)
    K=N-L+1
    '''ENCAMADO'''
    X=np.array([[signal[i+j] for j in range(0,L)] for i in range(0,K)]) # matriz de trayectoria
    Xh = X.transpose() #matriz de trayectoria transpuesta
    S = np.matmul(Xh,X) #Calculo de la matriz cuadrada S
    s=linalg.svd(S,compute_uv=False) # descomposición en valores singulares
    return s
def SSA(signal, L):
    '''PREPARACIÓN - Obtenemos el número de elementos de la señal y de vectores columna'''
    N = len(signal) #no. de puntos de la señal
    K = N-L+1  #no. de columnas de la matriz
    '''ENCAMADO - Obtención de la matriz de trayectoria'''
    X = np.array([[signal[i+j] for j in range(0,L)] for i in range(0,K)]) #Matriz de Trayectoria
    '''SVD - Descomposición'''
    U,s,V = linalg.svd(X) #Aplicación del módulo linalg para la descomposicón SVD
    G = len(s) #no. de eigentripletas
    l = s**2 #Valores propios ó varianzas parciales
    X_elem = [] #Lista auxiliar donde se almacenarán las eigentripletas
    gkList = np.zeros(shape=(G,N)) #array de 2D donde se registrarán las G componentes de N puntos a obtener 
    for k in range(0,G):
        Uk = U[:,k] #Vector U de la eigentripleta k
        Vk = V[k,:] #Vector V de la eigentripleta k 
        X_k = s[k]*np.outer(Uk,Vk) #Calculo de la eigentripleta k
        X_elem.append(X_k) #Añadimos la eigentripleta  k a la lista
        '''PROMEDIACIÓN DIAGONAL - de eigentripletas a vectores'''
        gk = [] #lista auxiliar donde se almacenarán las series de tiempo of values of time series
        for i in range(min(K-1,L-1),-max(K-1,L-1)-1,-1): # loop sobre las diagonales
            gki=np.mean(np.diag(np.fliplr(X_k),i)) #Promediación diagonal - valores sucesivos de la serie de tiempo
            gk.append(gki) #Añadimos la serie a la lista
        gkList[k]=gk #Registramos la serie en el array de series de tiempo
    '''MATRIZ DE CORRELACIÓN'''
    w=[] # Lista auxiliar donde se almacenarán los valores de los pesos
    #Cotas para el índice 
    LL=min(L,K)
    KK=max(L,K)
    for ll in range(1,LL+1): # primero tercio 
        w.append(ll)
    for ll in range(LL+1,KK+1): # segundo tercio
        w.append(LL)
    for ll in range(KK+1,N+1): # tercer tercio
        w.append(N-ll)
    kMin=kkMin=0 # Establecemos que se muestre la matriz de correlación solo para las primeras 20 componentes 
    if L >= 20:
        kMax=kkMax=20
    else:
        kMax=kkMax=L
    #Cálculo de la matris de correlación
    wMatriz=[[sum(w*gkList[k]*gkList[kk])/(math.sqrt(sum(w*gkList[k]*gkList[k]))*math.sqrt(sum(w*gkList[kk]*gkList[kk]))) for k in range(kMin,kMax)] for kk in range(kkMin,kkMax)]
    wMatriz=np.array(wMatriz)
    return (l, gkList, wMatriz) #Añadimos la matriz de correlación a los resultados

