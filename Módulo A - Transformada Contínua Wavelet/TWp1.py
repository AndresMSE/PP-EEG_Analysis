#!/usr/bin/env python
# coding: utf-8

# In[ ]:


'''En este archivo de Python 3 se almacenarán todas las funciones implementadas y discutidas en el notebook
Transformada Continua Wavelet - Fundamentos para su exportación en otras aplicaciones y notebooks'''

import scipy as sc  #Paq. p-análisis numérico 
import numpy as np  #Paq. p-cálculo numérico

def Suma (x): #La función recibira una variable de tipo array
    suma = 0 #Definimos el valor inicial de la suma en 0
    for i in range(len(x)): #iniciamos un ciclo for durante el cual se realizará la operación
        suma += abs(x[i])**2 #iterativamente sumamos el valor anterior de la variable más el elemento i del array
    return suma
def p_interior(f,g): #Definimos nuestra función para dos arrays f y g que deberán tener la misma longitud de puntos
    resultado = 0 #Nombramos a la variable donde almacenaremos el resultado con su valor inicial
    for i in range(len(f)): #Establecemos un bucle para las operaciones iterativas
        resultado += f[i]*g[i] #Obtenemos el producto de los elementos f_i y g_i y lo sumamos iterativamente
    return resultado #La función arroja el resultado del producto interior
def conv_tiempo(f,g): #Definimos nuestra función para dos arrays f y g  
    n,m = len(f),len(g) #Obtenemos las longitudes de las series
    k = int(m/2) #Obtenemos el número de elementos a añadir en f para el encamado
    '''Encamado de la serie f'''
    f_pad = np.pad(f,(k,k),'constant') #Encamamos la serie usando el módulo de numpy pad()
    #El primer parámetro es el array a encamar, el segundo es la longitud de puntos a añadir al inicio 
    #y al final (k,k), el último parámetro indica que los valores serán constantes igual a 0
    '''Proceso de convolución'''
    p=n+m-1 #Longitud del resultado de la convolución
    conv = np.zeros(p) #Definimos el array donde se almacenarán los resultados de los productos punto
    for ti in range(len(f_pad)-m): #Establecemos el ciclo iterativo para los productos interiores
        f_temp = f_pad[ti:ti+m] #Extraemos el segmento de la serie f correspondiente
        g_temp = np.flip(g) #Invertimos la serie g 
        conv[ti+k] = p_interior(f_temp,g_temp) #Realizamos el producto interior y lo almacenamos
    '''Recortado del resultado'''
    conv_cut = conv[k:-k+1]
    return conv_cut
def cuadrado(m): #Función para crear una serie de pulso cuadrado con longitud m 
    pulso = np.ones(m) #Definimos que el pulso tendrá m número de unos
    pulso = np.pad(pulso,(10,10),'constant') #Establecemos que el pulso tendrá 10 valores 0 a los costados
    return pulso
def conv_v2(f,g): #Definimos la función para dos arrays f y g
    'Realizamos un padding a las series para que tengan el mismo número de puntos que el del resultado de la conv. n+m-1'
    n,m =len(f),len(g)
    f_pad = np.pad(f,(0,m-1),'constant')
    g_pad = np.pad(g,(0,n-1),'constant')
    #En ambos casos añadimos puntos para que al final la longitud de ambas series sea igual a p=n+m-1
    '''Cambio al dominio de la frecuencia por FFT1'''
    f_fft = np.fft.fft(f_pad)
    g_fft = np.fft.fft(g_pad)
    g_fft = g_fft/g_fft.max() #Normalizamos el espectro de g para conservar las unidades de f en la convolución
    '''Teorema de la convolución'''
    conv = np.multiply(f_fft,g_fft)
    conv_res = np.fft.ifft(conv)
    '''Recorte del resultado para que tenga longitud n'''
    k = int(len(g)/2)
    conv_res = abs(conv_res[k:-k+1])
    return conv_res
def morlet(fw,n_c,twav,freqm=500): 
    #Componente oscilatoria 'Seno complejo'
    com_sin = np.exp(1j*2*np.pi*fw*twav)
    #Ventana de Gauss
    sigma = n_c/(2*np.pi*fw)
    gauss_c = np.exp(-twav**2/(2*sigma**2))
    #Normalization
    A_sigma = 1/(np.sqrt(sigma*np.sqrt(np.pi)))
    #Morlet wavelet
    wavelet = A_sigma*com_sin*gauss_c
    return wavelet
def gauss(twav,n_c,fw):
    #Ventana de Gauss
    sigma = n_c/(2*np.pi*fw)
    gauss_c = np.exp(-twav**2/(2*sigma**2))
    #Normalization
    A_sigma = 1/(np.sqrt(sigma*np.sqrt(np.pi)))
    return gauss_c*A_sigma

