'''Las variables que utiliza la función es una señal de valores y (dependiente), asociada a un vector 
de puntos x (independinte)'''
def FWHM(y,x): #Definimos la función que nos calculará la amplitud a media anchura FWHM de una señal x
    ynorm = y / y.max() #normalizamos la señal
    x1,x2 =0,0 #definimos el valor inicial de las variables internas de la función
    x1_i,x2_i = 0,0
    resolution = 1e-2 #nivel de resolución de punto flotante para la detección del valor de media anchura
    'Ciclos for donde se buscarán los valores a media anchura para antes y después del máximo'
    for i in range(len(x[0:len(x)//2])):
        if 0.5 <= y[i] < 0.5 + resolution:
            x1 = x[i] #Almacenamos el punto en x donde se ubica el valor y tal que es igual a 0.5+/- 1e-3
            x1_i = i #Almacenamos el índice para dicho valor x 
            break
        else: 
    for i in range(len(x[len(x)//2:])):
        if 0.5 <= y[i+len(x)//2] < 0.5 + resolution:
            x2 = x[i+len(x)//2]
            x2_i = i+len(x)//2
    delta = x2-x1 #calculamos la diferencia entre los valores x1 y x2, la FWHM
    return abs(delta),x1_i,x2_i #La función arroja el valor calculado de FWHM así como los índices de x1 y x2 