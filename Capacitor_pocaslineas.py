# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 23:51:04 2020

@author: Esteban Velásquez
"""
#Simulación potencial dos placas circulares de radio 13cm paralelas
import scipy.constants as const
import math
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
b=1 #este numero es la distancia entre las placas en mm
dist=0.1#se utiliza 0.1 mm como distancia entre cada punto
d=np.zeros([40+b*10+2,2601+600]) #esta es la matriz que representa el valor de cada punto en el espacio siendo el tamaño de las placas de ancho y la distancia de altura
d[20,300:2900]=50
d[20+10*b+1,300:2900]=-50 #Se introducen los valores del potencial inicial
k= np.zeros([40+b*10+2,2601+600])# nueva matriz donde se van a almacenar la nueva iteración
k[20,300:2900]=50
k[20+10*b+1,300:2900]=-50
#for i in range(0,6)
for times in range(0,1500):#acá se inicializa los ciclos que se van a usar
  for i in range(1, 40+b*10+1):#se recorre la matriz en el alto
    for j in range(1, 3200):# se recorre la matriz en el ancho
      k[i,j]=0.25*(d[i+1,j]+d[i-1,j]+d[i,j+1]+d[i,j-1])
  d=deepcopy(k)
  d[20,300:2900]=50
  d[20+10*b+1,300:2900]=-50
# Ya se procede a inicializar el campo electrico y resolverlo
Ey=np.zeros([40+(b*10)+2,3200])
for i in range(0,40+b*10+1):#recorre en la altura
    for j in range(1, 3200):#recorre en el ancho
        Ey[i,j]=-(d[i+1,j]-d[i,j])/0.0001
Ex=np.zeros([40+(b*10)+2,3200])
for i in range(0,40+b*10+1):#recorre en la altura
    for j in range(0, 3199):#recorre en el ancho
        Ex[i,j]=-(d[i,j+1]-d[i,j])/0.0001  
Etot=np.hypot(Ex, Ey) 
#otra forma de calcularlo es:
Vgrad=np.gradient(d)
Exgrad=-Vgrad[1]*10**4
Eygrad=-Vgrad[0]*10**4

   
# ahora con geometria procedemos a hallar el area de la porción cicular de cada pedazo de linea tomando en cuenta el efecto de borde       

Qt=0
Areatot=0
for i in range(301,1601):#tomo solo la primera mitad del circulo y por simetria es la mitad me voy moviendo de izquierda a derecha
    R=1300
    D=1300-(i-300)# denota el radio en el que esta en el momento de utilizar el area respectiva
    Theta=(180/math.pi)*math.acos(D/1300)# denota el angulo en el que me estoy ubicando en grados
    Theta2=(180/math.pi)*math.acos((D+1)/1300)# denota el angulo de la proxima posicion en grados
    Area=((R**2*math.pi*(Theta*2)/360)-(R*D*math.sin(Theta*math.pi/180)))-(((R)**2*math.pi*(Theta2*2)/360)-((D+1)*R*math.sin(Theta2*math.pi/180)))
    #El area denota el pedacito de area de ese 0.1 mm de ancho
    Qt=Qt+Ey[21,i]*Area/(10000**2)# ac'a se aproxima que el camp en el eje perpendicular al ancho del condensador es el mismo
    Areatot=Areatot+Area# Esto define el are total utilizada
#Ahora calculamos la capacitancia con esta aproximación de tomar en cuenta el efecto borde dado 2D
cap=2*Qt*const.epsilon_0/100
Areasobred=math.pi*0.13**2/(0.001*b)
Err=100*(cap-Areasobred*const.epsilon_0)/(Areasobred*const.epsilon_0)
Areatot=2*Areatot/(10000**2)#para verificar que si recoriio toda el area del circulo

#sin tomar en cuenta el efecto de borde tenemos:
Qt1=math.pi*0.13**2*Ey[21,1600]
cap1=Qt1*const.epsilon_0/100
Err1=100*(cap1-Areasobred*const.epsilon_0)/(Areasobred*const.epsilon_0)

#procedemos ahora a graficar el potencial
plt.matshow(Etot, aspect= 'auto',cmap = 'gist_ncar', interpolation='lanczos')
plt.xlabel('$x(100um)$')
plt.ylabel('$y(100um)$')
cbar = plt.colorbar()
cbar.set_label('Electric field(V/m)')

#procedo ahora a graficar el campo electrico vectorialmente
ny, nx = 40+(b*10)+2,3200
x = np.linspace(0, 3200, nx)
y = np.linspace(40+(b*10)+2,0 , ny)
X, Y = np.meshgrid(x, y)
color = 2 * np.log(np.hypot(Ex, Ey))
fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.streamplot(X, Y, Ex, Ey, color=color, linewidth=1, cmap=plt.cm.inferno,
              density=5, arrowstyle='->', arrowsize=1.5)
ax.set_xlabel('$x(100um)$')
ax.set_ylabel('$y(100um)$')
ax.set_xlim(0,3200)
ax.set_ylim(40+(b*10)+2,0)
fig.colorbar(im.lines)
ax.set_title('log Electric field')
plt.show()

#ahora grafico vectorialmente con la otra opción:
ny1, nx1 = 40+(b*10)+2,3201
x1 = np.linspace(0, 3201, nx1)
y1 = np.linspace(40+(b*10)+2,0 , ny1)
X1, Y1 = np.meshgrid(x1, y1)
color2 = 2 * np.log(np.hypot(Exgrad, Eygrad))
fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.streamplot(X1, Y1, Exgrad, Eygrad, color=color2, linewidth=1, cmap='gist_ncar',
              density=5, arrowstyle='->', arrowsize=1.5)
ax.set_xlabel('$x(100um)$')
ax.set_ylabel('$y(100um)$')
ax.set_xlim(0,3200)
ax.set_ylim(40+(b*10)+2,0)
fig.colorbar(im.lines)
ax.set_title('log Electric field')
plt.show()

#comprobación con el metodo de gradiente
Qt2=math.pi*0.13**2*Eygrad[21,1600]
cap2=Qt2*const.epsilon_0/100
Err2=100*(cap2-Areasobred*const.epsilon_0)/(Areasobred*const.epsilon_0)

#por ultimo muestra los errores
print(cap)
print(Err)

