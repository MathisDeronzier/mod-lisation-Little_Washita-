#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 22:09:19 2021

@author: mathis
"""
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import random
from mpl_toolkits import mplot3d



############## Calcul du kernel d'une fonction aléatoire #####################

def kernel(Datas):
#Attention ce calclule est très mal optimisé et par conséquent s'upscale très mal aux
#grosses données
    dim=Datas.shape
    M=dim[0]
    N=dim[1]
    T=dim[2]
    v=dim[3]
    K=np.zeros((M,N,v))
    for e in range(v):
        for m in range(-(M//2)+1,M//2):
            for n in range(N//2):
                T_1=0
                T_2=0
                H=T*(min(M,M+m)-max(0, m))*(N-n)
                for t in range(T):
                    for i in range(max(0, m),min(M,M+m)):
                        for j in range(n,N):
                            T_1+=Datas[i,j,t,e]
                            T_2+=Datas[i-m,j-n,t,e]
                T_1=T_1/H
                T_2=T_2/H
                #print("fin de la première boucle:\n T_1=",T_1, ", T_2=",T_2)        
                indice_i=M//2+m
                indice_j=N//2+n
                for t in range(T):
                    for i in range(max(0, m),min(M,M+m)):
                        for j in range(n,N):
                            K[indice_i,indice_j,e]+=(Datas[i,j,t,e]-T_1)*(Datas[i-m,j-n,t,e]-T_2)
                K[indice_i,indice_j,e]=K[indice_i,indice_j,e]/(H-1)
                K[indice_i-2*m,indice_j-2*n,e]=K[indice_i,indice_j,e]
                #print("fin  de la deuxième boucle: \n K["+str(indice_i)+","+str(indice_j)+","+ str(v)+"]="+str(K[indice_i,indice_j,e])+
                #      ", K["+str(indice_i-2*m)+","+str(indice_j-2*m)+","+ str(v)+"]="+str(K[indice_i-2*m,indice_j-2*n,e]))
    return K


Datas=np.arange(0,11**3*2).reshape((11,11,11,2))
K=kernel(Datas)


################ Premier plot ####################

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#On veut faire des pas de 32km par 32km
dim=K.shape

x = y = np.array([32*k for k in range(dim[0])])-32*(dim[0]//2)
X, Y = np.meshgrid(x, y)

ax.plot_surface(X, Y, K[:,:,1],rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('K(X,Y) Label')

plt.show()

