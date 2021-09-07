#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 12:41:30 2021


########## LEs fonctions CDFT ###############

@author: Mathis Deronzier
"""


""" Nous allons ici écire l'algortihme CDF-t 
-On suppose que les données que nous avons sont sous format numpy
-Nous allons faire plusieurs méthodes de CDFT

"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy.random as rng
from scipy import stats
from math import floor



"""On va écire une fonction CDF-t qui nous permmettra de faire la 
transformation simplement """

def Cum_Dist_Fun(X):
#On suppose X un vecteur de tirages aléatoires la fonction renvoie un tableau 
    n=len(X)
    Y=X
    Y.sort()
    m=Y[0]
    s=0
    Z=[]
    for i in range(n):
        if Y[i]>m:
            Z.append([m,s/n])
            m=Y[i]
        s+=1
        if i==n-1:
            Z.append([m,s/n])
    return np.array(Z)
            
def find(x,X):
#On suppose que la liste X est triée et l'on cherche la position dans la liste 
#telle que l'élément soit dedans
    a=0
    b=len(X)
    while b-a>1:
        m=(a+b)//2
        if x<X[m]:
            b=m
        else:
            a=m
    return a

def find_t(x,a,b):
    #On a a<=x<b on renvoie t vérifiant x=ta+(1-t)b
    return (b-x)/(b-a)
    

def QQ(X,Y,Xf):
#Cette fonction prend en entrée Les listes X, Y et la liste Xf à downscaler
    F_X=Cum_Dist_Fun(X)
    F_Y=Cum_Dist_Fun(Y)
    Y_pred=[]
    nx=F_X.shape[0]
    ny=F_Y.shape[0]
    for x in Xf:
        rx=find(x,F_X[:,0])
        if rx<(nx-1):
            tx=find_t(x, F_X[rx,0],F_X[rx+1,0])
            px=tx*F_X[rx,1]+(1-tx)*F_X[rx+1,1]
            ry=find(px,F_Y[:,1])
            ty=find_t(px, F_Y[ry,1],F_Y[ry+1,1])
            py=ty*F_Y[ry,0]+(1-ty)*F_Y[ry+1,0]
        else:
            py=F_Y[ny-1,0]
        Y_pred.append(py)
    return np.array(Y_pred)
    
def plot_cdf(X):
#Fonction affichant la CDF de l'échantillon X
    S=Cum_Dist_Fun(X)
    #fig=plt.figure()
    plt.plot(S[:,0],S[:,1])
    plt.xlabel("X")
    plt.ylabel("F_X")
    plt.title("Fonction de répartition empirique de X")
    #fig.savefig("/home/mathis/stage/latex/images/classification_deb_prr")
    plt.show()
    

def CDFt_classic(X,Y,Xf):
#On choisit la fonction CDF-t classique
    return QQ(X,Y)*(X.argmax())/(Y.argmax())



    
