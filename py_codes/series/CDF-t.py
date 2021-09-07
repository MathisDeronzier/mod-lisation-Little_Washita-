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
#de valeurs per
    n=len(X)
    Y=X.sort()
    m=Y[0]
    s=0
    Z=[]
    for i in Y:
        if i>m:
            Z.append([m,s/n])
            s=1
            m=i
        else:
            s+=1
    return np.array(Z)
            
def find(x,X):
#On suppose que la liste X est triée et l'on cherche la position dans la liste 
#telle que l'élément soit dedans
    a=0
    b=len(X)
    while b-a>1:
        m=floor((a+b)/2)
        if x<=X[m]:
            a=m
        else:
            m=b
    return a

def QQ(X,Y,Xf):
#Cette fonction prend en entrée deux fonctions de distribution F et G et 
#Sort H la tranformation du quantile-quantile 
#Comment calculer la fonction de transfert?
    F_X=Cum_Dist_Fun(X)
    F_Y=Cum_Dist_Fun(Y)
    Y_pred=[]
    for x in Xf:
        pos_X=find(x,F_X[0,:])
        rg_X=F_Y[1,pos_X+1]*(pos_X-F_X[0,pos_X])/(F_X[0,pos_X+1]-F_X[0,pos_X])+F_X[1,pos_X]*(F_X[0,pos_X+1]-pos_X)/(F_X[0,pos_X+1]-F_X[0,pos_X])
        pos_Y=find(rg_X,F_Y[1,:])
        Y_pred.append(F_Y[0,pos_Y+1]*(rg_X-F_Y[1,pos_Y])/(F_Y[1,pos_Y+1]-F_Y[1,pos_Y])+
                     F_Y[0,pos_Y]*(F_Y[1,pos_Y+1]-rg_X)/(F_Y[1,pos_Y+1]-F_Y[1,pos_Y]))
    return np.array(Y_pred)
    

def CDFt_classic(X,Y,Xf):
    return QQ(X,Y)*(X.argmax())/(Y.argmax())



    
