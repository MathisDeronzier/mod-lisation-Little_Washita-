#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 20:15:15 2021

@author: mathis
"""

########### Optimisation sous contrainte ###############

import numpy as np
from scipy.optimize import minimize

def dist(u,v,x):
    y=((x-v)-u*(x-v).dot(u))
    return y.dot(y)

X=np.arange(100).reshape((2,50))
    
def f(U):
    tot=0
    for x in X:
        tot+=min(dist(U[0],U[1],x), dist(U[2],U[3],x))
    return tot

def grad(f,U,eps):
    dim=U.shape()
    grad=np.zeros((dim[0],dim[1]))
    for i in range(dim[0]):
        for j in range(dim[1]):
            U_1=U
            U_1[i,j]=U_1[i,j]+eps
            grad[i,j]=(f(U_1)-f(U))/eps
    return grad

def armijo(f,U,eps):
    alpha=1
    G=grad(f,U,eps)
    while f(U-alpha*G)>-alpha*G.dot(G)/2:
        alpha/=2
    return U-alpha*G

def proj(U):
    U[0]=U[0]/(U[0].dot(U[0]))
    U[2]=U[2]/(U[2].dot(U[2]))
    
def uzawa(f,U_0,eps):
    #On pose une condition d'arrêt
    U_1=proj(armijo(f,U_0,eps))
    while f(U_1)-f(U_1)>1000*eps:
        U_0, U_1 = U_1, proj(armijo(f,U_0,eps))
    return U_1
        
    
#### Test de la fonction d'optimisation #######
    
import pandas as pd
import os

os.chdir("/home/mathis/stage/R_codes/series/")


Dates=pd.read_csv("dates_1979_2014.csv", index_col=0).to_numpy()
dates=np.where(Dates[:,0]>=1996)[0]

etp_dsI=pd.read_csv("etp_dsI.csv", index_col=0).to_numpy()
n=len(etp_dsI)

ETP_series=np.zeros((n,11))
pr_series=np.zeros((n,11))
deb_series=np.zeros((n,11))
etr_series=np.zeros((n,11))

    
times_debtimes_ref<- np.where( simulations_ref[,1]%%86400 ==0)
times_ds1= np.where( simulations_ds1[:,0]%%86400 ==0)
times_ds2= np.where( simulations_ds2[:,0]%%86400 ==0)
times_ds3= np.where( simulations_ds3[:,0]%%86400 ==0)
times_ds4= np.where( simulations_ds4[:,0]%%86400 ==0)
times_dsI= np.where( simulations_dsI[:,0]%%86400 ==0)
times_deg1= np.where( simulations_deg1[:,0]%%86400 ==0)
times_deg2= np.where( simulations_deg2[:,0]%%86400 ==0)
times_deg3= np.where( simulations_deg3[:,0]%%86400 ==0)
times_deg4= np.where( simulations_deg4[:,0]%%86400 ==0)
times_IPSL= np.where( simulations_IPSL[:,0]%%86400 ==0)

for i in range(5):
    if i==0:
        ETP_series[:,i]=pd.read_csv('evap_ref_1979_2014.csv', index_col=0).to_numpy()[dates]
        pr_series=[:,i]=pd.read_csv('precip_ref_1979_2014.csv', index_col=0).to_numpy()[dates]     
        deb_series=[:,i]=
        etr_series=[:,i]
        
            
    
    
    
    
    
    