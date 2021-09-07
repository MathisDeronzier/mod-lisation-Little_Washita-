#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 19:26:26 2021

@author: mathis
"""
import os
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d



os.chdir("/home/mathis/stage/py_codes/kernels2/")

#Datas=np.arange(0,11**3*2).reshape((11,11,11,2))
#K=kernel(Datas)

files= os.listdir()

mois=['janvier', 'février', 'mars', 'avril', 'mai', 'juin', 'juillet','août', 'septembre', 'octobre', 'novembre', 'décembre']

def plot_surface(file, Title="", grids=False, X=[], Y=[]):
    #récupération de la matrice de covariance
    K=pd.read_csv(file).to_numpy()[:,1:]
    #K=K/max(K)
    #Création de la figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    #On veut faire des pas de 32km par 32km
    dim=K.shape
    if grids:
        ax.plot_surface(X, Y, K,rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
        #fig.title("Calcule du noyau de covariance empirique pour "+file)
        ax.set_title(Title)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('K(X,Y)')
        fig.savefig("/home/mathis/stage/latex/images/"+file[:len(file)-4])
        plt.show()
    else:
        x = y = np.array([32*k for k in range(dim[0])])-32*(dim[0]//2)
        X, Y = np.meshgrid(x, y)
            
        ax.plot_surface(X, Y, K,rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        #fig.title("Calcule du noyau de covariance empirique pour "+file)
        ax.set_title(   Title)
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('K(X,Y)')
        ax.set_zlim3d(0,40)
        fig.savefig("/home/mathis/stage/latex/images/"+file[:len(file)-4])
        plt.show()


#for i in range(1,13):
#    plot_surface('kernel_evap_m'+ str(i)+'.csv', "Noyau de covariance spatial de l'évapotranspiration en "+mois[i-1])

for i in range(1,13):
    plot_surface('kernel_precip_m'+ str(i)+'.csv', "Noyau de covariance spatial des précipitations en "+mois[i-1])
################ Premier plot ##################

#On va maintenant partitionner les points selon nos droite de régression 

