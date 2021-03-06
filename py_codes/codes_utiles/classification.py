#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 20:15:15 2021

@author: mathis
"""
###########   Classification des points  ###############
########### Optimisation sous contrainte ###############

#On commence par définir les fonctions qui seront utiles

import numpy as np
import pandas as pd
import os
from scipy.optimize import minimize
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import random
from mpl_toolkits import mplot3d

########### Quelques fonctions utiles ##############

def dist(u,v,x):
    #La distance d'un point à la droite définie par le couple (u,v) voir rapport
    y=(x-v)
    return y.dot(y)-(y.dot(u))**2/(u.dot(u))

    
def f(U,X):
    #La fonction que l'on va chercher à optimiser
    tot=0   
    for x in X:
        tot+=min(dist(U[0],U[1],x), dist(U[2],U[3],x))
    return tot

def intersect(X,Y):
    #renvoie l'intersection de deux listes triées par ordre croissant
    n=len(X)
    m=len(Y)
    a=0
    b=0
    l=[]
    while a<n and b<m:
        if X[a]==Y[b]:
            l.append(X[a])
            a+=1
            b+=1
        else:
            if X[a]>Y[b]:
                b+=1
            else:
                a+=1
    return np.array(l)

def dif_dist(X,u1,v1,u2,v2):
    distances=np.zeros((X.shape[0],1))
    for i in range(X.shape[0]):
        distances[i]=dist(u1,v1,X[i])-dist(u2,v2,X[i])
    return distances


#### Test de la fonction d'optimisation #######
os.chdir('/home/mathis/stage/R_codes/series')
Dates=pd.read_csv("dates_1979_2014.csv", index_col=0).to_numpy()
dates=np.where(Dates[:,0]>=1996)[0]

etp_dsI=pd.read_csv("etp_dsI.csv", index_col=0).to_numpy()
n=len(etp_dsI)

ETP_series=np.zeros((n,11))
pr_series=np.zeros((n,11))
prr_series=np.zeros((n,11))
deb_series=np.zeros((n,11))
etr_series=np.zeros((n,11))

#Définition des séries de temps
times=np.zeros((n,11))
times[:,0]= np.where(pd.read_csv('simu_ref.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,1]= np.where( pd.read_csv('simu_ds1.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,2]= np.where( pd.read_csv('simu_ds2.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,3]= np.where( pd.read_csv('simu_ds3.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,4]= np.where( pd.read_csv('simu_ds4.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,5]= np.where( pd.read_csv('simu_dsI.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,6]= np.where( pd.read_csv('simu_deg1.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,7]= np.where( pd.read_csv('simu_deg2.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,8]= np.where( pd.read_csv('simu_deg3.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,9]= np.where( pd.read_csv('simu_deg4.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times[:,10]= np.where( pd.read_csv('simu_IPSL.csv', index_col=0).to_numpy()[:,0]%86400 ==0)[0]
times=times.astype(int)

for i in range(6):
    if i==0:
        ETP_series[:,0]=pd.read_csv('etp_ref.csv', index_col=0).to_numpy()[:,0]
        pr_series[:,0]=pd.read_csv('pr_ref.csv', index_col=0).to_numpy()[:,0]
        deb_series[:,0]=-pd.read_csv('simu_ref.csv', index_col=0).to_numpy()[times[:,0],1]
        prr_series[:,0]=pd.read_csv('simu_ref.csv', index_col=0).to_numpy()[times[:,0],2]
        etr_series[:,0]=pd.read_csv('simu_ref.csv', index_col=0).to_numpy()[times[:,0],17]
    if i<5 and i>0:
        ETP_series[:,i]=pd.read_csv('etp_ds'+str(i)+'.csv', index_col=0).to_numpy()[:,0]
        pr_series[:,i]=pd.read_csv('pr_ds'+str(i)+'.csv', index_col=0).to_numpy()[:,0]
        deb_series[:,i]=-pd.read_csv('simu_ds'+str(i)+'.csv', index_col=0).to_numpy()[times[:,i],1]
        prr_series[:,i]=pd.read_csv('simu_ds'+str(i)+'.csv', index_col=0).to_numpy()[times[:,i],2]
        etr_series[:,i]=pd.read_csv('simu_ds'+str(i)+'.csv', index_col=0).to_numpy()[times[:,i],17]
        
        ETP_series[:,i+5]=pd.read_csv('etp_deg'+str(i)+'.csv', index_col=0).to_numpy()[:,0]
        pr_series[:,i+5]=pd.read_csv('pr_deg'+str(i)+'.csv', index_col=0).to_numpy()[:,0]
        deb_series[:,i+5]=-pd.read_csv('simu_deg'+str(i)+'.csv', index_col=0).to_numpy()[times[:,i],1]
        prr_series[:,i+5]=pd.read_csv('simu_deg'+str(i)+'.csv', index_col=0).to_numpy()[times[:,i],2]
        etr_series[:,i+5]=pd.read_csv('simu_deg'+str(i)+'.csv', index_col=0).to_numpy()[times[:,i],17]
    if i==5:
        ETP_series[:,i]=pd.read_csv('etp_dsI.csv', index_col=0).to_numpy()[:,0]
        pr_series[:,i]=pd.read_csv('pr_dsI.csv', index_col=0).to_numpy()[:,0]
        deb_series[:,i]=-pd.read_csv('simu_dsI.csv', index_col=0).to_numpy()[times[:,i],1]
        prr_series[:,i]=pd.read_csv('simu_dsI.csv', index_col=0).to_numpy()[times[:,i],2]
        etr_series[:,i]=pd.read_csv('simu_dsI.csv', index_col=0).to_numpy()[times[:,i],17]
        
        ETP_series[:,i+5]=pd.read_csv('etp_IPSL.csv', index_col=0).to_numpy()[:,0]
        pr_series[:,i+5]=pd.read_csv('pr_IPSL.csv', index_col=0).to_numpy()[:,0]
        deb_series[:,i+5]=-pd.read_csv('simu_IPSL.csv', index_col=0).to_numpy()[times[:,i],1]
        prr_series[:,i+5]=pd.read_csv('simu_IPSL.csv', index_col=0).to_numpy()[times[:,i],2]
        etr_series[:,i+5]=pd.read_csv('simu_IPSL.csv', index_col=0).to_numpy()[times[:,i],17]
    

#Il faut décaler la série temporelle pour qu'il y ait une correlation  


## Analyse simple des deltas
        
# =============================================================================
# for decal in range(1):
#     points_prr=np.where(prr_series[0:n-decal,0]>1e-7)[0]
#     points_deb=np.where(deb_series[decal:n,0]>1e-7)[0]
#     points=intersect(points_prr, points_deb)
#     for i in range(1,5):
#         delta_prr=(prr_series[points,i]-prr_series[points,0])/prr_series[points,0]
#         delta_deb=(deb_series[points+decal,i]-deb_series[points+decal,0])/deb_series[points+decal,0]
#         points_inter_prr=np.where(delta_prr<2)[0]
#         points_inter_deb=np.where(delta_deb<2)[0]
#         points_inter=intersect(points_inter_prr, points_inter_deb)
#         delta_prr=delta_prr[points_inter]
#         delta_deb=delta_deb[points_inter]
#         X=np.column_stack((delta_prr, delta_deb))
#         U_0=np.array([[1,0,0,0,1,1,0,0]])
#         def F(U):
#             return f(U.reshape((4,2)),X)
#         #On ajoute les contraintes pour que ca soit sur la sphre unit 
#         #cons = ({'type': 'eq', 'fun': lambda x:  x[0]**2+x[1]**2-1,
#         #'jac' : lambda x: np.array([2*x[0], 2*x[1],0,0,0,0,0,0])},
#         #{'type': 'eq', 'fun': lambda x: x[4]**2+x[5]**2-1,
#         #'jac' : lambda x: np.array([0,0,0,0,2*x[4],2*x[5],0,0])})
#         #sol=minimize(F, U_0, constraints=cons).x
#         sol=minimize(F, U_0).x
#         print( "le vecteur de minimisation pour decal="+str(decal)+ " et i="+str(i)+" est ", sol )
#         def D1(x):
#             return x*sol[1]/sol[0] -sol[1]*sol[2]/sol[0]+sol[3]
#         def D2(x):
#             return x*sol[5]/sol[4] -sol[5]*sol[6]/sol[4]+sol[7]
#         
#         fig = plt.figure()
#         ax=fig.add_axes([0,0,1,1])
#         #ax.set_xlim(-1.5,1.5)
#         #ax.set_ylim(-1.5,1.5)
#         x=np.array([min(delta_prr),max(delta_prr)])
#         y1=np.array([D1(x[0]),D1(x[1])])
#         y2=np.array([D2(x[0]),D2(x[1])])        
#         ax.plot(x, y1, label="D1", color="blue") 
#         ax.plot(x, y2, label="D2", color="red") 
#         ax.scatter(delta_prr, delta_deb, color="black")
#         ax.set_xlabel("delta précipitations")
#         ax.set_ylabel("delta débits")
#         ax.set_title("tracé des Deltas pr/débits dwnsc("+ str((2*i-1)^2)+")")
#         fig.savefig("/home/mathis/stage/latex/images/deb_prr_ds"+str(i))
#         plt.show()
# =============================================================================


for decal in range(-1,0):
    if decal<0:
        points_pr=np.where(pr_series[-decal:n,0]>1e-4)[0]
        points_deb=np.where(deb_series[0:n+decal,0]>1e-7)[0]
        points_prr=np.where(deb_series[0:n+decal,0]!=0)[0]
        points=intersect(intersect(points_pr, points_deb),points_prr)
    else:
        points_pr=np.where(pr_series[0:n-decal,0]>1e-4)[0]
        points_deb=np.where(deb_series[decal:n,0]>1e-7)[0]
        points=intersect(points_pr, points_deb)
        points_prr=np.where(deb_series[0:n+decal,0]!=0)[0]
        points=intersect(intersect(points_pr, points_deb),points_prr)
    for i in range(4,5):
        if decal<0:
            delta_pr=(pr_series[points-decal,i]-pr_series[points-decal,0])/pr_series[points-decal,0]
            delta_prr=(prr_series[points,i]-prr_series[points,0])/prr_series[points,0]
            delta_deb=(deb_series[points,i]-deb_series[points,0])/deb_series[points,0]
        else:
            delta_pr=(pr_series[points,i]-pr_series[points,0])/pr_series[points,0]
            delta_deb=(deb_series[points+decal,i]-deb_series[points+decal,0])/deb_series[points+decal,0]
        
        points_inter_pr=np.where(delta_pr<3)[0]
        points_inter_deb=np.where(delta_deb<3)[0]
        points_inter=intersect(points_inter_pr, points_inter_deb)
        delta_pr=delta_pr[points_inter]
        delta_deb=delta_deb[points_inter]
        delta_prr=delta_prr[points_inter]
        X=np.column_stack((delta_pr, delta_deb))
        U_0=np.array([[1,0,0,0,1,1,0,0]])
        def F(U):
            return f(U.reshape((4,2)),X)
        #On ajoute les contraintes pour que ca soit sur la sphre unit 
        #cons = ({'type': 'eq', 'fun': lambda x:  x[0]**2+x[1]**2-1,
        #'jac' : lambda x: np.array([2*x[0], 2*x[1],0,0,0,0,0,0])},
        #{'type': 'eq', 'fun': lambda x: x[4]**2+x[5]**2-1,
        #'jac' : lambda x: np.array([0,0,0,0,2*x[4],2*x[5],0,0])})
        #sol=minimize(F, U_0, constraints=cons).x
        sol=minimize(F, U_0).x
        frac=len(points_inter)/n
        print( "le vecteur de minimisation pour decal="+str(decal)+ " et i="+str(i)+" est ", sol ,
              "la faction d'échantillons étudiée est ="+ str(frac))
        def D1(x):
            return x*sol[1]/sol[0] -sol[1]*sol[2]/sol[0]+sol[3]
        def D2(x):
            return x*sol[5]/sol[4] -sol[5]*sol[6]/sol[4]+sol[7]
        fig=plt.figure()
        dist=dif_dist(X,sol[0:2],sol[2:4],sol[4:6],sol[6:8])
        blue_point=np.where(dist>0)[0]
        red_point=np.where(dist<0)[0]
        plt.scatter(delta_pr[blue_point], delta_deb[blue_point], color="red")
        plt.scatter(delta_pr[red_point], delta_deb[red_point], color="blue")
        x=np.array([min(delta_pr),max(delta_pr)])
        y1=np.array([D1(x[0]),D1(x[1])])
        y2=np.array([D2(x[0]),D2(x[1])])        
        plt.plot(x, y1, label="D1", color="blue") 
        plt.plot(x, y2, label="D2", color="red") 
        plt.xlabel("delta précipitations")
        plt.ylabel("delta débits")
        plt.xlim((-1.1, 2.5))
        plt.ylim((-1.1, 2.5))
        #plt.title("tracé des Deltas précipitation réelles/débits dwnsc("+ str((2*i-1)**2)+") decal="+str(decal))
        fig.savefig("/home/mathis/stage/latex/images/classification_deb_pr")
        plt.show()
        ################### On affiche maintenant les points qui nous intéressent
        fig=plt.figure()
        plt.scatter(delta_prr[blue_point], delta_deb[blue_point], color="red")
        plt.scatter(delta_prr[red_point], delta_deb[red_point], color="blue")
        #x=np.array([min(delta_pr),max(delta_pr)])
        #y1=np.array([D1(x[0]),D1(x[1])])
        #y2=np.array([D2(x[0]),D2(x[1])])        
        #plt.plot(x, y1, label="D1", color="blue") 
        #plt.plot(x, y2, label="D2", color="red") 
        plt.xlabel("delta précipitations réelles")
        plt.ylabel("delta débits")
        plt.xlim((-1.1, 2.5))
        plt.ylim((-1.1, 2.5))
        #plt.title("tracé des Deltas précipitation/débits dwnsc("+ str((2*i-1)**2)+") decal="+str(decal))
        fig.savefig("/home/mathis/stage/latex/images/classification_deb_prr")
        plt.show()
         

#On a décalé les points de 

#U=uzawa(U_0, X, 1e-5)


# Il semblerait qu'ils soient tous indépendants, il faut alors vérifier les données 
    


    
