#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 14:43:49 2021

@author: mderon
"""
"""On va ici écrire un programme récupérant les données NARR pour chaque année
et crée un gros fichier css contenant les données de toutes les années"""

import os #Bib pour les commandes système 
import xarray as xr #Bib pour manipuler des datasets
import argparse #Bib pour entrer des paramêtres dans le code
import numpy as np 
from scipy.ndimage import convolve #Outil pour la convolution de matrice
import pandas as pd #Lib pour créer et manipuler des Dataframes


#La lat et lon qu'on veut récupérer
latitudes=[84,98]
longitudes=[183,207]


nb_j=[31,28,31,30,31,30,31,31,30,31,30,31] #Liste du nombre de jours par mois
nb_j_bis=[31,29,31,30,31,30,31,31,30,31,30,31] #Liste du nombre de jours par mois: année bissextile

################## Programmes pour récupérer les données ######################


def get_data(year_start,year_stop, lat = latitudes, lon = longitudes):
    #Fonction qui récupère les données (Datas) des années year_start à year_stop prendre (1979,2003) il y a des bugs après
    #Datas est un vecteur multidimensionel de la forme: Data(lat, lon, date, variable)
    #Sortie (Datas, vecteur des dates, les latitudes, les longitudes)
    d_tot=0
    for year in range(year_start,year_stop+1):
        os.chdir("/home/hydrogeol2/mderon/NARR_1979_2014/"+str(year))
        y_m=year_miss2(year,True)
        if y_m!=[]:
            print("Toutes les données pour l'année "+str(year)+" n'ont pas été telechargées, liste des dates manquantes:")
            print(y_m)
        d_tot+=bissextile(year)+365-len(y_m) #Le nombre de jours dans l'année
    Datas=np.zeros((lat[1]-lat[0],lon[1]-lon[0],d_tot,2))
    dates=np.zeros((d_tot,3))
    brebis_galeuses=[]
    i=0 #Compteur pour les dates
    for year in range(year_start,year_stop+1):
        print("year =",year)
        os.chdir("/home/hydrogeol2/mderon/NARR_1979_2014/"+str(year))
        l=os.listdir()
        nb_d=nb_j
        if bissextile(year):
            nb_d=nb_j_bis
        for month in range(1,13):
            for day in range(1, nb_d[month-1]+1):
                h=0
                for hour in range(8):
                    file_name=name_file(year,month,day,hour)
                    if file_name in l:
                        h+=1
                        try:
                            ds=xr.open_dataset(file_name)
                            Datas[:,:,i,0]=Datas[:,:,i,0]+ds["A_PCP_221_SFC_acc3h"].data[lat[0]:lat[1],lon[0]:lon[1]]
                            Datas[:,:,i,1]=Datas[:,:,i,1]+ds["PEVAP_221_SFC_acc3h"].data[lat[0]:lat[1],lon[0]:lon[1]]
                        except:
                            print("the document can't be open")
                            h-=1
                if h>0:
                    Datas[:,:,i,:]=Datas[:,:,i,:]*(8/h)
                    dates[i,0],dates[i,1],dates[i,2]=year,month,day
                    i+=1
                else:
                    brebis_galeuses.append([year,month,day])
        grid_lat = ds["gridlat_221"].data[lat[0]:lat[1],lon[0]:lon[1]]
        grid_lon = ds["gridlon_221"].data[lat[0]:lat[1],lon[0]:lon[1]]
    return(Datas,dates,grid_lat,grid_lon,np.array(brebis_galeuses))


def name_dates(year,month,day):
#Fonction renvoyant un format usuel pour les dates
    str_y=str(year)
    str_m="0"+str(month)
    if len(str(month))==2:
        str_m=str(month)
    str_d=str(day)
    if len(str(day))==1:
        str_d="0"+str(day)
    return (str_y+"-"+str_m+"-"+str_d)

###################### Ici les programmes de dégradation ######################


def axial_sym(conv):
    #symetrise une matrice: utile pour les convolutions
    dim=conv.shape
    sym=np.zeros(dim)
    for i in range(dim[0]):
        for j in range(dim[1]):
            sym[i,j]=conv[dim[0]-1-i,dim[1]-1-j]
    return sym

def agglomerate(Datas, conv):
    #fonction faisant le produit de convolution de Datas notre vecteur multidimensionnel et 
    #d'une matrice de convolution qu'on chosit généralement uniforme.
    dim_Datas=Datas.shape
    conv_Datas=np.zeros(dim_Datas)
    sym=axial_sym(conv)
    for i in range(dim_Datas[2]):
        for j in range(dim_Datas[3]):
             conv_Datas[:,:,i,j]=convolve(Datas[:,:,i,j], sym, mode='constant', cval=0.00)
    return conv_Datas
            
def agglomerated(conv_Datas, conv, start_point):
    #Renvoie la convolution donnée selon 1 point choisi (i,j) de coordonnées (lat[i],lat[j]), pour avoir différents maillages
    dim_Datas=conv_Datas.shape
    dim_conv=conv.shape
    width=dim_conv[1]//2
    height=dim_conv[0]//2
    if start_point[0]>(dim_Datas[0]-height) or  start_point[0]<height or start_point[1]>(dim_Datas[1]-width) or start_point[1]<width:
        print("erreur le point choisi ne peut être approximé")
        return(0)
    else:
        a = start_point[0]//dim_conv[0] + (dim_Datas[0]-start_point[0])//dim_conv[0]
        b = start_point[1]//dim_conv[1] + (dim_Datas[1]-start_point[1])//dim_conv[1]
        c = start_point[0]%dim_conv[0]
        d = start_point[1]%dim_conv[1]
        agglomerate_data = np.zeros((a,b,dim_Datas))
        for t in range(dim_Datas[2]):
            for i in range(a):
                for j in range(b):
                    agglomerate_data[i,j,t,0] = conv_Datas[c+dim_conv[0]*i,d+dim_conv[1]*j,t,0]
                    agglomerate_data[i,j,t,1] = conv_Datas[c+dim_conv[0]*i,d+dim_conv[1]*j,t,1]
        return agglomerate_data

def agglomerate_unif(Datas, n, lat, lon):
    dims=Datas.shape
    res=np.zeros(dims[2])
    for i in range(dims[2]): 
        sum=0
        for k in range(-n,n+1):
            for l in range (-n,n+1):
                sum+=Datas[lat+k,lon+l,i]
        res[i]=sum/((2*n+1)**2)
    return res


###################### Ici les programmes peu importants ######################

def name_file(y,m,d,h,option=False):
    #Récupère les noms des fichiers dans le répertoire
    str_y=str(y)
    str_m="0"+str(m)
    if m>9:
        str_m=str(m)
    str_d=str(d)
    if d<10:
        str_d="0"+str(d)
    str_h=str(3*h)
    if 3*h<10:
        str_h="0"+str(3*h)
    if option:
        return (str_y+str_m+str_d)
    else:
        return (str_y+str_m+str_d+str_h+".nc")

def bissextile(year):
    #Vérifie si l'année est bissextile ou non.
    return ((year%400==0)or((year%4==0)and(year%100!=0)))


def year_miss(year,option=False):
    #Fonction renvoyant les données manquantes pour l'année year l'option permet de choisir le format sous lequel on renvoit les données
    nb_d=nb_j
    if bissextile(year): nb_d=nb_j_bis
    files_missing=[]
    files=os.listdir("/home/hydrogeol2/mderon/NARR_1979_2014/"+str(year))
    for m in range(1,13):
        for j in range(1, nb_d[m-1]+1):
            for h in range(8):
                if name_file(year,m,j,h) not in files:
                        files_missing.append(name_file(year,m,j,h,option))
    return files_missing

def year_miss2(year,option=False):
    #Fonction renvoyant les données manquantes pour l'année year l'option permet de choisir le format sous lequel on renvoit les données
    nb_d=nb_j
    if bissextile(year): nb_d=nb_j_bis
    files_missing=[]
    files=os.listdir("/home/hydrogeol2/mderon/NARR_1979_2014/"+str(year))
    for m in range(1,13):
        for j in range(1, nb_d[m-1]+1):
            t=0
            for h in range(8):
                if name_file(year,m,j,h) in files:
                    t+=1
            if t==0:
                    files_missing.append(name_file(year,m,j,0,option))
    return files_missing

def good_name(file):
    return file[0:14]=='merged_AWIP32.'

###### Ici les programmes pour trouver les points de grille correspondants ########

def inside(a,interval):
    #fonction prenant un point et un intervalle renvoyant si oui ou non le point
    #appartient à l'intervalle
    return (a<interval[1] and a>interval[0])

def min_max_grid(grid_lat,grid_lon):
    lat=[33,36.5]
    lon=[-104,-96]
    shape=grid_lon.shape
    min_max_lat=[shape[0],0]
    min_max_lon=[shape[1],0]
    for i in range (shape[0]):
        for j in range(shape[1]):
            if (inside(grid_lat[i,j],lat) and inside(grid_lon[i,j],lon)):
                if i<min_max_lat[0]:
                    min_max_lat[0]=i
                if i>min_max_lat[1]:
                    min_max_lat[1]=i
                if j<min_max_lon[0]:
                    min_max_lon[0]=j
                if j>min_max_lon[1]:
                    min_max_lon[1]=j
    return [min_max_lat, min_max_lon]     


from numpy.linalg import norm
def closer(x,X):
    #Renvoie l'élément dans X le plus proche de x
    sol=0
    dist=norm(x-X[0])
    for i in range(1,len(X)):
        d=norm(x-X[i])
        if d<dist: dist,sol= d,i
    return sol

def find_closer(grid_lat,grid_lon,x):
    #Renvoie les coordonnées des points de grille qui sont les plus proche du point x=[lat,lon]
    sol=(0,0)
    dist=norm(x-np.array([grid_lat[0,0],grid_lon[0,0]]))
    for i in range(grid_lon.shape[0]):
        for j in range(grid_lon.shape[1]):
            d=norm(x-np.array([grid_lat[i,j],grid_lon[i,j]]))
            if d<dist: dist,sol= d,(i,j)
    return sol

def days(an1,an2):
    #Renvoie le nombre de jours de l'année an1 à an2 comprise
    tot_days=0
    for y in range(an1,an2+1):
        tot_days+=365
        if bissextile(y):
            tot_days+=1
    return tot_days

