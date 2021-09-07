# CvM pour Mathis 

# libraries needed for cvm
#library(goftest) On le fait avec une alternative

#install.packages("goftest")
library(goftest)


########## Cette fonction est appelee par la fonction "cvm" et ne doit pas etre appelee par l'utilisateur.

cramvonMis<-function(X,Y){
  n<-length(X)
  m<-length(Y)
  rangsXX<-rangs(X,X)
  rangsXY<-rangs(X,Y)
  rangsYY<-rangs(Y,Y)
  rangsYX<-rangs(Y,X)
  Rx<-0
  Ry<-0
  Rxy<-0
  for (i in 1:n){
    Rx<-Rx+ rangsXX[i]^2
    Ry<-Ry+ rangsXY[i]^2
    Rxy<-Rxy+rangsXX[i]*rangsXY[i]
  }
  for (i in 1:m){
    Rx<-Rx+ rangsYY[i]^2
    Ry<-Ry+ rangsYX[i]^2
    Rxy<-Rxy+rangsYY[i]*rangsYX[i]
  }
  return((m*n)/(m+n)^2*(Rx/n^2+Ry/m^2-2*Rxy/(m*n)))
}

########## Elle est a sourcer
CramerVonMisesTwoSamples = function (S1, S2,prec=FALSE){
  xS1 = sort(S1)
  M = length(xS1)
  xS2 = sort(S2)
  N = length(xS2)
  a = data.frame(val = xS1, rang = seq(M), ens = rep(1, M))
  #if(prec){
  #  sec<-which(a[,1]==0)
  #  n=length(sec)
  #  a[,2]<-a[,2]-n
  #  a[which(a[,1]==0),2]=1} 
  b = data.frame(val = xS2, rang = seq(N), ens = rep(2, N))
  #if(prec){
  #  sec<-which(b[,1]==0)
  #  n=length(sec)
  #  b[,2]<-b[,2]-n
  #  b[sec,2]=1} 
  c = rbind(a, b)
  c = c[order(c$val), ]
  c = data.frame(c, rangTot = seq(M + N))
  dtfM = c[which(c$ens == 1), ]
  dtfN = c[which(c$ens == 2), ]
  somN = sum((dtfN$rang - dtfN$rangTot)^2)
  somM = sum((dtfM$rang - dtfM$rangTot)^2)
  U = N * somN + M * somM
  CvM = ((U/(N * M))/(N + M)) - ((4 * M * N - 1)/(6 * (M + N)))
  return(CvM)
}

rangs<-function(X,Y){
  n<-length(X)
  m<-length(Y)
  ordX<-order(X)
  ordY<-order(Y)
  rang<-1
  rangsX<-rep(0,n)
  for(i in 1:n){
    while (rang<=m & X[ordX[i]]>=Y[ordY[rang]]){
      rang<-rang+1
    } 
    rangsX[ordX[i]]<-rang-1
  }
  return(rangsX)
}
################
cvm = function(TS1, TS2, alg=FALSE){
  # compute 2 things:
  # - value of the CvM "distance" T
  # and
  # - p-value of the CvM test performed between time series T1 and T2
  # pval must be < alpha to reject "HO: equality of distribution of the 2 datasets" with a confidence (1-alpha)*100 % (if alpha= 0.05, confidence = 95%)
  # pval > alpha => cannot reject equality => distributions are the same
  # pval < alpha => reject equality        => ditributions are significantly different
  
  L = length(TS1) # je fais l'hypothese que col1 et col2 ont la meme longueur. Si ce n'est pas la meme longueur mais proche ca ne change rien.
  if (alg){T = cramvonMis(TS1,TS2)}
  else{T = CramerVonMisesTwoSamples(TS1,TS2)}# from CDFt
  pval = pCvM(T, n = L, lower.tail = FALSE) # from goftest
  return(list(pval=pval,T=T))
}

# Ici je genere des donnees aleatoirement. En pratique, tu mettras tes donnees.
#Obs = rnorm(10000)
#Model = rnorm(10000)

# calucl de cvm
#CVM = cvm(Obs, Obs)
#CVM$T # donne la valeur de la "distance" entre pdfs
#CVM$pval # p-value of the CvM test

#alpha = 0.05 # Fixe confidence at 95%

#### conclusion
#if(CVM$pval>=alpha){
#  cat("pval > alpha => cannot reject equality => distributions are the same")
#}
#if(CVM$pval<alpha){
#  cat("pval < alpha => reject equality        => ditributions are significantly different (at the chosen confidence level)")
#}



