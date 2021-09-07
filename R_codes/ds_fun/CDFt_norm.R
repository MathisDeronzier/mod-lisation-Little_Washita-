
#Fonction normalisant l'échantillon X
norm<-function(X){
  return ((X-mean(X))/sqrt(var(X)))
}

transform<-function(X){
  return (X/sqrt(var(X)))
}

indice<-function(x,list){#Trouve le plus petit indice a tel que: x>list(a)
  a<-1
  b=length(list)
  while ((b-a)>1){
    m=floor((a+b)/2)
    mid=list[m]
    if (x<mid){b<-m}
    else{a<-m}
  }
  return(a)
}

attribute<-function(x,G,seq){#Interpolation linéaire du quantile-quantile
  a<-indice(x,G[1,])
  if(G[1,a]==x){
    return(seq[a])
  }
  else{
    t<-(x-G[1,a])/(G[1,a+1]-G[1,a])
    return((1-t)*seq[a]+t*seq[a+1])
  }
}
  

CDFt_norm<-function(X,Y,X1,f,g){
  X.n<-transform(X)
  Y.n<-transform(Y)
  X1.n<-transform(X1)
  range.x<-range(cbind(range(X.n),range(Y.n),range(X1.n)))
  x<-exp(seq(log(range.x[1]-0.1),log(range.x[2]+0.1),length.out = floor(4/3*length(X1))))
  G<-matrix(0,nrow=2, ncol=floor(2*length(X1)))
  G[1,]<-ecdf(X.n)(x)
  G[2,]<-ecdf(Y.n)(x)
  X1.cdf<-ecdf(X1.n)
  fun<-function(a){return(attribute(X1.cdf(a),G, x))}
  Y1.n<-apply(as.matrix(X1.n),1,fun)
  Y1<-norm(Y1.n)*g(var(X1))
  return()
}

#Fonction réalisant un quantile-quantile entre les lois de X et celle de Y à partir
#leur fonction de répartition empirique

#WDs<-c("C:/Users/mathisDeronzier/Mathématiques/Stage LSCE downscaling upscaling/stage/NARR datas/R_codes","/home/users/mderon/R_codes")
#wd<-WDs[1]
#setwd(wd)

#serie_prec_ref <- read.csv(paste(loc,"/series/precp_ref_1979_2014.csv",sep=""))[,2]
#degraded_prec1 <- read.csv(paste(loc,"/series/precp_deg1_1979_2014.csv",sep=""))[,2]

#DATES <- read.csv(paste(loc,"/series/dates_1979_2014.csv", sep=""))[,(2:4)]
#On fixe les dates de calibration et de projection
#cal <- which (DATES[,1]<1996)
#X<-degraded_prec1[cal]
#Y<-serie_prec_ref[cal]
#X1<-X

#X.n<-norm(X)
#Y.n<-norm(Y)
#X1.n<-norm(X1)
#range.x<-range(cbind(range(X.n),range(Y.n),range(X1.n)))
#x<-seq(range.x[1]-0.1,range.x[2]+0.1,length.out = floor(4/3*length(X1)))
#G<-matrix(0,nrow=2, ncol=floor(4/3*length(X1)))
#G[1,]<-ecdf(X.n)(x)
#G[2,]<-ecdf(Y.n)(x)
#plot(G[1,],G[2,], type="l")
#Y1.n<-matrix(0,nrow=length(X1), ncol=1)
#X1.cdf<-ecdf(X1.n)
#fun<-function(a){return(attribute(X1.cdf(a),G, x))}
#Y1.n<-apply(as.matrix(X1.n),1,fun)
#plot(x,ecdf(Y1.n)(x)-ecdf(Y1.n)(x), type="l")
#plot(ecdf(Y1.n), col= "blue")
#lines(ecdf(Y.n), col="red")
#cvm(Y.n,Y1.n)
#Y1<-norm(Y1.n)*f(mean(X1))+sqrt(g(var(X1)))

#G[2,1000]
#Y1
