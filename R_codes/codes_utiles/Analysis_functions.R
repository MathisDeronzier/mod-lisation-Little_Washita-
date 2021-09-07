#Fonction renvoyant la liste des temps de sécheresse d'affilé
freq_without_rain <- function(pr_serie){
  l<-1
  k<-0
  n<-length(pr_serie)
  drought<-FALSE
  n_drought<-matrix(0,n)
  for (i in 1:n){
    if (pr_serie[i]==0){
      if(drought){
        k <- k+1
      }
      else{
        drought <- TRUE
        k <- 1
      }
    }
    else{
      if(drought){
        n_drought[l] <- k
        l <- l+1
        drought<-FALSE
      }
    }
  }
  return(n_drought[1:(l-1)])
}

#Fonction affichant les débit.
#Fonction regardant les queues des quantiles 
#LA classique rmse

rmse<-function(predicted, real){
  return(sqrt(mean((predicted-real)^2)))
}

### On crée ici la librairie permettant de faire CDFt
### sur un tableau glissant


bis<-function(year){
  return ((year%%400==0) || ((year%%4==0) && (year%%100 !=0)))
}

nbj_y<-function(year){
  if (bis(year)){
    return(c(31,29,31,30,31,30,31,31,30,31,30,31))
  }
  else{return(c(31,28,31,30,31,30,31,31,30,31,30,31))}
}

n_to_day<-function(n){
  nb_j<-c(31,29,31,30,31,30,31,31,30,31,30,31)
  s<-0
  i<-1
  while (s+nb_j[i]<n){
    s<-s+nb_j[i]
    i<-i+1
  }
  return (c(0,i,n-s))
}

day_to_n<-function(date){
  nb_j<-c(31,29,31,30,31,30,31,31,30,31,30,31)
  s<-0
  for (i in 1:(date[2]-1)){
    s<-s+nb_j[i]
  }
  s<-s+date[3]
  return(s)
}

#Fonction renvoyant le nombre de jours entre deux dates
dist<-function(date1,date2){
  date1<-as.matrix(date1)
  nb_j<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  d1<-date1[3]
  d2<-date2[3]
  m1<-date1[2]
  m2<-date2[2]
  if(m1==m2){
    return(abs(d1-d2))
  }
  if (m1>m2){
    m3<-m2
    d3<-d2
    m2<-m1
    d2<-d1
    m1<-m3
    d1<-d3
  }
  s<-0
  if (m2-m1>6){
    for (i in 1:(12+m1-m2)){
      s=s+nb_j[(m2+i-2)%%12+1]
    }
    return(s+d1-d2)
  }
  else{
    for (i in 1:(m2-m1)){
      s<-s+nb_j[m1+i-1]
    }
    return(s+d2-d1)
  }
}
###################### Fonctions utiles pour CDFt ##############################
#Cette fonction donne les rangs des X_i quand on l'applique pour (X,X) donc aussi
#utile pour Cramer Von Mises

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

rangs2<-function(serie){
  ord<-order(serie)
  n<-length(serie)
  minloc<-serie[ord[1]]
  rangs_serie<-rep(0,n)
  rang<-1
  for (i in 1:n){
    if(minloc<serie[ord[i]]){
      rang<-rang+1
      minloc<-serie[ord[i]]
      rangs_serie[ord[i]]<-rang
      cat("i=",i,"  rang=", rang, " serie[ord[",i,"]]=", serie[ord[i]], "\n", sep="")
    }
    else{
      rangs_serie[ord[i]]<-rang
      cat("i=",i,"  rang=", rang, " serie[ord[",i,"]]=", serie[ord[i]],  "\n", sep="")
    }
  }
  return(rangs_serie)
}

#Pour prédire la moyenne en fonction de celle obtenue

sous_part<-function(series){
  p<-1/4
  n<-length(pr_series[,1])
  B<-rbinom(n,1,p)
  return(series[which(B==1),])
}

norm<-function(X){
  return((X-mean(X))/var(X))
}
