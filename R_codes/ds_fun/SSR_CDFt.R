
SSR_CDFt <- function(ObsRp, DataGp, DataGf, th=NaN, npas = round((3*length(ObsRp)/4)), dev = 2){
  
  if(is.na(th)){
    th_O = min(ObsRp[which(ObsRp>0)],na.rm=TRUE)
    th_M = min(DataGp[which(DataGp>0)], DataGf[which(DataGf>0)], na.rm=TRUE)
  }
  else{
    th_O = th
    th_M = th
  }
  
  ### st for stoch simulations : from 0 to Unif [0,th]
  ObsRp_st = ObsRp
  DataGp_st = DataGp
  DataGf_st = DataGf
  WObs = which(ObsRp<=th_O)
  ObsRp_st[WObs] = runif(length(WObs),0,th_O)
  WGp = which(DataGp<=th_M)
  DataGp_st[WGp] = runif(length(WGp),0,th_M)
  WGf = which(DataGf<=th_M)
  DataGf_st[WGf] = runif(length(WGf),0,th_M)
  ### 
  
  
  ### Normalization based on the 90th quantile
  Q90O = quantile(ObsRp_st, probs=0.9, na.rm=TRUE)
  Q90Gp = quantile(DataGp_st, probs=0.9, na.rm=TRUE)
  DataGp2 = DataGp_st * (Q90O/Q90Gp)
  DataGf2 = DataGf_st * (Q90O/Q90Gp)
  
  ##
  
  FRp=ecdf(ObsRp_st)
  FGp=ecdf(DataGp2)
  FGf=ecdf(DataGf2)
  
  a=abs(max(DataGf_st, na.rm=TRUE)-max(DataGp_st, na.rm=TRUE))
  m=0
  M=max(ObsRp_st, DataGp_st, DataGf_st, na.rm=TRUE)+dev*a
  x=seq(m,M,length.out=npas)
  FGF=FGf(x)
  FGP=FGp(x)
  FRP=FRp(x)
  
  FGPm1.FGF=quantile(DataGp2,probs=FGF, na.rm=TRUE)
  
  FRF=FRp(FGPm1.FGF)
  
  ######################################
  # FRf=FRp with shift for x<min(DataGf)
  
  if(min(ObsRp_st, na.rm=TRUE)<min(DataGf2, na.rm=TRUE)){
    i=1
    while(x[i]<=quantile(ObsRp_st,probs=FRF[1],na.rm=TRUE)){
      i=i+1
    }
    
    j=1
    while(x[j]<min(DataGf2, na.rm=TRUE)){
      j=j+1
    }
    
    k=i
    while(j>0 && k>0){
      FRF[j]=FRP[k]
      j=j-1
      k=k-1
    }
    
    ##########    
    
    if(j>0){
      for(k in j:1){
        FRF[k]=0
      }
    }
    
  }
  
  ######################################
  # FRf=FRp with shift for x>max(DataGf)
  
  if(FRF[length(x)]<1){
    i=length(x)
    QQ=quantile(ObsRp_st,probs=FRF[length(x)], na.rm=TRUE)
    while(x[i]>=QQ){
      i=i-1
    }
    i=i+1
    
    j=length(x)-1
    while(j>0 && FRF[j]==FRF[length(x)]){
      j=j-1
    }
    
    if(j==0){
      stop("In CDFt, dev must be higher\n")
    }
    
    dif=min((length(x)-j),(length(x)-i))
    FRF[j:(j+dif)]=FRP[i:(i+dif)]
    k=j+dif
    
    if(k<length(x)){
      FRF[k:(length(x))]=1
    }
    
  }
  
  
  ######################################################################################
  ### Quantile-matching based on the new large-scale CDF and downscaled local-scale CDF.
  
  ############
  NaNs.indices = which(is.na(DataGf2))
  No.NaNs.indices = which(!is.na(DataGf2))
  
  qntl = array(NaN, dim=length(DataGf2))
  qntl[No.NaNs.indices] = FGf(DataGf2[No.NaNs.indices])
  
  xx = array(NaN, dim=length(DataGf2))
  xx = approx(FRF,x,qntl,yleft=x[1],yright=x[length(x)],ties='mean')
  
  ##############################################
  
  
  #################
  #################
  
  EmpGf2_th = (ecdf(DataGf2))(th_O)
  EmpGp2m1.Gf2_th = quantile(DataGp2, probs=EmpGf2_th, na.rm=TRUE)
  FRf_th = EmpFp.Gp2m1.Gf2_th = (ecdf(ObsRp_st))(EmpGp2m1.Gf2_th) 
  
  xx$y[which(qntl<=(FRf_th))] = 0 
  xx$y[which(xx$y<=th_O)] = 0 
  
  #################
  #################
  
  
  
  #######################################################################################
  FGp=ecdf(DataGp)
  FGf=ecdf(DataGf)
  FGP=FGp(x)
  FGF=FGf(x)
  
  return(list(x=x,FRp=FRP,FGp=FGP,FGf=FGF,FRf=FRF,FRf_0=FRf_th,DS=xx$y))
  
}

