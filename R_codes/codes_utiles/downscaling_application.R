downscal_appli<-function(loc,prog, D_proj,liste.f.g=0){
  source(paste(loc,"/ds_fun/CDFt.R",sep=""))
  source(paste(loc,"/ds_fun/SSR_CDFt.R",sep=""))
  source(paste(loc,"/ds_fun/CDFt_norm.R",sep=""))
  DATES <- read.csv(paste(loc,"/series/dates_1979_2014.csv", sep=""))[,(2:4)]
  #On fixe les dates de calibration et de projection
  cal <- which (DATES[,1]<1996)
  if(D_proj){
    proj <- which (DATES[,1]>=1996)
  }
  else{
    proj <- cal
  }
  DATES_cal <- DATES[cal,]
  DATES_proj <- DATES[proj,]
  
  # Lecture des données 
  
  serie_prec_ref <- read.csv(paste(loc,"/series/precp_ref_1979_2014.csv",sep=""))[,2]
  serie_etp_ref <- read.csv(paste(loc,"/series/evap_ref_1979_2014.csv",sep=""))[,2]
  
  degraded_prec1 <- read.csv(paste(loc,"/series/precp_deg1_1979_2014.csv",sep=""))[,2]
  degraded_etp1 <- read.csv(paste(loc,"/series/evap_deg1_1979_2014.csv",sep=""))[,2]
  
  degraded_prec2 <- read.csv(paste(loc,"/series/precp_deg2_1979_2014.csv",sep=""))[,2]
  degraded_etp2 <- read.csv(paste(loc,"/series/evap_deg2_1979_2014.csv",sep=""))[,2]
  
  degraded_prec3 <- read.csv(paste(loc,"/series/precp_deg3_1979_2014.csv",sep=""))[,2]
  degraded_etp3 <- read.csv(paste(loc,"/series/evap_deg3_1979_2014.csv",sep=""))[,2]
  
  degraded_prec4 <- read.csv(paste(loc,"/series/precp_deg4_1979_2014.csv",sep=""))[,2]
  degraded_etp4 <- read.csv(paste(loc,"/series/evap_deg4_1979_2014.csv",sep=""))[,2]
  
  IPSL_prec <- 1e5*read.csv(paste(loc,"/series/IPSL_pr_1979_2014.csv",sep=""))[,1]
  IPSL_etp <- 1e5*read.csv(paste(loc,"/series/IPSL_etp_1979_2014.csv",sep=""))[,1]
  
  #On construit nos series de calibration et de projection 
  
  #Series de precipitation 
  
  ref_prec_cal <- serie_prec_ref[cal]
  deg1_prec_cal <- degraded_prec1[cal]
  deg2_prec_cal <- degraded_prec2[cal]
  deg3_prec_cal <- degraded_prec3[cal]
  deg4_prec_cal <- degraded_prec4[cal]
  IPSL_prec_cal <- IPSL_prec[cal]
  
  ref_prec_proj <- serie_prec_ref[proj]
  deg1_prec_proj <- degraded_prec1[proj]
  deg2_prec_proj <- degraded_prec2[proj]
  deg3_prec_proj <- degraded_prec3[proj]
  deg4_prec_proj <- degraded_prec4[proj]
  IPSL_prec_proj <- IPSL_prec[proj]
  
  #Series de projection
  
  ref_etp_cal <- serie_etp_ref[cal]
  deg1_etp_cal <- degraded_etp2[cal]
  deg2_etp_cal <- degraded_etp2[cal]
  deg3_etp_cal <- degraded_etp3[cal]
  deg4_etp_cal <- degraded_etp4[cal]
  IPSL_etp_cal <- IPSL_etp[cal]
  
  ref_etp_proj <- serie_etp_ref[proj]
  deg1_etp_proj <- degraded_etp2[proj]
  deg2_etp_proj <- degraded_etp2[proj]
  deg3_etp_proj <- degraded_etp3[proj]
  deg4_etp_proj <- degraded_etp4[proj]
  IPSL_etp_proj <- IPSL_etp[proj]
  
  # Initialisation des vecteurs qui vont stocker les resultats de CDFt
  
  n<-length(proj) #Longueure de la periode de projection
  
  etp_ds1 = array(NaN, dim=n)
  etp_ds2 = array(NaN, dim=n)
  etp_ds3 = array(NaN, dim=n)
  etp_ds4 = array(NaN, dim=n)
  etp_dsI = array(NaN, dim=n)
  
  pr_ds1 = array(NaN, dim=n)
  pr_ds2 = array(NaN, dim=n)
  pr_ds3 = array(NaN, dim=n)
  pr_ds4 = array(NaN, dim=n)
  pr_dsI = array(NaN, dim=n)
  
  #On definit notre fonction de downscaling
  eval(parse(text=paste("fun <- ",prog)))
  eval(parse(text=paste("fun2 <- ","CDFt")))
  if (prog=="CDFt" || prog =="SSR_CDFt" ){
    #boucle sur les 12 mois (à améliorer)
    for(month in 1:12){
      cat("running month =",month,"\n")
      # selection des jours appartenant au mois month pour la calibration
      Imonthcal <- which(DATES_cal[,2]==month)
      #On definit le nombre de pas de discretisation pour le cdf
      NPAS = round((3*length(Imonthcal)/4))
      Imonthproj <- which(DATES_proj[,2]==month)
      #Application de CDF-t
      #Evapo-transpiration
      C_etp1 = CDFt(ref_etp_cal[Imonthcal], deg1_etp_cal[Imonthcal], deg1_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etp2 = CDFt(ref_etp_cal[Imonthcal], deg2_etp_cal[Imonthcal], deg2_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etp3 = CDFt(ref_etp_cal[Imonthcal], deg3_etp_cal[Imonthcal], deg3_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etp4 = CDFt(ref_etp_cal[Imonthcal], deg4_etp_cal[Imonthcal], deg4_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etpI = CDFt(ref_etp_cal[Imonthcal], IPSL_etp_cal[Imonthcal], IPSL_etp_proj[Imonthproj], npas = NPAS)$DS
      
      #Precipitation
      C_pr1 = fun(ref_prec_cal[Imonthcal], deg1_prec_cal[Imonthcal], deg1_prec_proj[Imonthproj], npas = NPAS)$DS
      C_pr2 = fun(ref_prec_cal[Imonthcal], deg2_prec_cal[Imonthcal], deg2_prec_proj[Imonthproj], npas = NPAS)$DS
      C_pr3 = fun(ref_prec_cal[Imonthcal], deg3_prec_cal[Imonthcal], deg3_prec_proj[Imonthproj], npas = NPAS)$DS
      C_pr4 = fun(ref_prec_cal[Imonthcal], deg4_prec_cal[Imonthcal], deg4_prec_proj[Imonthproj], npas = NPAS)$DS
      C_prI = fun(ref_prec_cal[Imonthcal], IPSL_prec_cal[Imonthcal], IPSL_prec_proj[Imonthproj], npas = NPAS)$DS
      
      # Seuillage a 0 des valeurs negatives de precipitation
      C_pr1[which(C_pr1<0)] = 0
      C_pr2[which(C_pr2<0)] = 0
      C_pr3[which(C_pr3<0)] = 0
      C_pr4[which(C_pr4<0)] = 0
      C_prI[which(C_pr4<0)] = 0
      
      # Remplissage des resultats 
      etp_ds1[Imonthproj] = C_etp1
      etp_ds2[Imonthproj] = C_etp2
      etp_ds3[Imonthproj] = C_etp3
      etp_ds4[Imonthproj] = C_etp4
      etp_dsI[Imonthproj] = C_etpI
      
      pr_ds1[Imonthproj] = C_pr1
      pr_ds2[Imonthproj] = C_pr2
      pr_ds3[Imonthproj] = C_pr3
      pr_ds4[Imonthproj] = C_pr4
      pr_dsI[Imonthproj] = C_prI
      
      # supprimer les objets 
      #rm(C_etp1,C_etp2,C_etp3, C_pr1, C_pr2 , C_pr3, Imonthcal, Imonthproj,)
      #gc()
    }
  }
  
  if (prog=="CDFt_norm"){
    #Il faut créer les fonctions f et g ici.
    #boucle sur les 12 mois (à améliorer)
    for(month in 1:12){
      cat("running month =",month,"\n")
      # selection des jours appartenant au mois month pour la calibration
      Imonthcal <- which(DATES_cal[,2]==month)
      #On definit le nombre de pas de discretisation pour le cdf
      NPAS = round((3*length(Imonthcal)/4))
      Imonthproj <- which(DATES_proj[,2]==month)
      #Application de CDF-t
      #Evapo-transpiration
      C_etp1 = CDFt(ref_etp_cal[Imonthcal], deg1_etp_cal[Imonthcal], deg1_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etp2 = CDFt(ref_etp_cal[Imonthcal], deg2_etp_cal[Imonthcal], deg2_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etp3 = CDFt(ref_etp_cal[Imonthcal], deg3_etp_cal[Imonthcal], deg3_etp_proj[Imonthproj], npas = NPAS)$DS
      C_etp4 = CDFt(ref_etp_cal[Imonthcal], deg4_etp_cal[Imonthcal], deg4_etp_proj[Imonthproj], npas = NPAS)$DS
      
      #Precipitation
      C_pr1 = fun(ref_prec_cal[Imonthcal], deg1_prec_cal[Imonthcal], deg1_prec_proj[Imonthproj], liste.f.g[[1]],liste.f.g[[5]])$DS
      C_pr2 = fun(ref_prec_cal[Imonthcal], deg2_prec_cal[Imonthcal], deg2_prec_proj[Imonthproj], liste.f.g[[2]],liste.f.g[[6]])$DS
      C_pr3 = fun(ref_prec_cal[Imonthcal], deg3_prec_cal[Imonthcal], deg3_prec_proj[Imonthproj], liste.f.g[[3]],liste.f.g[[7]])$DS
      C_pr4 = fun(ref_prec_cal[Imonthcal], deg4_prec_cal[Imonthcal], deg4_prec_proj[Imonthproj], liste.f.g[[4]],liste.f.g[[8]])$DS
      
      # Seuillage a 0 des valeurs negatives de precipitation
      C_pr1[which(C_pr1<0)] = 0
      C_pr2[which(C_pr2<0)] = 0
      C_pr3[which(C_pr3<0)] = 0
      C_pr4[which(C_pr4<0)] = 0
      
      # Remplissage des resultats 
      etp_ds1[Imonthproj] = C_etp1
      etp_ds2[Imonthproj] = C_etp2
      etp_ds3[Imonthproj] = C_etp3
      etp_ds4[Imonthproj] = C_etp4
      
      pr_ds1[Imonthproj] = C_pr1
      pr_ds2[Imonthproj] = C_pr2
      pr_ds3[Imonthproj] = C_pr3
      pr_ds4[Imonthproj] = C_pr4
      
      # supprimer les objets 
      #rm(C_etp1,C_etp2,C_etp3, C_pr1, C_pr2 , C_pr3, Imonthcal, Imonthproj,)
      #gc()
      
    }
  }
  
  #On crée les dataframes
  
  #L'evapotranspiration
  etp_ref <- data.frame(ref_etp_proj)
  etp_ds1 <- data.frame(etp_ds1)
  etp_ds2 <- data.frame(etp_ds2)
  etp_ds3 <- data.frame(etp_ds3)
  etp_ds4 <- data.frame(etp_ds4)
  etp_dsI <- data.frame(etp_dsI)
  etp_IPSL <- data.frame(IPSL_etp_proj)
  
  etp_deg1 <- data.frame(deg1_etp_proj)
  etp_deg2 <- data.frame(deg2_etp_proj)
  etp_deg3 <- data.frame(deg3_etp_proj)
  etp_deg4 <- data.frame(deg4_etp_proj)
  
  #La precipitation
  pr_ref <- data.frame(ref_prec_proj)
  pr_ds1 <- data.frame(pr_ds1) 
  pr_ds2 <- data.frame(pr_ds2)
  pr_ds3 <- data.frame(pr_ds3) 
  pr_ds4 <- data.frame(pr_ds4)
  pr_dsI <- data.frame(pr_dsI)
  pr_IPSL <- data.frame(IPSL_prec_proj)
  
  pr_deg1 <- data.frame(deg1_prec_proj)
  pr_deg2 <- data.frame(deg2_prec_proj)
  pr_deg3 <- data.frame(deg3_prec_proj)
  pr_deg4 <- data.frame(deg4_prec_proj)
  
  #On les enregistre
  write.csv(etp_ref, paste(loc,"/series/etp_ref.csv",sep=""))
  write.csv(etp_ds1, paste(loc,"/series/etp_ds1.csv",sep=""))
  write.csv(etp_ds2, paste(loc,"/series/etp_ds2.csv",sep=""))
  write.csv(etp_ds3, paste(loc,"/series/etp_ds3.csv",sep=""))
  write.csv(etp_ds4, paste(loc,"/series/etp_ds4.csv",sep=""))
  write.csv(etp_dsI, paste(loc,"/series/etp_dsI.csv",sep=""))
  write.csv(etp_IPSL, paste(loc,"/series/etp_IPSL.csv",sep=""))
  write.csv(etp_deg1, paste(loc,"/series/etp_deg1.csv",sep=""))
  write.csv(etp_deg2, paste(loc,"/series/etp_deg2.csv",sep=""))
  write.csv(etp_deg3, paste(loc,"/series/etp_deg3.csv",sep=""))
  write.csv(etp_deg4, paste(loc,"/series/etp_deg4.csv",sep=""))
  
  
  write.csv(pr_ref, paste(loc,"/series/pr_ref.csv",sep=""))
  write.csv(pr_ds1, paste(loc,"/series/pr_ds1.csv",sep=""))
  write.csv(pr_ds2, paste(loc,"/series/pr_ds2.csv",sep=""))
  write.csv(pr_ds3, paste(loc,"/series/pr_ds3.csv",sep=""))
  write.csv(pr_ds4, paste(loc,"/series/pr_ds4.csv",sep=""))
  write.csv(pr_dsI, paste(loc,"/series/pr_dsI.csv",sep=""))
  write.csv(pr_IPSL, paste(loc,"/series/pr_IPSL.csv",sep=""))
  write.csv(pr_deg1, paste(loc,"/series/pr_deg1.csv",sep=""))
  write.csv(pr_deg2, paste(loc,"/series/pr_deg2.csv",sep=""))
  write.csv(pr_deg3, paste(loc,"/series/pr_deg3.csv",sep=""))
  write.csv(pr_deg4, paste(loc,"/series/pr_deg4.csv",sep=""))
}