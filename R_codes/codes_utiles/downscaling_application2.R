downscal_appli_2<-function(loc,prog, D_proj){
  source(paste(loc,"/ds_fun/CDFt.R",sep=""))
  source(paste(loc,"/ds_fun/SSR_CDFt.R",sep=""))
  DATES <- read.csv(paste(loc,"/series/dates_1979_2014.csv", sep=""))[,(2:4)]
  #On fixe les dates de calibration et de projection
  cal <- which (DATES[,1]<1996)
  if(D_proj){
    proj <- which (DATES[,1]>=1996)
  }
  else{
    proj <- which (DATES[,1]<1996)
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
  
  # On construit nos series de calibration et de projection 
  
  #Series de precipitation 
  
  ref_prec_cal <- serie_prec_ref[cal]
  deg1_prec_cal <- degraded_prec1[cal]
  deg2_prec_cal <- degraded_prec2[cal]
  deg3_prec_cal <- degraded_prec3[cal]
  deg4_prec_cal <- degraded_prec4[cal]
  
  ref_prec_proj <- serie_prec_ref[proj]
  deg1_prec_proj <- degraded_prec1[proj]
  deg2_prec_proj <- degraded_prec2[proj]
  deg3_prec_proj <- degraded_prec3[proj]
  deg4_prec_proj <- degraded_prec4[proj]
  
  #Series d'etp
  
  ref_etp_cal <- serie_etp_ref[cal]
  deg1_etp_cal <- degraded_etp2[cal]
  deg2_etp_cal <- degraded_etp2[cal]
  deg3_etp_cal <- degraded_etp3[cal]
  deg4_etp_cal <- degraded_etp4[cal]
  
  ref_etp_proj <- serie_etp_ref[proj]
  deg1_etp_proj <- degraded_etp2[proj]
  deg2_etp_proj <- degraded_etp2[proj]
  deg3_etp_proj <- degraded_etp3[proj]
  deg4_etp_proj <- degraded_etp4[proj]
  
  # Initialisation des vecteurs qui vont stocker les resultats de CDFt
  
  n<-length(proj) #Longueure de la periode de projection
  
  etp_ds1 = array(NaN, dim=n)
  etp_ds2 = array(NaN, dim=n)
  etp_ds3 = array(NaN, dim=n)
  etp_ds4 = array(NaN, dim=n)
  
  pr_ds1 = array(NaN, dim=n)
  pr_ds2 = array(NaN, dim=n)
  pr_ds3 = array(NaN, dim=n)
  pr_ds4 = array(NaN, dim=n)
  
  m=25 #On commence par définir la taille du tableau
  #On definit notre fonction de downscaling
  eval(parse(text=paste("fun <- ",prog)))
  #eval(parse(text=paste("fun2 <- ","CDFt")))
  if (prog=="CDFt" || prog =="SSR_CDFt" ){
    #boucle sur les chaque jour, tableau glissant
    for(days in 1:366){
      day<-n_to_day(days)
      if(days%%30==0){cat("running day =",day,"\n")}
      # selection des jours appartenant au mois month pour la calibration
      D<-function(date){return (dist(date,day))}
      #DATES_dist_cal<-apply(DATES_cal,1,D)# Cette merde ne fonctionne apparrement pas
      DATES_dist_cal<-matrix(0,nrow=length(cal),ncol=1)
      for (i in 1:length(cal)){
        DATES_dist_cal[i]<-D(DATES_cal[i,])
      } 
      dayscal <- which(DATES_dist_cal<m)
      #On definit le nombre de pas de discretisation pour le cdf
      NPAS = round((3*length(dayscal)/4))
      daysproj <- which((DATES_proj[,2]==day[2] & DATES_proj[,3]==day[3]))
      #Application de CDF-t
      #Evapo-transpiration
      C_etp1 = CDFt(ref_etp_cal[dayscal], deg1_etp_cal[dayscal], deg1_etp_proj[daysproj], npas = NPAS)$DS
      C_etp2 = CDFt(ref_etp_cal[dayscal], deg2_etp_cal[dayscal], deg2_etp_proj[daysproj], npas = NPAS)$DS
      C_etp3 = CDFt(ref_etp_cal[dayscal], deg3_etp_cal[dayscal], deg3_etp_proj[daysproj], npas = NPAS)$DS
      C_etp4 = CDFt(ref_etp_cal[dayscal], deg4_etp_cal[dayscal], deg4_etp_proj[daysproj], npas = NPAS)$DS
      
      #Precipitation
      C_pr1 = fun(ref_prec_cal[dayscal], deg1_prec_cal[dayscal], deg1_prec_proj[daysproj], npas = NPAS)$DS
      C_pr2 = fun(ref_prec_cal[dayscal], deg2_prec_cal[dayscal], deg2_prec_proj[daysproj], npas = NPAS)$DS
      C_pr3 = fun(ref_prec_cal[dayscal], deg3_prec_cal[dayscal], deg3_prec_proj[daysproj], npas = NPAS)$DS
      C_pr4 = fun(ref_prec_cal[dayscal], deg4_prec_cal[dayscal], deg4_prec_proj[daysproj], npas = NPAS)$DS
      
      # Seuillage a 0 des valeurs negatives de precipitation
      C_pr1[which(C_pr1<0)] = 0
      C_pr2[which(C_pr2<0)] = 0
      C_pr3[which(C_pr3<0)] = 0
      C_pr4[which(C_pr4<0)] = 0
      
      # Remplissage des resultats 
      etp_ds1[daysproj] = C_etp1
      etp_ds2[daysproj] = C_etp2
      etp_ds3[daysproj] = C_etp3
      etp_ds4[daysproj] = C_etp4
      
      pr_ds1[daysproj] = C_pr1
      pr_ds2[daysproj] = C_pr2
      pr_ds3[daysproj] = C_pr3
      pr_ds4[daysproj] = C_pr4
      
      # supprimer les objets 
      #rm(C_etp1,C_etp2,C_etp3, C_pr1, C_pr2 , C_pr3, dayscal, daysproj,)
      #gc()
      
    }
  }
  
  #On cree les dataframes
  #L'evapotranspiration
  etp_ds1 <- data.frame(etp_ds1)
  etp_ds2 <- data.frame(etp_ds2)
  etp_ds3 <- data.frame(etp_ds3)
  etp_ds4 <- data.frame(etp_ds4)
  
  #La precipitation
  pr_ds1 <- data.frame(pr_ds1) 
  pr_ds2 <- data.frame(pr_ds2) 
  pr_ds3 <- data.frame(pr_ds3) 
  pr_ds4 <- data.frame(pr_ds4) 
  
  #On les enregistre
  write.csv(etp_ds1, paste(loc,"/series/etp_ds1.csv",sep=""))
  write.csv(etp_ds2, paste(loc,"/series/etp_ds2.csv",sep=""))
  write.csv(etp_ds3, paste(loc,"/series/etp_ds3.csv",sep=""))
  write.csv(etp_ds4, paste(loc,"/series/etp_ds4.csv",sep=""))
  
  write.csv(pr_ds1, paste(loc,"/series/pr_ds1.csv",sep=""))
  write.csv(pr_ds2, paste(loc,"/series/pr_ds2.csv",sep=""))
  write.csv(pr_ds3, paste(loc,"/series/pr_ds3.csv",sep=""))
  write.csv(pr_ds4, paste(loc,"/series/pr_ds4.csv",sep=""))
}