
#Programme permettant de créer les fichier dowscalés


# source de la fonction CDFt. Chemin a modifier.
localisations= c(".","/home/hydrogeol2/mderon/NARR_1979_2014")
loc = localisations[1]
source(paste(loc,"/ds_fun/SSR_CDFt.R",sep=""))
source(paste(loc,"/ds_fun/CDFt.R",sep=""))

DATES <- read.csv(paste(loc,"/series/dates_1979_2014.csv", sep=""))[,(2:4)]
#On fixe les dates de calibration et de projection


cal <- which (DATES[,1]<1996)
proj <- which (DATES[,1]>=1996)
DATES_cal <- DATES[cal,]
DATES_proj <- DATES[proj,]

dim(DATES_cal)
dim(DATES_proj)

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


#Serie de reference (non degradees)

ref_prec_cal <- serie_prec_ref[cal]
ref_prec_proj <- serie_prec_ref[proj]

ref_etp_cal <- serie_etp_ref[cal]
ref_etp_proj <- serie_etp_ref[proj]

#1ere degradation (x9)
deg1_prec_cal <- degraded_prec1[cal]
deg1_prec_proj <- degraded_prec1[proj]

deg1_etp_cal <- degraded_etp2[cal]
deg1_etp_proj <- degraded_etp2[proj]

#2eme degradation (x25)
deg2_prec_cal <- degraded_prec2[cal]
deg2_prec_proj <- degraded_prec2[proj]

deg2_etp_cal <- degraded_etp2[cal]
deg2_etp_proj <- degraded_etp2[proj]

#3eme degradation (x49)
deg3_prec_cal <- degraded_prec3[cal]
deg3_prec_proj <- degraded_prec3[proj]

deg3_etp_cal <- degraded_etp3[cal]
deg3_etp_proj <- degraded_etp3[proj]

#4eme degradation (x81)
deg4_prec_cal <- degraded_prec4[cal]
deg4_prec_proj <- degraded_prec4[proj]

deg4_etp_cal <- degraded_etp4[cal]
deg4_etp_proj <- degraded_etp4[proj]

# Initialisation des vecteurs qui vont stocker les resultats de CDFt

n<-length(cal) #Longueure de la periode de projection ici meme que cal

etp_ds1 = array(NaN, dim=n)
etp_ds2 = array(NaN, dim=n)
etp_ds3 = array(NaN, dim=n)
etp_ds4 = array(NaN, dim=n)

pr_ds1 = array(NaN, dim=n)
pr_ds2 = array(NaN, dim=n)
pr_ds3 = array(NaN, dim=n)
pr_ds4 = array(NaN, dim=n)

# boucle sur les 12 mois
for(month in 1:12){
  cat("running month =",month,"\n")
  # selection des jours appartenant au mois month pour la calibration
  Imonthcal = which(DATES_cal[,2]==month)
  #On definit le nombre de pas de discretisation pour le cdf
  NPAS = round((3*length(Imonthcal)/4))
  Imonthproj = which(DATES_proj[,2]==month)
  
  #Application de CDF-t
  #Evapo-transpiration
  C_etp1 = CDFt(ref_etp_cal[Imonthcal], deg1_etp_cal[Imonthcal], deg1_etp_cal[Imonthcal], npas = NPAS)$DS
  C_etp2 = CDFt(ref_etp_cal[Imonthcal], deg2_etp_cal[Imonthcal], deg2_etp_cal[Imonthcal], npas = NPAS)$DS
  C_etp3 = CDFt(ref_etp_cal[Imonthcal], deg3_etp_cal[Imonthcal], deg3_etp_cal[Imonthcal], npas = NPAS)$DS
  C_etp4 = CDFt(ref_etp_cal[Imonthcal], deg4_etp_cal[Imonthcal], deg4_etp_cal[Imonthcal], npas = NPAS)$DS
  
  #Precipitation
  C_pr1 = CDFt_SSR(ref_prec_cal[Imonthcal], deg1_prec_cal[Imonthcal], deg1_prec_cal[Imonthcal], npas = NPAS)$DS
  C_pr2 = CDFt_SSR(ref_prec_cal[Imonthcal], deg2_prec_cal[Imonthcal], deg2_prec_cal[Imonthcal], npas = NPAS)$DS
  C_pr3 = CDFt_SSR(ref_prec_cal[Imonthcal], deg3_prec_cal[Imonthcal], deg3_prec_cal[Imonthcal], npas = NPAS)$DS
  C_pr4 = CDFt_SSR(ref_prec_cal[Imonthcal], deg4_prec_cal[Imonthcal], deg4_prec_cal[Imonthcal], npas = NPAS)$DS
  
  # Seuillage a 0 des valeurs negatives de prÃ©cipitation
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
  #rm(C_etp1,C_etp2,C_etp3, C_pr1, C_pr2 , C_pr3,
  #   Imonthcal, Imonthproj,)
  #gc()
  
}

# verif de base
summary(pr_ds1)
summary(pr_ds2)
summary(pr_ds3)
summary(pr_ds4)

summary(etp_ds1)
summary(etp_ds2)
summary(etp_ds3)
summary(etp_ds4)

summary(serie_prec_ref)
boxplot(serie_prec_ref, pr_ds1)
plot(ecdf(pr_ds1))
lines(ecdf(serie_prec_ref), col=2)

range(serie_etp_ref)
range(etp_ds3)
boxplot(serie_etp_ref, etp_ds1)
plot(ecdf(etp_ds1))
lines(ecdf(serie_etp_ref), col=2)

#On cree les dataframes
#L'evapotranspiration
etp_ds1 <- data.frame(pr_ds1)
etp_ds2 <- data.frame(pr_ds2)
etp_ds3 <- data.frame(pr_ds3)
etp_ds4 <- data.frame(pr_ds4)

#La precipitation
pr_ds1 <- data.frame(pr_ds1) 
pr_ds2 <- data.frame(pr_ds2) 
pr_ds3 <- data.frame(pr_ds3) 
pr_ds4 <- data.frame(pr_ds4) 

#On les enregistre
sum(freq_without_rain(ref_prec_cal))
sum(freq_without_rain(pr_ds4[,1]))
sum(freq_without_rain(deg1_prec_cal))  

write.csv(etp_ds1, paste(loc,"/series/etp_ds1.csv",sep=""))
write.csv(etp_ds2, paste(loc,"/series/etp_ds2.csv",sep=""))
write.csv(etp_ds3, paste(loc,"/series/etp_ds3.csv",sep=""))
write.csv(etp_ds4, paste(loc,"/series/etp_ds4.csv",sep=""))

write.csv(pr_ds1, paste(loc,"/series/pr_ds1.csv",sep=""))
write.csv(pr_ds2, paste(loc,"/series/pr_ds2.csv",sep=""))
write.csv(pr_ds3, paste(loc,"/series/pr_ds3.csv",sep=""))
write.csv(pr_ds4, paste(loc,"/series/pr_ds4.csv",sep=""))

