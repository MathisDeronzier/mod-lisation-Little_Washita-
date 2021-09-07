### exemple de lecture de fichier nc (netcdf)

# library necessaire
library(ncdf4)

# sur obelix

# Nom du fichier
DIR_nc = "/home/hydrogeol2/fpicour/Data_Mathis/"
file.nc = paste(DIR_nc,"CM6/Historical/evspsblpot_Eday_IPSL-CM6A-LR_land-hist_r1i1p1f1_gr_18500101-20141231.nc", sep="")

## ouverture du fichier
cat("Reading",file.nc,"\n")
nc_id = nc_open(file.nc, write=FALSE) # ouverture du fichier

## recuperation des donnees evt (evspsblpot)
evp = (ncvar_get(nc_id,"evspsblpot"))

## recuperation des donnees lat et lon
LAT = ncvar_get(nc_id,"lat")
LON = ncvar_get(nc_id,"lon")

# fermeture du fichier
nc_close(nc_id)


