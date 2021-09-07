#Our convolution will depend if we have even or odd matrices. In the odd case 
#we'll have to rename our lat and long names

convolution_mat <- function(mat, conv,dim.1,dim.2,even){
  conv.mat<-mat
  width<-floor(dim(conv)[1]/2)
  if (even){
    for (i in width:(dim.1-width+1)){
      for(j in width:(dim.2-width)+1){
        conv.mat[i,j]<- tr(t(mat[(i-width+1):(i+width),(j-width+1):(j+width)])%*%conv)
      }
    }
  }
  else{
    for (i in (1+width):(dim.1-width)){
      for(j in (1+width):(dim.2-width)){
        conv.mat[i,j]<- tr(t(mat[(i-width):(i+width),(j-width):(j+width)])%*%conv)
      }
    }
  }
  return(conv.mat)
} 

convolution<-function(data, conv){
  conv.data<-data
  dim.data<-dim(data)
  even<-(dim(conv)[1]%%2==0)
  width<-floor(dim(conv)[1]/2)
  for ( k in (1:dim.data[3])){
    conv.data[,,k] <- convolution_mat(data[,,k],conv,dim.data[1],dim.data[2]) 
  }
  lat <- apply(X = as.matrix(dimnames(data)[[1]]), 1, FUN = as.numeric)
  lon <- apply(X = as.matrix(dimnames(data)[[2]]), 1, FUN = as.numeric)
  conv.lat<-lat
  conv.lon<-lon
  if (even){ 
    for (i in width:(dim.1-width+1)){
      conv.lat[i]<- conv%*% as.matrix(lat[1,(i-width+1):(i+width)])
    }
    for (i in width:(dim.1-width+1)){
      conv.lon[i]<- t(conv)%*%as.matrix(lat[1,(i-width+1):(i+width)])
    }
  }
  dimnames(conv.data)[[1]]<-lat
  dimnames(conv.data)[[2]]<-lon
  return (conv.data)
}

degradation<-function(data,n){
  conv<-matrix(1/(n*n), nrow=n, ncol=n)
  data.conv<-convolution(data,conv)
  return (data.conv)
}