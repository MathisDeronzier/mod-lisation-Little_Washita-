cdf_pr_series<-ecdf(etp_series[,1])
x<-seq(0,1,length.out = 100)
plot(x,cdf_pr_series(x),col="red",type="l")
cdf_pr_series<-ecdf(etp_series[,2])
lines(x,cdf_pr_series(x), col="blue")
  seq(0,1,length.out = 100)
cdf_pr_series(seq(0,3,length.out = 100))


colors<-rainbow(4)
plot(x,CDF_pr_series$ref(x),xlab="x = Qtte precip", ylab="y = F_X", type="l",col="blue", main= "Fonctions de distributions")
for (i in 2:2){
  eval(parse(text=paste("y <- CDF_pr_series[[",i,"]]","(x)",sep="")))
  lines(x,y, col = colors[i-1])
}
colors

summary(pr_series)


plot(x,CDF_etp_series$ref(x),xlab="x = Qtte precip", ylab="y = F_X", type="l",col="blue", main= "Fonctions de distributions")
for (i in 2:5){
  eval(parse(text=paste("y <- CDF_etp_series[[",i,"]]","(x)",sep="")))
  lines(x,y, col = colors[i-1])
  }


plot(x,CDF_etp_series$ref(x),xlab="x = Qtte precip", ylab="y = F_X", type="l",col="blue", main= "Fonctions de distributions")
for (i in 2:5){
  eval(parse(text=paste("y <- CDF_etp_series[[",i,"]]","(x)",sep="")))
  lines(x,y, col = colors[i-1])
}

layout((1:2),1,2)
deg<-c(1,9,25,49,81)
plot(deg,pr_CVM_stat, main= "Distance Cramer von Mise séries prec deg", x_lab="deg" , y_lab="CvM", type="l")
plot(deg,etp_CVM_stat, main= "Distance Cramer von Mise séries etp deg", x_lab="deg" , y_lab="CvM", type="l")


ord<-order(Mu[,1])
plot(Mu[ord,1],col="red",type="l")
lines(Mu[ord,2],col="blue")


