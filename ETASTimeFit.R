
library("PtProcess")


  data<-gsub(" ", "", paste(dataname,".csv"))
  XX<-read.csv(data,head=T)
  XX<-XX[XX$mag>=M0,]
 
  xEtas=c()
  XX$time<-gsub("[Z]","",XX$time)
  XX$time<-gsub("[T]"," ",XX$time)
  XX$Y <- as.numeric(substr(XX$time,1,4))
  XX$M <- as.numeric(substr(XX$time,6,7))
  XX$D <- as.numeric(substr(XX$time,9,10))
  XX$Hr <- as.numeric(substr(XX$time,12,13))
  XX$Min <- as.numeric(substr(XX$time,15,16))
  XX$Sec <- as.numeric(substr(XX$time,18,19))
  xEtas$time<-as.numeric(difftime(XX$time,XX$time[1],units="days"))
  xEtas$magnitude<-(XX$mag-min(XX$mag))##Relative magbitude
  xEtas<-as.data.frame(xEtas)
  xEtas<-xEtas[A:B,]
  
  
  #For KK uncomment lines 25-41 and comment 5-21
  # data<-gsub(" ", "", paste(dataname,".csv"))
  # XX<-read.csv(data,head=T)
  # XX<-XX[XX$magnitude>=M0,]
  # 
  # xEtas=c()
  # XX$origintime<-gsub("[Z]","",XX$origintime)
  # XX$origintime<-gsub("[T]"," ",XX$origintime)
  # XX$Y <- as.numeric(substr(XX$origintime,1,4))
  # XX$M <- as.numeric(substr(XX$origintime,6,7))
  # XX$D <- as.numeric(substr(XX$origintime,9,10))
  # XX$Hr <- as.numeric(substr(XX$origintime,12,13))
  # XX$Min <- as.numeric(substr(XX$origintime,15,16))
  # XX$Sec <- as.numeric(substr(XX$origintime,18,19))
  # xEtas$time<-as.numeric(difftime(XX$origintime,XX$origintime[1],units="days"))
  # xEtas$magnitude<-(XX$magnitude-min(XX$magnitude))##Relative magbitude
  # xEtas<-as.data.frame(xEtas)
  # xEtas<-xEtas[A:B,]
  




#----------------------------------------------------------
#    Define Full Model Object

library("PtProcess")

params <- runif(5,0,1)
params[2]=5*params[2]
initial <- log(params[1:5])



#----------------------------------------------------------
#    Define Full Model Object

##The following is essentially the same as what was in the paper
dmagn_mark <- function(x, data, params){
  #  Gamma distribution
  #  exponential density when params[7]=0
  if (params[7]>0){
    lambda <- etas_gif(data, x[,"time"], params=params[1:5])
    y <- dgamma(x[,"magnitude"], shape=1+sqrt(lambda)*params[7],
                rate=params[6], log=TRUE)
  } else y <- dexp(x[,"magnitude"], rate=params[6], log=TRUE)
  return(y)
}

rmagn_mark <- function(ti, data, params){
  #  Gamma distribution
  #  exponential density when params[7]=0
  if (params[7]>0){
    lambda <- etas_gif(data, ti, params=params[1:5])
    y <- rgamma(1, shape=1+sqrt(lambda)*params[7],
                rate=params[6])
  } else y <- rexp(1, rate=params[6])
  return(list(magnitude=y))
}

TT <- c(0, xEtas$time[length(xEtas$time)]) #Where to integrate over

#params <- c(10, 5, 2, 20, 11, 1/mean(dataX$magnitude), 0)
x <- mpp(data=xEtas, gif=etas_gif,
         mark=NULL,
         params=params, TT=TT,
         gmap=expression(params[1:5]),
         mmap=expression(params[1:5]))


#----------------------------------------------------------
#    Fit Model with Exponential Magnitude Distribution


expmap <- function(y, p){
  #   for exponential distribution
  y$params[1:5] <- exp(p)
  return(y)
}
params <- runif(5,0,1)
params[2]=5*params[2]
initial <- log(params[1:5])
z <- optim(initial, neglogLik, object=x, pmap=expmap,
           control=list(trace=1, maxit=100))

initial <- z$par
z <- nlm(neglogLik, initial, object=x, pmap=expmap,
         print.level=2, iterlim=500, typsize=initial)

#    write estimates to a new model object x0
x0 <- expmap(x, z$estimate)
Min<-c(z$minimum)
for(i in 1:9){
  params <- runif(5,0,1)
  params[2]=5*params[2]
  initial <- log(params[1:5])
  z2 <- optim(initial, neglogLik, object=x, pmap=expmap,
              control=list(trace=1, maxit=100))
  
  initial <- z2$par
  z2 <- nlm(neglogLik, initial, object=x, pmap=expmap,
            print.level=2, iterlim=500, typsize=initial)
  Min<-c(Min,z2$minimum)
  if(z2$minimum<z$minimum){
    #    write estimates to a new model object x0
    z<-z2
    
  }
}
print(logLik(x0))


file=gsub(" ", "", paste("A",A,"B",B,data,"TimeETAS.RData"))


save.image(file)


