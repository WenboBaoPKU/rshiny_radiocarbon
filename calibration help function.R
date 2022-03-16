##### load the prior data ######################################
intcal<-ccurve(1)
ca.curve<-as.data.frame(intcal)
# There are some bugs in the following function ,but it can work now.
#################################################### 
Mu_To_Theta_Min <-function(rdata,error,cc=ca.curve)  {
  
  a<-which(ca.curve[,2]<=rdata-7*error)
  ncp<-1
  crossp<-1
  if(length(a)==0){
    ncp=0
    print('there is no cross point..')
    corssp<-0
  }else{
    for(i in 1:(length(a)-1)){
      if(a[i+1]-a[i]>1){
        print(i+min(a)-1)
        crossp<-i+min(a)-1
        ncp=2
        print('there are more than one cross points..')
        break
      }
    }
    if(ncp==1){
      print('there is only one cross point..')
      print(max(a))
      crossp<-max(a)
    }
  }
  return(ca.curve[crossp,1])
}
#---------------------------------------------------
Mu_To_Theta_Max <- function(rdata,error,cc=ca.curve){
  
  a<-which(ca.curve[,2]<=rdata+7*error)
  ncp<-1
  crossp<-1
  if(length(a)==0){
    ncp=0
    print('there is no cross point.')
    corssp<-9501
  }else{
    for(i in 1:(length(a)-1)){
      if(a[i+1]-a[i]>1){
        print(i+min(a)-1)
        crossp<-i+min(a)-1
        ncp=2
        print('there are more than one cross points.')
        break
      }
    }
    if(ncp==1){
      print('there is only one cross point.')
      print(max(a))
      crossp<-max(a)
    }
  }
  
  return(ca.curve[crossp,1])
 
}
####################################################
# Theta_To_Mu<-function(theta,cc=ca.curve){
# 
#   mu<-as.numeric()
# 
#   for(i in 1:length(theta)){
#     #
#     id <- sum(cc[,1] <= theta[i])
#     b <- (cc[id+1,2]-cc[id,2])/(cc[id+1,1]-cc[id,1])
#     mu[i] <- b*(theta[i]-cc[id+1,1])+cc[id+1,2]
#   }
# 
#   return(mu)
# }
####################################################
Theta.MuER<- function(theta, cc=ca.curve, rule=1) {
  
  mu <- approx(cc[,1], cc[,2], theta, rule=rule)$y
  er <- approx(cc[,1], cc[,3], theta, rule=rule)$y
  return(data.frame(theta,mu, er))
}
####################################################
# Theta_To_Sigma <- function(theta=theta,cc=ca.curve){
#   
#   sigma <- as.numeric()
#   
#   for(i in 1:length(theta)){
#     #
#     id <- sum(cc[,1] <= theta[i])
#     b <- (cc[id+1,3]-cc[id,3])/(cc[id+1,1]-cc[id,1])
#     sigma[i] <- b*(theta[i]-cc[id+1,1])+cc[id+1,3]
#   }
#   return(sigma)
# } 
####################################################
# Make_Post <- function(theta,cc=ca.curve){
#   
#   mu<-Theta_To_Mu(theta=theta,cc=ca.curve)
#   sigma<-Theta_To_Sigma(theta=theta,cc=ca.curve)
#   
#   post<-data.frame(theta,mu,sigma)
#   return(post)
#   
# }
####################################################
Grid_Post <- function(post,rdata,error){
  post$LL<-sapply(1:nrow(post),
                  function(i)dnorm(rdata,mean=post$mu[i],sd=sqrt(error^2+post$er[i]^2),log=TRUE))
  
  post$prob<-post$LL+dunif(post$theta,1,50000,log=TRUE)
  #head(post)
  post$prob<-exp(post$prob)
  post$prob<-post$prob/sum(post$prob)
  
  return(post)
  
}
####################################################
PlotCurve <- function(post,pm,psd,rdata,error){
  ##### draw the posterior density ͼ######
  xrange <- post$theta
  ydens <- post$prob
 
  opar<-par(no.readonly = TRUE)
  par(fig=c(0,1,0,0.8))
  plot(x=xrange,y=ydens,
       type='n',
       col='blue',
       yaxt='n',
       xlim=c(range(xrange)[2]+50,range(xrange)[1]-50),
       ylim=c(0,max(ydens)*1.1),
       xlab='BP Calendar Age (Year)',ylab='Probability')
  
  xx<-c(c(xrange[1],xrange),rev(c(xrange[1],xrange)))
  yy<-c(c(ydens[1],ydens),rep(0,times=length(ydens)+1))
  polygon(x=xx,y=yy,col='darkgrey',border = 'white')

  text(x=range(xrange)[1],
       y=max(ydens),
       family='serif',
       cex=1,
       paste(' Mean :',pm,'\n','Sd :',psd))
  
  ##### draw the calibrate curve #########
  mu <- post$mu
  er <- post$er
  
  
  par(fig=c(0,1,0.50,1),new=TRUE)
  plot(x=xrange,
       y=mu,
       type='n',
       #main='IntCal20',
       
       xlim=c(range(xrange)[2]+50,range(xrange)[1]-50),
       #ylim=c(min(mu)-5,max(mu)+10),
       col='green',xlab='',ylab='14C Age',
       xaxt='n'
  )

  xx<-c(c(xrange[1],xrange),rev(c(xrange[1],xrange)))
  yy<-c(c(mu[1],mu+er),rev(c(mu[1],mu-er)))
  polygon(xx,yy,col='#FFB6C1',border = 'white')
  
  text(x=range(xrange)[1],
       y=(max(mu)-min(mu))*0.9+min(mu),
       #family='serif',
       cex=0.6,
       paste(' IntCal20'))
  
  text(x=range(xrange)[1],
       y=(max(mu)-min(mu))*0.7+min(mu),
       family='serif',
       cex=1,
       paste('Radiocarbon Age\n',rdata,'+/-',error))
  abline(h=rdata,col='grey',lty=6)
  abline(h=rdata+error/2,col='pink',lty=6)
  abline(h=rdata-error/2,col='pink',lty=6)
  par(opar)
  
}




