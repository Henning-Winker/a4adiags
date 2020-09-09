#' Convert length frequency to age composition data
#' 
#' Converts a data.frame of multi-annual length frequency into age composition data
#' 
#' @param Ldat a matrix with two columns, where each row specifies the lower and upper length for a given length bin (the lowest should be -Inf, and highest Inf)
#' @param linf Infinitve length
#' @param k Brody growth coefficient
#' @param t0 theoretical age at zero Length 
#' @param cvL, CV of mean length-at-age 
#' @param ages desired age range
#' @param aW parameter of length-weight in mm
#' @param bW parameter of length-weight relationship
#' @param sd.fix if TRUE a fix sd at length-at-age is assumed corresponding to 0.5*Linf
#' @return AgeComp_at, a matrix of expected age-composition data, where columns are samples in year t, and cells are the count of samples with a given age and year
#' @export

a4aLAK = function(Ldat, linf=100, k=0.2, t0= -0.5, cvL=0.15 ,ages=0:5,sd.fix=FALSE,plot=FALSE,aW=NULL,bW=NULL){
  # Calculate expected length at age
  L_a = linf*(1-exp(-k*(ages-t0)))
  if(L_a[1]<0) L_a[1] = 0.2*L_a[2]  
  # Make key to convert age to length
  nages = length(ages) 
  nyrs = ncol(Ldat[,-1])
  nL = nrow(Ldat)   
  Lbin = Ldat[,1]
  bin = median(Ldat[-1,1]-Ldat[-nrow(Ldat),1])/2
  lak = matrix(NA, ncol=length(ages), nrow(Ldat))
  if(sd.fix==T){sdL = cvL*rep(Linf,nages)*0.5} else {sdL = cvL*L_a}
  for( a in 1:ncol(lak)){
    for( l in 1:nrow(lak)){
      lak[l,a] = pnorm(Ldat[l,1]+bin, mean=L_a[a], sd=sdL[a]) - pnorm(Ldat[l,1]-bin, mean=L_a[a], sd=sdL[a])
    }}
 
  # length to age key
  l2a = t(lak)
  l2a = l2a / outer(rep(1,nrow(l2a)),colSums(l2a))
   
  

  # Calculate 
  AgeComps = data.frame(age=ages)
  for(yr in 1:nyrs){
  AgeComps = cbind(AgeComps,l2a %*% Ldat[,yr+1] )
  }
  
  W_a = NULL
  if(is.null(aW)==F & is.null(bW)==F){
    La05 = linf*(1-exp(-k*(ages+0.5-t0)))
    if(La05[1]<=0) La05[1] = 0.2 * La05[2]
    W_a = aW*(La05)^bW/1000
  }
  names(AgeComps) <- c("age",paste(names(Ldat[,-1])))
  
  if(plot==TRUE){
  a4apar(mfrow=c(2,1),plot.cex=0.9)
  plot(Lbin,lak[,1],ylim=c(0,max(lak)*1.1),xlim=range(Lbin),type="n",ylab="Probability")
  Lall = base::rowSums(Ldat[,-1])
  Lall = Lall/max(Lall)*max(lak)
  polygon(c(Lbin,rev(Lbin)),c(Lall,rep(0,nL)),col=grey(0.5,0.5),border=grey(0.5,0.5))  
  for(a in 1:nages) lines(Lbin,lak[,a],col=1,lwd=1)  
  bp = barplot(base::rowSums(AgeComps[,-1]),xlab="Age",ylab="Density",ylim=c(0,1.1*max(base::rowSums(AgeComps[,-1]))))
  box()
  axis(1,at=bp[,1],labels=c(ages))
 }
  
  Return = list()
  Return$AgeComps = AgeComps
  Return$InputLength = Ldat
  Return$AgeLengthKey = lak
  Return$W_a = W_a
  return(Return)
}

