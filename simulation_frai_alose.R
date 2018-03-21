library(msm)
library(lubridate)
library(MASS)
rm(list=ls())

################################
########## FONCTIONS ########### 
################################

#######################################################################################
# This is the main function, generating splashes for a given number of females (nfemales).
# The season starts on "first.day".
# Each female starts spawning at a random number of days after season starts, drawn in a negative binomial of parameters "mu.days" and "theta.days".
# Each female performs a random number of splashes, drawn in a Poisson of parameter "lambda.nbulls"
# Splashes occur either isolately (with probability p.isolated.bull) or in volleys of several splashes on the same night
# Successive nights of activity are separated by a random number of days, drawn in a Poisson of parameter "rest.days"
# Splashes are performed at random hours. The timing of the first splash of the night follows a Gaussian with average "mean.hour" and standard deviation "sd.hour"
# Subsequent splashes in a volley are separated by a random delay, drawn in a Gamma with shape "volley.shape" and rate "volley.rate"
# Spawning area is fragmented in "nplaces" spawning sites, and each female occupies a central place chosen randomly in a uniform.
# Splashes occur in random places (identical within volley), drawn in a Gaussian centered on the central place, with a standard deviation "dispersal". 
# The function returns a data frame with one line per splash, indicating time and place.
simul.bull<-function(nfemales=10,lambda.nbulls=12.5,first.day=20170524,mu.days=3.9,theta.days=2.32,
                     p.isolated.bull=0.2,rest.days=4,mean.hour=82800, sd.hour=3600, nplaces=10,
                     volley.shape=1.0492316876,volley.rate=0.0014073256,dispersal=1){
  full.data= data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("place.bulls","date.bulls","hour.bulls","time.bulls","std.time.bulls"))),stringsAsFactors=F)
  for(f in 1:nfemales){
    nbulls<-rpois(1,lambda.nbulls)
    #Draw a random starting date from a negative binomial
    start.spawn<-ymd(first.day)+days(rnegbin(1,mu=mu.days,theta=theta.days))  
    # Splashes are performed at random nights, with a lag drawn in a Poisson distribution
    # Occurrence within a volley is drawn in a Bernoulli distribution
    date.bulls<-start.spawn+cumsum(ifelse(rbinom(nbulls,1,p.isolated.bull),rpois(nbulls,lambda=rest.days),0))  
    # The hour of the first splash of the night is drawn from the Cassou-Leins Normal distribution
    # The place of the first splash of the night is drawn from a truncated Normal centered around a central place
    hour.bulls<-seconds(round(rnorm(nbulls,mean.hour,sd.hour)))
    central.place<-round(runif(1,1,nplaces))
    place.bulls<-round(rtnorm(nbulls,mean=central.place,sd=dispersal,lower=1,upper=nplaces))
    # The delay to other splashes in a volley is drawn from a Gamma distribution
    # All splashes of a volley occur in the same place
     for(i in 2:length(hour.bulls)){
      if(date.bulls[i]==date.bulls[i-1]){ 
        hour.bulls[i]<-hour.bulls[i-1]+seconds(round(rgamma(1,volley.shape,volley.rate)))
        place.bulls[i]<-place.bulls[i-1]
      }
    }
    time.bulls<-date.bulls+hour.bulls
    tz(time.bulls)<-"Europe/Paris"
    std.time.bulls<-ymd(20160501)+hour.bulls
    full.data<-rbind(full.data,data.frame(place.bulls,date.bulls,hour.bulls,time.bulls,std.time.bulls))
  }  
  return(full.data)
}
############################################################################################
# We could add a function for social facilitation
# For dates: if there are more than s splashes in my patch less than d days before my next planned splash, I anticipate my splash
# For hours: if there are more than s splashes in my patch less than m minutes before my next planned splash, I anticipate my splash

######################################################################################################
# This function uses a dataset with 1 line per splash (with its date, hour and place)
# to plot 3 distributions of the number of splashes 1) per spawning ground, 2) per day, 3) per hour accross days
distr.plot<-function(data){
  windows()
  par(mfrow=c(1,3))
  hist(data$place.bulls,freq=T,main="Nb bulls/frayere",xlab="Nb de frayere")
  hist(data$date.bulls,freq=T,breaks=30, main="Saison sur toutes frayeres",xlab="Date")
  hist(data$std.time.bulls,freq=T,breaks=20, main="Distribution horaire",xlab="Heure")
}
#######################################################################################################
# This function takes a full record of all splashes occurring in a river (all.splashes) and a sampling plan (ech)
# and returns the observations yielded
# all.splashes is a data frame with one row per splash, with its place and time
# ech is a data frame with 3 variables: site number, starting time and ending time of the session
sampling<-function(all.splashes,ech){
  if(missing(ech)){
    #build the sampling plan using start.date, stop.date, n.sites, sampling.intensity, sampling.skew
  }
  sampled.data<-subset(all.splashes,subset=as.logical(colSums(sapply(all.splashes$time.bulls, '%within%', interval(dmy_hms(ech$start,tz="Europe/Paris"),dmy_hms(ech$end,tz="Europe/Paris")),simplify=T)*sapply(all.splashes$place.bulls,'==',ech$site))))
  return(sampled.data)
}

######################################################################################
# Fonction mod?lisant l'estimation d'effectif par la m?thode classique de Cassou-Leins et Cassou-Leins
# Arguments = donn?es observ?es, m?thode d'extrapolation horaire (nuit par nuit et site par site ou en m?langeant toute la saison et tous les sites)
#N_classic<-function(observed.data,extrapolation){
#  Cassou_Leins<-c(0.41,0.49,0.79,1.55,2.84,3.66,4.55,6.5,7.74,8.71,9.46,9.64,9.53,7.75,6.5,5.41,3.92,3.78,2.81,2.06,1.08,0.52,0.21,0.09)
#  fitdistr(observed.data$n.quarter,densfun="normal") #ajuste une loi normale sur le nombre de bulls par quart d'heure
#}

#######################
### MAIN CODE ####
#######################
sampling.plan<-read.csv("G:/Mes Documents/St_Pee/Alose/Simulations/fictional sampling.csv",header=T,sep=";")
actual.nfemales=NULL
actual.nbulls=NULL
observed.nbulls=NULL
for(r in 1:100){
 actual.nfemales[r]<-floor(rnorm(1,100,20))
 full.data<-simul.bull(nfemales=actual.nfemales[r],nplaces = 2)
 actual.nbulls[r]<-nrow(full.data)
 observed.data<-sampling(full.data,sampling.plan)
 observed.nbulls[r]<-nrow(observed.data)
}
windows()
plot(as.data.frame(cbind(actual.nfemales,actual.nbulls,observed.nbulls)))


