library(msm)
library(lubridate)
library(MASS)
rm(list=ls())

################################
########## FONCTIONS ########### 
################################


# This is the main function, generating as series of n splashes for one female.
# The series starts at a random date (at which the first batch of oocytes is ripe).
# Then, splashes are performed at random dates, accounting for volleys within night.
# Then, splashes are performed at random hours, independent for different nights but dependent within night.
# Then, splashes occur in random places, identical within volley.
# The function returns a data frame with one line per splash, indicating time and place.
simul.bull<-function(n){
  #Draw a random starting date from a negative binomial
  start.spawn<-ymd(20170524)+days(rnegbin(1,mu=3.9,theta=2.32))  
  # Splashes are performed at random nights, with a lag drawn in a Poisson distribution
  # Occurrence within a volley is drawn in a Bernoulli distribution
  date.bulls<-start.spawn+cumsum(ifelse(rbinom(n,1,0.2),rpois(n,lambda=4),0))  
  # The hour of the first splash of the night is drawn from the Cassou-Leins Normal distribution
  # The place of the first splash of the night is drawn from a truncated Normal centered around a central place
  hour.bulls<-seconds(round(rnorm(n,82800,3600)))
  central.place<-round(runif(1,1,10))
  place.bulls<-round(rtnorm(n,mean=central.place,sd=1,lower=1,upper=10))
  # The delay to other splashes in a volley is drawn from a Gamma distribution
  # All splashes of a volley occur in the same place
   for(i in 2:length(hour.bulls)){
    if(date.bulls[i]==date.bulls[i-1]){ 
      hour.bulls[i]<-hour.bulls[i-1]+seconds(round(rgamma(1,1.0492316876,0.0014073256)))
      place.bulls[i]<-place.bulls[i-1]
    }
  }
  time.bulls<-date.bulls+hour.bulls
  tz(time.bulls)<-"Europe/Paris"
  std.time.bulls<-ymd(20160501)+hour.bulls
  return(data.frame(place.bulls,date.bulls,hour.bulls,time.bulls,std.time.bulls))
}

# We could add a function for social facilitation
# For dates: if there are more than s splashes in my patch less than d days before my next planned splash, I anticipate my splash
# For hours: if there are more than s splashes in my patch less than m minutes before my next planned splash, I anticipate my splash


# This function uses a dataset with 1 line per splash (with its date, hour and place)
# to plot 3 distributions of the number of splashes 1) per spawning ground, 2) per day, 3) per hour accross days
distr.plot<-function(data){
  windows()
  par(mfrow=c(1,3))
  hist(data$place.bulls,freq=T,main="Nb bulls/frayere",xlab="Nb de frayere")
  hist(data$date.bulls,freq=T,breaks=30, main="Saison sur toutes frayeres",xlab="Date")
  hist(data$std.time.bulls,freq=T,breaks=20, main="Distribution horaire",xlab="Heure")
}

# This function applies a sampling plan to the simulated data 
# and returns the data that would have been observed under this sampling plan
sampling<-function(full.data, sampling.matrix=NULL){
  #S'il n'est pas fourni, on d?finit le plan d'?chantillonnage sous forme de matrice jour * 1/4 d'heure, 
  #la cellule contenant le num?ro du site ?chantillonn? (0 indique qu'on n'?tait nulle part)
  if (is.null(sampling.matrix)){
    x<-as.character(date("2016/05/01")+c(0:60)*days(1))
    y<-c("23:0","23:15","23:30","23:45","0:0","0:15","0:30","0:45","1:0","1:15","1:30","1:45","2:0","2:15","2:30","2:45","3:0","3:15","3:30","3:45","4:0","4:15","4:30","4:45")
    z<-c("site1","site2","site3","site4","site5","site6","site7","site8","site9","site10")
    sampling.matrix=matrix(data=round(runif(length(x)*length(y),min=0,max=10)),byrow=T
                          ,length(x),length(y),dimnames=list(x,y))
  }
  sampling.array<-array(data=rbinom(length(x)*length(y)*length(z),1,0.5)
                        ,dim=c(length(x),length(y),length(z)),dimnames=list(x,y,z))
  
  #transforme la matrice en tableau de donn?es avec une ligne par intervalle de temps (length(x)*length(y) lignes au total) pour appliquer le filtre d'?chantillonnage au full.data
  sampling.frame<-data.frame(site=as.vector(t(sampling.matrix)),date=as.POSIXlt(c(rep(x[1],4),rep(x[-1],each=length(y)),as.character(rep(date(x[length(x)])+days(1),20)))),heure=rep(y,length(x)))
  sampling.frame$start<-as.POSIXlt(paste(sampling.frame$date,sampling.frame$heure),format="%Y-%m-%d %H:%M",tz="Europe/Paris")
  sampling.frame$end<-sampling.frame$start+minutes(15)
  sampling.frame$quarter<-which(y==sampling.frame$heure)
  
  # il faut tester si chaque bull est ?chantillonn? (site correspondant observ? au moment correspondant) -> variable sampled = TRUE si echantillonn?, FALSE sinon
  sampled.mat<-matrix(nrow=nrow(full.data),ncol=nrow(sampling.frame))
  for(i in 1:nrow(sampling.frame))
    sampled.mat[,i]<-(full.data$place.bulls == sampling.frame$site[i]) * (full.data$time.bulls %within% interval(sampling.frame$start[i],sampling.frame$end[i]))

  full.data$sampled<-as.logical(rowSums(sampled.mat))
  observed.data<-full.data[full.data$sampled==T,]
  observed.data$quarter<-paste(hour(floor_date(observed.data$std.time.bulls,unit="15 mins")),":",minute(floor_date(observed.data$std.time.bulls,unit="15 mins")),sep="")
  observed.data$n.quarter<-match(observed.data$quarter,y)
  return(observed.data)
}

# Fonction mod?lisant l'estimation d'effectif par la m?thode classique de Cassou-Leins et Cassou-Leins
# Arguments = donn?es observ?es, m?thode d'extrapolation horaire (nuit par nuit et site par site ou en m?langeant toute la saison et tous les sites)
N_classic<-function(observed.data,extrapolation){
  Cassou_Leins<-c(0.41,0.49,0.79,1.55,2.84,3.66,4.55,6.5,7.74,8.71,9.46,9.64,9.53,7.75,6.5,5.41,3.92,3.78,2.81,2.06,1.08,0.52,0.21,0.09)
  fitdistr(observed.data$n.quarter,densfun="normal") #ajuste une loi normale sur le nombre de bulls par quart d'heure
  
}

#######################
### CODE PRINCIPAL ####
#######################

full.data= data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("place.bulls","date.bulls","hour.bulls","time.bulls","std.time.bulls"))),stringsAsFactors=F)

N=100  #nombre de femelles
nbulls<-rpois(N,12.5)  #nombre de bulls par femelle

for(i in 1:N)
 full.data<-rbind(full.data,simul.bull(nbulls[i]))
distr.plot(full.data)


observed.data<-sampling(full.data)
distr.plot(observed.data)

