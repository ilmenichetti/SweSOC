# ICBM accessory functions
# PET function
# author: Lorenzo Menichetti, on specifications from Martin Bolinder
# Year 2018, August

#GAI functions with "GAImax values calculated from aboveground NPP", M. Bolinder
PET<-function(latitude, altitude, humidity, temperature, windspeed, sun, sun.mode, date){

  require(lubridate)

  #convert the date in days of the year (lubridate is needed for this step)
  day<-yday(date)

  #to be sure humidity does not go above 100
  humidity[humidity>100]<-100

  #calculate penman=Eto from cloud etc;
  vind2m =windspeed*4.87/(log(67.8*10-5.42));     #eq47;

  #parameters;
  P=101.3*((293-0.0065*altitude)/293)**5.26;    #eq7;
  lambda=2.45;
  cp=0.001013;
  epsilon=0.622;
  gamma=cp*P/(epsilon*lambda); #eq8;

  #vapour pressure deficit;
  es=0.6108*exp(17.27*temperature/(temperature+237.3)); #eq12;
  ea=humidity/100*es; #eq17;
  vpd=es-ea;

  #calulation of shortwave (or solar) radiation Rs from cloudiness;
  Gsc=0.082;
  pi=3.141592654;
  dr=1+0.033*cos(2*pi/365*day);  #eq23;
  de=0.409*sin(2*pi/365*day-1.39);  #eq24;
  radians=pi/180*latitude;  #eq22;
  ws=acos(-tan(radians)*tan(de));  #eq25;
  Ra= 24*60/pi*Gsc*dr*(ws*sin(radians)*sin(de)+cos(radians)*cos(de)*sin(ws));    #eq21;
  Radevap=Ra/lambda;
  Daylight=24/pi*ws;  #eq34;

  ## where the function differs depending on the sunlight mode
  if(sun.mode=="Rsolar"){cat("\r sun mode: Rsolar - calculating from directly measured Rsolar","\n")
    Rs=sun #provided as input;
  }else if (sun.mode=="sunlight") {  #
    cat("\r sun mode: sunlight - Rsolar is nor provided, extrapolating from sunlight","\n")
    Rs=(0.25+0.5*sun)*Ra*1.05
  }else if (sun.mode=="cloudiness") {
    sunlight=1-sun
    cat("\r sun mode: cloudiness - Rsolar is nor provided, extrapolating from cloudiness","\n")
    Rs=(0.25+0.5*sunlight)*Ra*1.05
  }else {
    cat("\r please provide either Rsolar, sunlight (from 0 to 1) or cloudiness (from 0 to 1)","\n")
    }  #eq35 if solar rad not measured;
  #1.05 correction according to calibration;

  #calulation of net radiation and longwave radiation;
  Rso=(0.75+0.00002*altitude)*Ra;  #eq37;
  albedo=0.23;
  Rns=(1-albedo)*Rs;   #eq38;
  TmeanK=273.15+temperature;
  SBolz=4.903*10^(-9);
  RsRso=pmin(1,Rs/Rso);
  Rnl=SBolz*TmeanK**4*(0.34-0.14*ea**0.5)*(1.35*RsRso-0.35); #eq39;
  Rn=Rns-Rnl;

  delta=(4098*(0.6108*exp(17.27*temperature/(temperature+237.3)))/((temperature+237.3)**2));  #eq13;

  Etotal=0.408*delta*Rn + gamma*900/(temperature+273)*vind2m*vpd; #eq6 page 65;
  Etonam=delta+gamma*(1+0.34*vind2m); #eq6 page 65;
  Et0= Etotal/ Etonam; #eq6 page 65;
  Et0[Et0<0]=0

  if(any(is.na(Et0))){
  NAs<-which(is.na(Et0))
  cat("WARNING- there are some NAs:", NAs, "
      Compensating by Stineman interpolation")
  Et0<-na_interpolation(Et0, option="stine")
  }

  return(data.frame(date, ET0=Et0))

}

