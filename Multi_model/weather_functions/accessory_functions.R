

# porosity
#'Internal function for determining the soil porosity from texture.
#'
#' If clay is present the function uses  Toth et al., 2015, otherwise  Kätterer et al., 2006
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param sand sand content\%
#' @param clay clay content \% (optional)
#' @param SOC SOC content \%
#'
#' @return a single numerical value with the soil porosity
#'
#' @references
#' Kätterer, T., O. Andrén, and P-E. Jansson. 2006. “Pedotransfer Functions for Estimating Plant Available Water and Bulk Density in Swedish Agricultural Soils.” Acta Agriculturae Scandinavica, Section B - Plant Soil Science 56 (4): 263–76. https://doi.org/10.1080/09064710500310170.
#'
#' #' @export
poros <-
  function(sand, clay, SOC)
  { ... }

poros<-function(SOC, clay){
    porosity<- 0.3697 + clay + 0.0403 * SOC -0.0098 * clay * SOC
  #K?tterer et al 2006
  return(porosity)
}


#WP
#'Internal function for determining the soil wilting point
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param sand sand content \%
#' @param clay clay content \% (optional)
#' @param SOC SOC content\%
#'
#' @return a single numerical value with the soil wilting point
#'
#' @references
#' Kätterer, T., O. Andrén, and P-E. Jansson. 2006. “Pedotransfer Functions for Estimating Plant Available Water and Bulk Density in Swedish Agricultural Soils.” Acta Agriculturae Scandinavica, Section B - Plant Soil Science 56 (4): 263–76. https://doi.org/10.1080/09064710500310170.
#'
#' @export
WP <-
  function(sand, clay, SOC)
  { ... }

WP<-function(sand, clay, SOC){

  WP<-0.0086+0.4473*clay-0.0157*SOC*clay+0.0123*SOC*sand
  #Kätterer et al 2006
  return(WP)
}


#FC
#' internal function for determining the soil field capacity
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param sand sand \%
#' @param SOC SOC \%
#'
#' @return a single numerical value with the soil field capacity
#'
#' @references
#' Kätterer, T., O. Andrén, and P-E. Jansson. 2006. “Pedotransfer Functions for Estimating Plant Available Water and Bulk Density in Swedish Agricultural Soils.” Acta Agriculturae Scandinavica, Section B - Plant Soil Science 56 (4): 263–76. https://doi.org/10.1080/09064710500310170.
#'
#' @export
FC <-
  function(sand, SOC)
  { ... }

FC<-function(sand, SOC){

  FC<-0.4384-0.3839*sand+0.0796*SOC*sand
  #Kaetterer et al 2006
  return(FC)
}


#add.years
#'
#'Extrapolates future climate
#'
#'This function extrapolate future climate based on the average of the past 10 years (or less if less are presente)
#'averaging day by day
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param dataframe the weather data frame
#' @param new.years the number of years you want to extend the weather data frame
#'
#' @return a table, same structure of \link{template}
#'
#' @examples
#'
#'add.years(dataframe, 5)
#'
#' @export
add.years<-function(dataframe, new.years)
  {... }

#function to extrapolate missing weather dates from average
add.years<-function(dataframe, new.years){

  dataframe$date<-as.Date(as.character(dataframe$date))
  add_years<-new.years
  date<-as.Date(dataframe$date)
  last_years<-tail(unique(year(date)),10)
  date_vector<-seq(ymd(as.Date(tail(date,1)+1)), ymd(as.Date(paste(tail(last_years,1)+add_years,"-01-01", sep=""))), by = "day")

  newmat<-mat.or.vec((length(date_vector)-1),dim(dataframe)[2])
  colnames(newmat)<-colnames(dataframe)
  newmat<-as.data.frame(newmat)
  which_ones<-year(dataframe$date) %in% last_years

  for(i in 1:(length(date_vector)-1)){
    DATE_TODAY<-date_vector[i]
    which_days<-day(dataframe[which_ones,]$date)==day(DATE_TODAY)
    which_months<-month(dataframe[which_ones,]$date)==month(DATE_TODAY)
    which_months_day<-which_days & which_months
    newmat[i,]$date<-as.character(DATE_TODAY)
    newmat[i,]$year<-year(DATE_TODAY)
    newmat[i,]$month<-month(DATE_TODAY)
    newmat[i,]$day<-day(DATE_TODAY)
    newmat[i,]$air_temp_deg_C<-mean(dataframe[which_months_day,]$air_temp_deg_C)
    newmat[i,]$precipitation_mm<-mean(dataframe[which_months_day,]$precipitation_mm)
    newmat[i,]$windspeed<-mean(dataframe[which_months_day,]$windspeed)
    newmat[i,]$humidity<-mean(dataframe[which_months_day,]$humidity)
    newmat[i,]$Rsolar<-mean(dataframe[which_months_day,]$Rsolar)
  }

  newdataframe<-rbind(dataframe,newmat)

  return(newdataframe)
}






#fill.na
#'
#'Fill NAs in the data with linear interpolation
#'
#'This function fills NAs based on the average of the past 10 years (or less if less are presente)
#'averaging day by day
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param dataframe the weather data frame
#'
#' @return a table, same structure of \link{dataframe}
#'
#' @examples
#'
#'add.years(dataframe, 5)
#'
#' @export
fill.na<-function(dataframe)
{... }


#function to extrapolate missing weather dates from average
fill.na<-function(dataframe){

  dataframe$date<-as.Date(as.character(dataframe$date))

  dataframe.na<-dataframe
    dataframe.na$air_temp_deg_C<-na.approx(dataframe$air_temp_deg_C)
    dataframe.na$precipitation_mm<-na.approx(dataframe$precipitation_mm)
    dataframe.na$windspeed<-na.approx(dataframe$windspeed)
    dataframe.na$humidity<-na.approx(dataframe$humidity)
    dataframe.na$Rsolar<-na.approx(dataframe$Rsolar)

    which.na<-which(is.na(dataframe.na[,6:9]))
    if(length(which.na)>0){cat("impossible to interpolate, NAs in ", which.na,". you mighthave NAs in the initial or final values")}


  return(dataframe.na)
}







#is.leapyear
#'
#'Identifies if a year is leap or not
#'
#'Well... pretty simple one, main heading is self-explanatory
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param year a year
#'
#' @examples
#'
#'is.leapyear(1996)
#'
#' @export
is.leapyear<-function(year)
{... }

#load the function to indentify leap years
is.leapyear<-function(year){#http://en.wikipedia.org/wiki/Leap_year
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))}



#soiltemp
#' internal function for determining the soil temperature from the air temperature
#'
#' @name soiltemp
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param L soil depth (mm)
#' @param GAI green area index daily values
#' @param date date vector (daily steps)
#' @param temperature air temperature (°C)
#' @param LAI (optional) LAI. If not present LAI is calculated only according to LAI=0.8*GAI, otherwise LAI is used directly
#'
#' @details
#' The function calculates first the surface temperature. If the temperature is below zero:
#' \deqn{T_{surface_i}=0.20 \cdot T}
#' And if the temperature is above zero
#' \deqn{T_{surface_i}=T_i \cdot (0.95+0.05 \cdot exp(-0.4 \cdot (LAI_i-3))}
#' And then calculates the soil temperature according to:
#' \deqn{T_{soil_{i+1}}=T_{soil_i} + (T_{surface_i} - T_{soil_i}) \cdot 0.24 \cdot e^{(-Z_{depth} \cdot 0.017)} \cdot exp(-0.15 \cdot GAI_i)}
#' And where \deqn{Z_{depth}=\frac{L}{20}}
#' The LAI is calculated as \dfunc{LAI= 0.8 \cdot GAI}
#'
#' @return a vector with the daily soil temperature values
#'
#' @references
#' Kätterer, T., and O. Andrén. 2009. “Predicting Daily Soil Temperature Profiles in Arable Soils in Cold Temperate Regions from Air Temperature and Leaf Area Index.” Acta Agriculturae Scandinavica, Section B - Plant Soil Science 59 (1): 77–86. https://doi.org/10.1080/09064710801920321.
#'
#' @encoding UTF-8
#'
#' @export
soiltemp <-
  function(L, GAI, date, temperatur)
  { ... }
soiltemp<-function(L, GAI, date, temperature, LAI=NULL){

  #This functions comes from K?tterer and Andr?n (2008), where L=thickness of topsoil (mm), Zdepth is defined as the mean depth of the topsoil layer (cm), i., midpoint soil depth & and the relationship using LAI = 0.8 x GAI can be seen in Fig. 1 of this publication
  Zdepth=L/20;


  if(is.null(LAI)){LAI=0.8*GAI}
  Tsurface<-c()
  soilT<-c()
  soilT[1]<-0
  for(i in 1:length(temperature)){
    if (temperature[i]<0) {Tsurface[i]=0.20*temperature[i]}
    if (temperature[i]>=0) {Tsurface[i]=temperature[i]*(0.95+0.05*exp(-0.4*(LAI[i]-3)))}
    soilT[i+1]=soilT[i] + (Tsurface[i] - soilT[i])*0.24*exp(-Zdepth*0.017)*exp(-0.15*GAI[i])
  }

  return(head(soilT,-1))
}




#waterbalance
#'Internal function for the water balance model
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param twilt wilting point (0 to 1)
#' @param tfield field capacity (0 to 1)
#' @param precipitation daily precipitations (mm)
#' @param GAI gren area index daily values
#' @param date date vector
#' @param ET0 Evapotranspiration (calculated based on PET and GAI)
#' @param L soil depth (mm)
#'
#' @return The function returns a data frame with water balance and date (days)
#'
#' @details
#' The formulas come mainly from  Allen et al., 1998 <https://www.fao.org/3/x0490e/x0490e00.htm> and it is used to simulate the soil water balance.
#' The calculation is done through multiple steps, iterated for each timestep:
#'
#' \emph{Step 1: Soil water W is initialized assuming saturation, based on the depth L and volumetric capacity}
#'  \deqn{W[1] = \Theta_f \cdot L }
#' \emph{Step 2: The single crop coefficient Kc is calculated based on GAI}
#'  \deqn{K_c=1.3-0.5 \cdot exp(-0.17 \cdot GAI)}
#' \emph{Step 3: calculation of crop evapotranspiration (ETc) under standard condition}
#'  \deqn{ET_c=ET_0 \cdot K_c}
#' \emph{Step 4:  the intercepted water It is calculated based on crop ET, GAI and precipitation P}
#'  \deqn{It=min(P,ET_c,0.2 \cdot GAI)}
#' \emph{Step 5:  potential evapotraspiration is calculated}
#'  \deqn{ E_{pot}=(ET_c-It)}
#' \emph{Step 6:  Calculation of the percolation. Water (W_b, water bypass) is lost when above field capacity, but allowing saturation for one day}
#'  \deqn{W_b = max(0, W-(\Theta_f \cdot L))}
#' \emph{Step 7:  Soil evaporation reduction coefficient}
#'  \deqn{Kr=(1-(0.95 \cdot tfield-\Theta)/(0.95 \cdot tfield-\alpha \cdot twilt))^2}
#' Subsequent conditions are applied so that Kr cannot be above one, and the values before the minimum Kr are also zero.
#' \emph{Step 8:  Actual evapotraspiration is calculated}
#'  \deqn{E_{act}=E_{pot} \cdot Kr }
#' \emph{Step 9:  The water balance is calculated (stepwise)}
#'  \deqn{ W[i+1]=W[i]+P[i]-E_{act}[i]-It-W_b[i]}
#'
#' @md
#'
#' @examples
#'
#'
#' @export
waterbalance <-
  function(twilt, tfield, precipitation, GAI, date, ET0, L)
  { ... }

waterbalance<-function(twilt, tfield, precipitation, GAI, date, ET0, L, alpha=0.7){

  length_sim<-length(precipitation) # define the lenght of the simulation based on the length of the input file

  water<-c()
  Eact<-c()
  bypass<-c() # this is the vector where to store the percolation

  water[1] = tfield*L #setting initial water content to max

  for (i in 1:length_sim){

    # calculating the single crop coiefficient Kc based on GAI
    kc=1.3-0.5*exp(-0.17*GAI[i]);

    # calculation of crop evapotranspiration (ETc) under standard condition
    ETc=ET0[i]*kc;

    inter=min(precipitation[i],ETc,0.2*GAI[i]) #intercepted water
    Epot=(ETc-inter) #potential evapotraspiration

    #percolation option 1, direct percolation. Water is lost when above field capacity
    #if (water[i] >= tfield*L)  {water[i]=tfield*L} #percolation. If the water is more than soil total field capacity, water is lost
    #percolation option 2, saturation is allowed for one day
    bypass[i] = max(0, water[i]-(tfield*L)) #this is the water that will percolate the following day

    theta=water[i]/L;
    Kr=(1-(0.95*tfield-theta)/(0.95*tfield-alpha*twilt))^2
    if (Kr>1){Kr=1}
    Eact[i]=Epot*Kr #actual evapotranspiration
    water[i+1]=water[i]+precipitation[i]-Eact[i]-inter-bypass[i]
  }

  if(any(water<0)){
    cat("WARNING: some water content values are below zero, forcing them to zero...but have a look at the data just in case")
    water[water<0]=0}

  result<-data.frame(head(water,-1), date, Eact)
  result$date<-as.Date(result$date)
  colnames(result)<-c("water", "date", "Eact")
  return(result)
}

