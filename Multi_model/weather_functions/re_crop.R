### To compile the manual
### devtools::build_manual()


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
#' Tóth, B., M. Weynants, A. Nemes, A. Makó, G. Bilas, and G. Tóth. 2015. “New Generation of Hydraulic Pedotransfer Functions for Europe: New Hydraulic Pedotransfer Functions for Europe.” European Journal of Soil Science 66 (1): 226–38. https://doi.org/10.1111/ejss.12192.
#' Kätterer, T., O. Andrén, and P-E. Jansson. 2006. “Pedotransfer Functions for Estimating Plant Available Water and Bulk Density in Swedish Agricultural Soils.” Acta Agriculturae Scandinavica, Section B - Plant Soil Science 56 (4): 263–76. https://doi.org/10.1080/09064710500310170.
#' @export
poros <-
  function(sand, clay, SOC)
  { ... }

poros<-function(sand, clay, SOC){

  if(!missing(clay)){
    porosity<-0.4115+0.0409*SOC-0.6089*clay*sand-0.0031*SOC^2*clay+0.2276*clay^3}else{
    porosity<-0.3843+SOC*0.0448+SOC*sand*-0.0204
  }
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
#' Tóth, B., M. Weynants, A. Nemes, A. Makó, G. Bilas, and G. Tóth. 2015. “New Generation of Hydraulic Pedotransfer Functions for Europe: New Hydraulic Pedotransfer Functions for Europe.” European Journal of Soil Science 66 (1): 226–38. https://doi.org/10.1111/ejss.12192.
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
#' Tóth, B., M. Weynants, A. Nemes, A. Makó, G. Bilas, and G. Tóth. 2015. “New Generation of Hydraulic Pedotransfer Functions for Europe: New Hydraulic Pedotransfer Functions for Europe.” European Journal of Soil Science 66 (1): 226–38. https://doi.org/10.1111/ejss.12192.
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


#re_temperature
#' internal function for determining the dependence of decomposition over soil temperature
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param soilT soil temperature daily values
#'
#' @details in case of NAs in the inputs, the function attempts to fill them with a Stineman interpolation
#'
#' @return a vector with the daily water reduction values
#'
#' @references
#'
#' @export
re_temperature <-
  function(soilT)
  { ... }

#dependence of decompositon over temperature, Katterer
re_temperature<-function(soilT){

  tmin=-3.78;
  re_temp<-c()
  for(i in 1:length(soilT)){
  re_temp[i]=(soilT[i]-tmin)^2/(30-tmin)^2
  }

  if(any(is.na(re_temp))){
    NAs<-which(is.na(re_temp))
    cat("WARNING- there are some NAs:", NAs, "
        Compensating by Stineman interpolation")
    water<-na_interpolation(water, option="stine")
  }
  re_temp[soilT<tmin]<-0
  return(re_temp)
}



#re_water
#' internal function for determining the dependence of decompositono over soil moisture
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#' @param twilt wilting point (0 to 1)
#' @param tfield field capacity (0 to 1)
#' @param water water balance (in mm)
#' @param porosity soil porosity /0 to 1)
#' @param L soil depth (in mm)
#'
#' @return a vector with the daily water reduction values
#'
#' @references
#' Moyano FE, Manzoni S, Chenu C. Responses of soil heterotrophic respiration to moisture availability: An exploration of processes and models. Soil Biol Biochem. Elsevier Ltd; 2013;59: 72–85. doi:10.1016/j.soilbio.2013.01.002
#'
#' @export
re_water <-
  function(twilt, tfield, water, porosity, L)
  { ... }

#dependence of decompositono over soil moisture
re_water<-function(twilt, tfield, water, porosity, L){

  rs=0.5;
  topt = 0.2+1.26*tfield**2.03;
  tth= 0.0965*log(twilt) + 0.3;
  theta=water/L;
  re_wat<-c()
  for(i in 1:length(water)){
    if (theta[i]<tth) {re_wat[i]=0}; end;
    if (tth<=theta[i] & theta[i]<min(topt,tfield)){
      re_wat[i]=(theta[i]-tth)/(min(topt,tfield)-tth)}
    if (min(topt,tfield)<= theta[i] & theta[i] <= max(topt,tfield)) {
      re_wat[i]=1}
    if (theta[i] >max(topt,tfield)) {
      re_wat[i]=1-(1-rs)*(theta[i]-max(topt,tfield))/(porosity-max(topt,tfield))}
    }
  re_wat[re_wat>1]=1
  re_wat[re_wat<0]=0
  return(re_wat)
}



