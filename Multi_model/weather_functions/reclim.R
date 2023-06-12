
### command to build .pdf documentation
###devtools::build_manual()

### re_clim main calculations
# author : Lorenzo Menichetti
# ilmenichetti@gmail.com
# lorenzo.menichetti@slu.se

#load the libraries
library("zoo")
library("lubridate")
library("imputeTS")
library("RColorBrewer")
library("anytime")

### Roxygen documentation

#reclim
#'Wrapper for the functions calculating the ICBM climate scaling factors
#'
#'This functions runs the re_clim calculation on a dataset (composed by two different tables, one daily for weather and one annual for aboveground biomass, please refer to the \link{template} for the data structure)
#'The function is a wrapper, performing a few data checks and running functions to calculate several parameters
#'and hopefully runs without the user having to bother too much with intermediate steps.
#'
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}, Martin Bolinder, Olaf Andrén, Thomas Kätterer
#'
#' @param weather data matrix of weather data, must be exactly in the format of the  \link{template} attached as example and contain the following headers:
#'  ("date", "year", "month", "day", "air_temp_deg_C", "precipitation_mm", "windspeed_kmh", "humidity_percent", "Rsolar_lang")
#'
#' @param aboveground data matrix of weather data, must be exactly in the format of the  \link{template} attached as example and contain the following headers:
#'  ("year", "crop_description", "crop_id", "treat", "variance", "seeding", "harvest", "harvest2", "tillage", "minimum_cover", "total_dm_kg_ha", "total_dm_kg_ha2" )
#'  "harvest2" and "total_dm_kg_ha2" are optional and used in case of a double cut for leys
#' @param sun.mode mode of sun data, can be either "Rsolar" (expressed in Langleys) or "cloudiness" (expressed in percent of sunny time per day)
#' @param latitude the latitude, in degrees
#' @param altitude altitude in meters
#' @param depth depth considered in centimeters
#' @param sand sand, in \%. This is needed if porosity, wilting point and field capacity are not specified
#' @param clay clay, in \%. This is needed if porosity, wilting point and field capacity are not specified
#' @param ave_SOC average SOC over the whole period, in \%. This is needed if porosity, wilting point and field capacity are not specified
#' @param porosity soil porosity, as 0 to 1. If speciefied with wilting point and field capacity there's no need for other soil edaphic properties.
#' @param wilting_point wilting point, as mm over the total. If speciefied with porosity and field capacity there's no need for other soil edaphic properties.
#' @param field_capacity field capacity , as mm over the total. If speciefied with wilting point and porosity there's no need for other soil edaphic properties.
#'
#' @return
#'
#' @examples
#'
#'reclim_out<-reclim(weather=weather_testdata,
#'                   aboveground=aboveground_testdata,
#'                   latitude=44,
#'                   altitude=20,
#'                   sand=22,
#'                   clay=36,
#'                   ave_SOC=1.2,
#'                   depth=20,
#'                   sun.mode="Rsolar")
#'
#' @export
reclim <-
  function(weather, aboveground, sun.mode, latitude, altitude, depth, sand=NULL,
           clay=NULL, ave_SOC=NULL, porosity=NULL,
           wilting_point=NULL, field_capacity=NULL)
  { ... }


reclim <- function(weather, #table
                   aboveground, #table
                   sun.mode, #solar calculation mode, "Rsolar" or cloudiness
                   latitude, #degrees
                   altitude, #meters
                   depth,#centimeters
                   sand=NULL,# sand in percent
                   clay=NULL, #clay in percent
                   ave_SOC=NULL,#average SOC in percent over the whole time series (for soil phisical properties calculation)
                   porosity=NULL,
                   wilting_point=NULL,
                   field_capacity=NULL) {

  #get the levels of the treatments, forcing it to factor just in case
  treat_levels<-levels(as.factor(aboveground$treat))

  cat(paste("Hi there! we're going to try our best to calculate your ICBM climatic reduction factors...","\n","please forgive me for any issue, and write me an email at ilmenichetti@gmail.com, we can maybe try to smooth it out together","\n","\n"))

  cat("performing basic data integrity check on weather data","\n")

  ## perform some basic integrity tests on the input data
  # test the dates
    # make sure the date is in the correct format
    weather$date<-anydate((weather$date))
    # check that the first date is a 1st January
    start_on_january<-format(as.Date(weather$date[1]), "%d-%m")=="01-01"
    if(!start_on_january){ stop("either your weather time series does not start on the first of January or is not in the right format (YYYY-mm-dd, separaed by hyphens), please check and correct")}
    # check that the last date is a 31st december
    end_on_december<-format(tail(weather$date,1), "%d-%m")=="31-12"
    if(!end_on_december){ stop("your weather time series does not end on the last of December, correct")}
    #check if there are nas in the dates
    if(any(is.na(weather[,1:4]))){ stop("you have NAs in the dates, correct")}
    #check if there are gaps in the dates
    date_vec<-seq(weather$date[1], tail(weather$date,1), by = "day")
    holes_in_dates<- length(weather$date) == length(date_vec)
    which_missing<-which(!date_vec %in% weather$date)
    which_missing_message<-paste(which_missing, collapse =",")
    if(!holes_in_dates){stop(paste("you have missing days in your weather time series, positions nr",which_missing_message,", correct"))}

  #test the weather time series
    # check eventual NAs in the air temperature
    if(any(is.na(weather$air_temp_deg_C))){stop("some NAs detected in air temperature time series, correct")}
    if(any(weather$air_temp_deg_C>60)){warning("Hey, some of your temperature values are above 60 degrees C...u sure?")}
    # check eventual NAs in the precipitations
    if(any(is.na(weather$precipitation_mm))){stop("some NAs detected in precipitation series, correct")}
    if(any(weather$precipitation_mm>200)){warning("Hey, some of your temperature values are above 200mm (in a day)...u sure?")}
    # fill eventual NAs in the windspeed
    if(any(is.na(weather$windspeed_kmh))){stop("some NAs detected in windspeed series, correct")}
    if(any(weather$windspeed_kmh>160)){warning("Hey, some of your windspeed values are above 160 km/h... u sure?")}
    # check eventual NAs in the humidity
    if(any(is.na(weather$humidity_percent))){stop("some NAs detected in humidity series, correct")}
    if(any(weather$humidity_percent>100)){stop("some of your humidity values are above 100 percent, correct")}

    # check eventual NAs in the solar radiation data
    	if(sun.mode=="Rsolar"){
    	  if(any(is.na(weather$Rsolar_lang))){stop("some NAs detected in solar radiation series, correct")}
    	  if(any(weather$Rsolar_lang>200)){warning("Hey, some of your solar radiation values are above 200 Langley... u sure?")}
    	}	else if (sun.mode == "cloudiness"){
    	  if(any(is.na(weather$cloudiness))){stop("some NAs detected in couldiness series, correct")}
    	  if(any(weather$cloudiness>100)){stop("some of your cloudiness values are above 100 percent, correct")}
    	}	else { stop("you did not specify correctly the mode for solar radiation (Rsolar or cloudiness)")}

  cat("## basic data integrity check on weather data cleared","\n")


  cat("performing basic data integrity check on aboveground data","\n")

  # perform some basic integrity tests on the input data
    #test that there are no gap years
    treatment1<-treat_levels[1]
    years_vec<-seq(from=aboveground$year[1], to=tail(aboveground$year,1), by=1)
    holes_in_years<- length(aboveground[aboveground$treat==treatment1,]$date) == length(years_vec)
    which_missing_AG<-which(years_vec %in% aboveground[aboveground$treat==treatment1,]$date)
    which_missing_AG_message<-paste(which_missing_AG, collapse=",")
    if(holes_in_years){ stop(paste("you have missing years in your aboveground time series, positions nr",which_missing_AG_message,", correct"))}
    #check eventual NAs in the crop_id
    if(any(is.na(aboveground$crop_id))){stop("some NAs detected in crop_id series, correct")}
    #check eventual NAs in the treat
    if(any(is.na(aboveground$treat))){stop("some NAs detected in treat series, correct")}
    #check eventual NAs in the total_dm_kg_ha
    if(any(is.na(aboveground$total_dm_kg_ha))){stop("some NAs detected in total_dm_kg_ha series, correct")}
    #check eventual wrong codes in the crop_id
    GAI_crop_list<-c("spring_small_grains", "spring_oil_seeds","winter_small_grains", "winter_oil_seeds","root_crop", "fodder","fodder_maize", "ley1", "ley", "missing")
    aboveground_error<-!(aboveground$crop_id %in% GAI_crop_list)
    which_wrong_crop<-which(aboveground_error)
    which_wrong_crop_message<-paste(which_wrong_crop, collapse =",")
    if(any(aboveground_error)){stop(paste("it looks like there are some errors in how you identified the crop IDs, check. Errors are in positions", which_wrong_crop_message ))}

  cat("## basic data integrity check on aboveground data cleared","\n")

  cat("performing data integrity check on the consistency of the two dataset","\n")

    # check that the datasets cover the same period
    years_in_weather_only<- !(unique(format(weather$date, "%Y")) %in% years_vec)
    years_in_aboveground_only<- !(years_vec %in% unique(format(weather$date, "%Y")))
    which_years_in_weather_only<-paste(unique(format(weather$date, "%Y"))[which(years_in_weather_only)], collapse=",")
    which_years_in_aboveground_only<-paste(years_vec[which(years_in_aboveground_only)], collapse=",")

    if(any(years_in_weather_only) | any(years_in_aboveground_only)){
      stop(paste("the datasets cover a different period.","\n","Years in weather only:",which_years_in_weather_only,"\n", "Years in aboveground only:", which_years_in_aboveground_only, "\n","\n", "Please reconciliate, the periods need to be the same","\n","(tip: you can decide to extend with means, but you need to make some assumptions there)"))
    }

  cat("## datasets appear to cover the same period","\n")

  cat("checking if basic soil phisical data are available","\n")

  #check if porosity, wilting point are there, if not check if the variables for inferring them are there
  if(is.null(porosity) | is.null(wilting_point) | is.null(field_capacity)){
    if(is.null(sand) | is.null(clay) | is.null(ave_SOC)){
      stop("You must specify EITHER porosity AND wilting_point AND field_capacity OR sand AND clay AND ave_SOC")}
    }

  if(is.null(sand) | is.null(clay) | is.null(ave_SOC)){
    if(is.null(porosity) | is.null(wilting_point) | is.null(field_capacity)){
      stop("You must specify EITHER porosit AND, wilting_point AND field_capacity OR sand AND clay AND ave_SOC")}
    }

  cat("## basic soil phisical data are present","\n")


  ## calculating the climatic parameters, loop for treatments

  #PET
  cat("PET calculation","\n")

  PET_calc<-PET(humidity=weather$humidity,
                windspeed=weather$windspeed,
                temperature=weather$air_temp_deg_C,
                latitude=latitude,
                altitude=altitude,
                sun=weather$Rsolar,
                sun.mode="Rsolar", # can be Rsolar, sunlight or cloudiness, latter two between 0 and 1
                date=as.Date(weather$date))

  #calculate porosity, wilting point and field capacity if missing
  if(!exists("porosity")){
    cat("porosity calculation","\n")
    porosity<-poros(sand=sand/100, clay=clay/100, SOC=ave_SOC)
    }
  if(!exists("wilting_point")){
    cat("wilting point calculation","\n")
    wilting_point<-WP(sand=sand/100, clay=clay/100, SOC=ave_SOC)
    }
  if(!exists("field_capacity")){
    cat("field capacity calculation","\n")
    field_capacity<-FC(sand=sand/100, SOC=ave_SOC)
    }


  #create the tables to store the values for each treatment
  GAI_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  water_bal_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  actevapo_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  soilT_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  re_wat_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  re_temp_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  re_crop_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  re_x1_tab<-mat.or.vec(length(treat_levels), length(weather$date))
  rownames(GAI_tab)<-treat_levels
  rownames(water_bal_tab)<-treat_levels
  rownames(actevapo_tab)<-treat_levels
  rownames(soilT_tab)<-treat_levels
  rownames(re_wat_tab)<-treat_levels
  rownames(re_temp_tab)<-treat_levels
  rownames(re_crop_tab)<-treat_levels
  rownames(re_x1_tab)<-treat_levels


  for(i in 1:length(levels(as.factor(aboveground$treat)))){ #loop for treatments starts

  treatment<-treat_levels[i]
  selected_aboveground<-aboveground[aboveground$treat==treatment,]

  cat(paste("\n","performing calculations for treatment",treatment,"\n"))

  cat(paste("GAI calculation for treatment",treatment,"\n"))
  GAI_calc<-GAI(yield=selected_aboveground$total_dm_kg_ha,
                year=selected_aboveground$year,
                crop=selected_aboveground$crop_id,
                variance=selected_aboveground$variance,
                seeding=selected_aboveground$seeding,
                harvest=selected_aboveground$harvest,
                tillage=selected_aboveground$tillage,
                minimum_cover=selected_aboveground$minimum_cover,
                yield2=selected_aboveground$total_dm_kg_ha2,
                harvest2=selected_aboveground$harvest2)



  #soil temperature
  cat(paste("Soil T calculation for treatment",treatment,"\n"))
  soilT<-soiltemp(temperature=weather$air_temp_deg_C,
                  L=depth*10,
                  GAI=GAI_calc$GAI)
  soilT_tab[i,]<-soilT

  #water balance simulation
  cat(paste("Water balance calculation for treatment",treatment,"\n"))
  water_bal<-waterbalance(twilt=wilting_point,
                          tfield=field_capacity,
                          ET0=PET_calc$ET0,
                          precipitation=weather$precipitation_mm,
                          GAI=GAI_calc$GAI,
                          date=GAI_calc$date,
                          L=depth*10)
  water_bal_tab[i,]<-water_bal$water
  actevapo_tab[i,]<-water_bal$Eact

  #re_wat
  cat(paste("re_wat calculation for treatment",treatment,"\n"))
  re_wat<-re_water(twilt=wilting_point,
                 tfield=field_capacity,
                 porosity=porosity,
                 water=water_bal$water,
                 L=depth*10)
  re_wat_tab[i,]<-re_wat

  #re_temp
  cat(paste("re_temp calculation for treatment",treatment,"\n"))
  re_temp<-re_temperature(soilT)
  re_temp_tab[i,]<-re_temp

  #re_crop
  re_x1=re_wat*re_temp
  #re_crop=re_x1/0.10516
  re_crop=re_x1/0.1056855 #updated on new value, August 2021
  re_crop_tab[i,]<-re_crop
  re_x1_tab[i,]<-re_x1


  } #end of loop for treatemtns


  results_daily<-data.frame(PET_calc, t(GAI_tab), t(water_bal_tab), t(soilT_tab),
                            t(re_wat_tab), t(re_temp_tab), t(re_crop_tab), t(re_x1_tab), t(actevapo_tab))

  #create the name vector for the result table
  shortnamelist<-c()
  shortnamelist2<-c()
  shortnamelist3<-c()
  shortnamelist4<-c()
  shortnamelist5<-c()
  shortnamelist6<-c()
  shortnamelist7<-c()
  shortnamelist8<-c()
  for(j in 1:length(treat_levels)){
    shortnamelist[j]<-paste("GAI_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist2[j]<-paste("water_bal_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist3[j]<-paste("soil_t_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist4[j]<-paste("re_wat_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist5[j]<-paste("re_temp_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist6[j]<-paste("re_crop_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist7[j]<-paste("re_x1_treat",treat_levels[j], sep=".")
  }
  for(j in 1:length(treat_levels)){
    shortnamelist8[j]<-paste("actevapo_treat",treat_levels[j], sep=".")
  }
  colnames(results_daily)<-c("date","PET",shortnamelist,shortnamelist2, shortnamelist3, shortnamelist4,shortnamelist5,shortnamelist6,shortnamelist7, shortnamelist8)


  #results_list<-list(results_daily, PET_calc, (GAI_tab), (water_bal_tab), (soilT_tab), (re_wat_tab), (re_temp_tab), (re_crop_tab), GAI_calc$crop,  porosity, wilting_point, field_capacity)
  results_list<-list(results_daily, PET_calc, (GAI_tab), (water_bal_tab), (soilT_tab), (re_wat_tab), (re_temp_tab), (re_crop_tab), (re_x1_tab), GAI_calc$crop,  porosity, wilting_point, field_capacity, actevapo_tab)
  #names(results_list)<-c("results_daily", "PET", "GAI", "water_bal", "soil_temp", "re_wat", "re_temp", "re_crop", "crop_id",  "porosity", "wilting_point", "field_capacity")
  names(results_list)<-c("results_daily", "PET", "GAI", "water_bal", "soil_temp", "re_wat", "re_temp", "re_crop", "re_x1", "crop_id",  "porosity", "wilting_point", "field_capacity", "actevapo")
  return(results_list)

}


# function to aggregate the calculated values by year
reclim_annual <- function(results_daily){
  DF_yearly<-aggregate(data.frame(results_daily), by=list(year(results_daily$date)), FUN="mean")

  return(DF_yearly)
  }


### collection of diagnostic plots
reclim_diagnostic_plots<-function(results_list, output_to_file=F){

  treatments<-rownames(results_list$GAI)

  crop_palette<-brewer.pal(8, "Dark2")

  if(output_to_file==TRUE){

  for(j in 1:length(treatments)){
    filename1<-paste("GAI_treatment_",j,".png", sep="")
    png(file=filename1, res=300, height=1500, width=2300)
  #plot the effect of the GAI function{
    levels<-levels(results_list$crop_id)
    date_range<-range(results_list$results_daily$date)
    plot(results_list$results_daily$date, results_list$GAI[j,], type="l",  lty=2, xlab="Time", ylab="GAI", main=treatments[j])
    for(i in 1:length(levels)){
      GAI_loop<-results_list$GAI[j,]
      GAI_loop[!results_list$crop_id==levels[i]]<-NA
      lines(results_list$results_daily$date, GAI_loop, col=crop_palette[i])
    }
    legend("topleft", levels, col=crop_palette, lty=1, pch=NA, bty="n")
    dev.off()

    filename2<-paste("meteo_",j,".png", sep="")
    png(file=filename2, res=300, height=3000, width=2300)
    par(mfrow=c(2,1))
    plot(results_list$results_daily$date, results_list$water_bal[j,], type="l",  lty=1, xlab="Time", ylab="Soil active water content (%)", main=treatments[j], col="darkblue")
    plot(results_list$results_daily$date, results_list$soil_temp[j,], type="l",  lty=1, xlab="Time", ylab=expression(paste("Soil temperature (", degree,"C)")), main=treatments[j], col="firebrick")
    dev.off()

    filename3<-paste("reduction_factors_",j,".png", sep="")
    png(file=filename3, res=300, height=4500, width=2300)
    par(mfrow=c(3,1))
    plot(results_list$results_daily$date, results_list$re_wat[j,], type="l",  lty=1, xlab="Time", ylab="re_wat", main=treatments[j], col="darkblue", ylim=c(0,1))
    plot(results_list$results_daily$date, results_list$re_temp[j,], type="l",  lty=1, xlab="Time", ylab="re_temp", main=treatments[j], col="firebrick", ylim=c(0,1.1))
    plot(results_list$results_daily$date, results_list$re_crop[j,], type="l",  lty=1, xlab="Time", ylab="re_crop", main=treatments[j], col="darkgreen")
    dev.off()
  }

  } else {
    for(j in 1:length(treatments)){

    #plot the effect of the GAI function{
    levels<-levels(results_list$crop_id)
    date_range<-range(results_list$results_daily$date)
    plot(results_list$results_daily$date, results_list$GAI[j,], type="l",  lty=2, xlab="Time", ylab="GAI", main=treatments[j])
    for(i in 1:length(levels)){
      GAI_loop<-results_list$GAI[j,]
      GAI_loop[!results_list$crop_id==levels[i]]<-NA
      lines(results_list$results_daily$date, GAI_loop, col=crop_palette[i])
    }
    legend("topleft", levels, col=crop_palette, lty=1, pch=NA, bty="n")

    readline(prompt="Press [enter] to continue")
    par(mfrow=c(2,1))
    plot(results_list$results_daily$date, results_list$water_bal[j,], type="l",  lty=1, xlab="Time", ylab="Soil active water content (%)", main=treatments[j], col="darkblue")
    plot(results_list$results_daily$date, results_list$soil_temp[j,], type="l",  lty=1, xlab="Time", ylab=expression(paste("Soil temperature (", degree,"C)")), main=treatments[j], col="firebrick")

    readline(prompt="Press [enter] to continue")
    par(mfrow=c(3,1))
    plot(results_list$results_daily$date, results_list$re_wat[j,], type="l",  lty=1, xlab="Time", ylab="re_wat", main=treatments[j], col="darkblue", ylim=c(0,1))
    plot(results_list$results_daily$date, results_list$re_temp[j,], type="l",  lty=1, xlab="Time", ylab="re_temp", main=treatments[j], col="firebrick", ylim=c(0,1.1))
    plot(results_list$results_daily$date, results_list$re_crop[j,], type="l",  lty=1, xlab="Time", ylab="re_crop", main=treatments[j], col="darkgreen")
    readline(prompt="Press [enter] to continue")
    dev.off()
    }
  }

 cat("if you want to save the plots in your active folder remember to use the option: output_to_file=TRUE")
}




# datasets
#' template data used to test the reclim package
#'
#' @name templates
#' @docType data
#' @author Lorenzo Menichetti \email{ilmenichetti@@gmail.com}
#'
#'
#' @format two data frames.
#' Data are fictional, generated from real measured data but with a great deal of noise added.
#' First dataframe (daily climatic values), with 5114 rows and 9 variables:
#' \describe{
#'   \item{date}{date vector, YYYY-mm-dd}
#'   \item{year}{year, integer}
#'   \item{month}{month, integer}
#'   \item{day}{day, integer}
#'   \item{air_temp_deg_C}{mean daily air temperature (°C)}
#'   \item{precipitation_mm}{cumulated daily precipitation, in mm}
#'   \item{windspeed_kmh}{mean daily wind speed, in km/h}
#'   \item{humidity_percent}{mean daily air humidity, in \%}
#'   \item{Rsolar}{solar radiation, in this case in Langleys. Can be also cloudiness, in hours of sun per day}
#'   ...
#' }
#' Second dataframe (annual crop data), with 28 rows and 12 variables:
#' \describe{
#'   \item{year}{date vector, YYYY-mm-dd}
#'   \item{crop_description}{year, integer}
#'   \item{crop_id}{month, integer}
#'   \item{treat}{day, integer}
#'   \item{variance}{variance of the biomass distibution function, in days}
#'   \item{seeding}{seeding day, in day of year}
#'   \item{harvest}{first harves day, in day of year}
#'   \item{harvest2}{second harvest day (if present, for lay), in day of year}
#'   \item{tillage}{tillage day, in day of year}
#'   \item{minimum_cover}{minimum biomass on the ground all the time, kg of dry mass per ha}
#'   \item{total_dm_kg_ha}{aboveground biomass of first harvest, kg of dry mass per ha}
#'   \item{total_dm_kg_ha2}{aboveground biomass of second harvest (if present, for lay), kg of dry mass per ha}
#'   ...
#' }
#'
#' @examples
#' data(aboveground_testdata)
#' data(weather_testdata)
#'
NULL
