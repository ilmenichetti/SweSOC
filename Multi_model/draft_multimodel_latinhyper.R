
library(readODS)
library(viridis)

#reading the site data
sites_full<-read_ods("./Dataset/LTE_data.ods", sheet = "Sites")
treatments_full<-read_ods("./Dataset/LTE_data.ods", sheet = "Treatments")
yields_full<-read_ods("./Dataset/LTE_data.ods", sheet = "Yields")
amendments_full<-read_ods("./Dataset/LTE_data.ods", sheet = "Amendments")
SOC_full<-read_ods("./Dataset/LTE_data.ods", sheet = "SOC")

#a bit of preprocessing
treatments_full$site<-as.factor(treatments_full$site)
sites_names<-sites_full$site




#reading the weather files for every site
weather<-list()
for(i in 1:length(sites_names)){
weather[[i]]<-read_ods("./Dataset/LTE_data.ods", sheet = sites_names[i])
}


#plotting the sites SOC
figure_side<-ceiling(sqrt(length(sites_names)))
png("SOC_by_site.png", width = 1800*figure_side, height = 1600*figure_side, res=300)
par(mfrow=c(figure_side, figure_side))
for(i in 1:length(sites_names)){
  treatments_site<-treatments_full[treatments_full$site==sites_names[i],]$treat
  treatment_palette<-viridis(length(treatments_site))
  plot(SOC_full[SOC_full$site==sites_names[i] & SOC_full$treat==1,]$year, SOC_full[SOC_full$site==sites_names[i] & SOC_full$treat==1,]$SOC, type="l", ylim=c(0,100*1000),
       col=treatment_palette[1], ylab = expression(paste("SOC (kg C h", a^-1,")")), xlab="year", main=sites_names[i])
  for (j in 2:length(treatments_site)){
    lines(SOC_full[SOC_full$site==sites_names[i] & SOC_full$treat==j,]$year, SOC_full[SOC_full$site==sites_names[i] & SOC_full$treat==j,]$SOC, col=treatment_palette[j])
    }
  legend("topleft", as.character(treatments_site), col = treatment_palette, bty="n", lty=1)
  }
dev.off()



# Soil temperature and water balance calculations
source("./Multi_model/weather_functions/PET.R")
source("./Multi_model/weather_functions/GAI.R")
source("./Multi_model/weather_functions/accessory_functions.R")

#depth for the calculations
depth=20

#creating the list (of lists) where to store the results
site_weather_data<-list()

for(j in 1:length(sites_full$site_ID)){

  #selecting a subset of data to work with the loop with only site_ID==j
  sites<-sites_full[sites_full$site_ID==j,]
  treatments<-treatments_full[treatments_full$site==sites$site,]
  yields<-yields_full[yields_full$site==sites$site,]
  amendments<-amendments_full[amendments_full$site==sites$site,]
  SOC<-SOC_full[SOC_full$site==sites$site,]



      ## Soil water balance

      # PET
      # calculates one PET for each day, but only one site
      PET_calc<-PET(humidity=weather[[j]]$humidity,
                    windspeed=weather[[j]]$windspeed,
                    temperature=weather[[j]]$air_temp_deg_C,
                    latitude=sites$latitude,
                    altitude=sites$altitude,
                    sun=weather[[1]]$Rsolar,
                    sun.mode="Rsolar", # can be Rsolar, sunlight or cloudiness, latter two between 0 and 1
                    date=as.Date(weather[[j]]$date))



      # GAI
      # calculates one GAI for each day AND each treatment, but only one site

      treat_N<-length(unique(yields$treat))

      GAI_calc<-list()
      water_bal<-list()
      soilT<-list()
      GAI_calc_annual<-list()
      water_bal_annual<-list()
      soilT_annual<-list()
      theta_annual<-list()

      porosity<-c()
      wilting_point<-c()
      field_capacity<-c()

      for (k in 1:treat_N){

      # #porosity
      porosity[k]<-poros(clay=sites$clay/100, SOC=mean(treatments[treatments$treat==k,]$SOC_ave))


      #wilting point
      wilting_point[k]<-WP(sand=sites$sand/100, clay=sites$clay/100, SOC=mean(treatments[treatments$treat==k,]$SOC_ave))

      #field capacity
      field_capacity[k]<-FC(sand=sites$sand/100, SOC=mean(treatments[treatments$treat==k,]$SOC_ave))


      treatment<-unique(yields$treat)[k]
      selected_yields<-yields[yields$treat==treatment,]

      GAI_calc[[k]]<-GAI(yield=selected_yields$total_yield,
                    year=selected_yields$year,
                    crop=selected_yields$crop_id,
                    variance=selected_yields$variance,
                    seeding=selected_yields$seeding,
                    harvest=selected_yields$harvest,
                    tillage=selected_yields$tillage,
                    minimum_cover=selected_yields$minimum_cover)
      GAI_calc_daily=cbind(GAI_calc[[k]], year=year(GAI_calc[[k]]$date))
      GAI_calc_annual[[k]]<-ddply(GAI_calc_daily, c("year"), summarise,
                                   GAI_annual = mean(GAI),
                                   LAI_annual = mean(LAI))



      # Soil water balance
      water_bal[[k]]<-waterbalance(twilt=wilting_point[k],
                              tfield=field_capacity[k],
                              ET0=PET_calc$ET0,
                              precipitation=weather[[j]]$precipitation_mm,
                              GAI=GAI_calc[[k]]$GAI,
                              date=GAI_calc[[k]]$date,
                              L=depth*10)

      water_bal_daily=cbind(water_bal[[k]], year=year(water_bal[[k]]$date))
      water_bal_annual[[k]]<-ddply(water_bal_daily, c("year"), summarise,
                           water_annual = mean(water),
                           Eact_annual = mean(Eact))
      theta_annual[[k]]=data.frame(year=water_bal_annual[[k]]$year, theta=water_bal_annual[[k]]$water_annual/(field_capacity[k]*depth*10))

      ## Soil temperature
      soilT[[k]]<-soiltemp(temperature=weather[[j]]$air_temp_deg_C,
                      L=depth*10,
                      GAI=GAI_calc[[k]]$GAI)
      soilT_daily=as.data.frame(cbind(soilT[[k]], year=year(water_bal[[k]]$date)))
      soilT_annual[[k]]<-ddply(soilT_daily, c("year"), summarise,
                            soilT = mean(V1))

      }

  site_weather_data[[j]]<-list(water_bal, soilT, PET_calc, theta_annual, soilT_annual)

} #closing the loop for sites



# Save the weather data object to a file
saveRDS(site_weather_data, file = "./Dataset/site_weather_data.rds")









# read the estimated climatic site data
site_weather_data<-readRDS(file = "./Dataset/site_weather_data.rds")
j=1
k=1
# the data structure is a list of lists:
# site_data[[j][[s]][[k]]
# where j is the site ID and s is the type (water or temperature, s=4 water annual, s=5 temp annual)

source("./Multi_model/reduction_functions/reduction_functions.R")
source("./Multi_model/SOC_models/SOC_models.R")


#defining the cost function
ICBM_cost<-function(parameters, soilT, soilmoist, Is, Ir, Ia, measured_years, measured_SOC, init, sim_steps=1){

  require(deSolve)

    #weather scaling
    re_temp<-lloydtaylor_temp(soilT,A=1)
    re_moist<-gompertz_moist(soilmoist, c=1)
    re=re_temp*re_moist
    #running the simulation
    time_sim=length(re)
    predicted=ICBM_desolve(Is=Is, Ir=Is, Ia=Ia, re=re, parameters=parameters, init=init, time_sim)

    which_measurement<-(measured_years-measured_years[1])+1
    predicted_vec<-predicted$SOC[which_measurement]

    rmse<-sqrt(mean((measured_SOC[!is.na(predicted_vec)] - predicted_vec[!is.na(predicted_vec)])^2))
    return(rmse)
  }



#iterating over each site and treatment
#The models are judged based on the overall RMSE but the RMSE for each site and treatment is saved in a list of vectors j*k
RMSE_complete<-list()
for(j in 1:length(sites_names)){
  site_treatments<-treatments_full[treatments_full$site==sites_names[j],]$treat
  RMSE_site<-c()
    for (k in 1:length(site_treatments)){
    SOC_sub<-SOC_full[SOC_full$site==sites_names[j] & SOC_full$treat==k,]
    RMSE_site[k]<-ICBM_cost(soilT=site_weather_data[[j]][[5]][[k]]$soilT,
                  soilmoist=site_weather_data[[j]][[4]][[k]]$theta,
                  Is=rep(500, 54),
                  Ir=rep(500, 54),
                  Ia=rep(500, 54),
                  parameters=c(ky=0.8,
                                  ko=0.0065,
                                  hs=0.1,
                                  hr=0.1,
                                  ha=0.1),
                  init<-c(Yr=SOC_sub$SOC[1]*0.05,Ys=SOC_sub$SOC[1]*0.05,Ya=SOC_sub$SOC[1]*0.05,O=SOC_sub$SOC[1]*0.85),
                  measured_SOC=SOC_sub$SOC,
                  measured_years=SOC_sub$year)

    } #closing k loop
  names(RMSE_site)<-site_treatments
  RMSE_complete[[j]]<-RMSE_site

} #closing j loop
names(RMSE_complete)<-sites_names
