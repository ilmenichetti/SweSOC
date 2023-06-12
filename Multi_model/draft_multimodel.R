
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
      porosity<-c()
      wilting_point<-c()
      field_capacity<-c()

      for (k in 1:treat_N){

      # #porosity
      # porosity[k]<-poros(sand=sites$sand/100, clay=sites$clay/100, SOC=mean(SOC[SOC$treat==k,]$SOC))

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

      # Soil water balance
      water_bal[[k]]<-waterbalance(twilt=wilting_point[k],
                              tfield=field_capacity[k],
                              ET0=PET_calc$ET0,
                              precipitation=weather[[j]]$precipitation_mm,
                              GAI=GAI_calc[[k]]$GAI,
                              date=GAI_calc[[k]]$date,
                              L=depth*10)

      ## Soil temperature
      soilT[[k]]<-soiltemp(temperature=weather[[j]]$air_temp_deg_C,
                      L=depth*10,
                      GAI=GAI_calc[[k]]$GAI)


      }

  site_weather_data[[j]]<-list(water_bal, soilT, PET_calc)

} #closing the loop for sites



# Save the weather data object to a file
saveRDS(site_weather_data, file = "./Dataset/site_weather_data.rds")








### ICBM with Bayesiantools
library(BayesianTools)
library(SoilR)
library(hydroGOF)
library(imputeTS)

# read the estimated climatic site data
site_weather_data<-readRDS(file = "./Dataset/site_weather_data.rds")

# the data structure is a list of lists:
# site_data[[j][[s]][[k]]
# where j is the site ID and s is the type (water or temperature)



#studying briefly the various scaling functions
str(site_weather_data[[j]][[3]])

fT.Arrhenius(site_weather_data[[j]][[2]][[k]]-273.15)
fT.Daycent2(site_weather_data[[j]][[2]][[k]])
fT.LandT(site_weather_data[[j]][[2]][[k]])
fT.linear(site_weather_data[[j]][[2]][[k]])

sites<-sites_full[sites_full$site_ID==j,]



gomperts_moist<-fW.Gompertz((site_weather_data[[j]][[1]][[k]][,1])/(depth*10))
moyano_moist<-fW.Moyano(site_weather_data[[j]][[1]][[k]][,1]/(depth*10))
candy_moist<-fW.Candy(site_weather_data[[j]][[1]][[k]][,1]/(depth*10), PV=poros(sand=sites$sand/100, clay=sites$clay/100, SOC=mean(treatments[treatments$treat==k,]$SOC_ave)))
century_moist<<-fW.Century(PPT=weather[[j]]$precipitation_mm, PET=site_weather_data[[j]][[3]][,2]) #Century has one single value for all treatments
rothc_moist<-fW.RothC(P=weather[[j]]$precipitation_mm, E=site_weather_data[[j]][[3]][,2], S.Thick = 20, pClay = sites$clay, pE = 1, bare = FALSE)$b #Century has one single value for all treatments


gomperts_moist_standard<-fW.Gompertz(0.5)
moyano_moist_standard<-fW.Moyano(0.5)
candy_moist_standard<-fW.Candy(0.5, 0.5)



gomperts_moist_standard_long<-fW.Gompertz(seq(0, 1, by=0.01))
moyano_moist_standard_long<-fW.Moyano(seq(0, 1, by=0.01))
candy_moist_standard<-fW.Candy(seq(0, 1, by=0.01), 0.5)

plot(gomperts_moist_standard_long, type="l", ylim=c(0,1))
lines(moyano_moist_standard_long)
lines(candy_moist_standard)

dates_vec<-as.Date(weather[[j]]$date)
year(dates_vec)

gomperts_moist_ave<-aggregate(gomperts_moist, list(year(dates_vec)), FUN=mean)
moyano_moist_ave<-aggregate(moyano_moist, list(year(dates_vec)), FUN=mean)
candy_moist_ave<-aggregate(candy_moist, list(year(dates_vec)), FUN=mean)
century_moist_ave<-aggregate(century_moist, list(year(dates_vec)), FUN=mean,  na.rm=TRUE)
rothc_moist_ave<-aggregate(rothc_moist, list(year(dates_vec)), FUN=mean)

zscores_moists<-data.frame(gomperts_moist_ave$x, moyano_moist_ave$x, candy_moist_ave$x ,century_moist_ave$x,
                         rothc_moist_ave$x)

names(zscores_moists)<-c("gompertz", "moyano", "candy", "century", "rothc")


temp_test=seq(-10,50)
standardtemp=20

arrhenius_norm=fT.Arrhenius(temp_test-273.15)/fT.Arrhenius(standardtemp-273.15)
daycent2_norm=fT.Daycent2(temp_test)/fT.Daycent2(standardtemp)
LloydTaylor_norm=fT.LandT(temp_test)/fT.LandT(standardtemp)
kirschbaum_norm=fT.KB(temp_test)/fT.KB(standardtemp)
linear_norm=fT.linear(temp_test)/fT.linear(standardtemp)
RothC_norm=fT.RothC(temp_test)/fT.RothC(standardtemp)
Q10_norm=fT.Q10(temp_test)/fT.Q10(standardtemp)
standcarb_norm=fT.Standcarb(temp_test)/fT.Standcarb(standardtemp)
demeter_norm=fT.Demeter(temp_test)/fT.Demeter(standardtemp)
century2_norm=fT.Century2(temp_test)/fT.Century2(standardtemp)
century1_norm=fT.Century1(temp_test)/fT.Century1(standardtemp)
scaled_temps<-data.frame(linear_norm, RothC_norm, kirschbaum_norm ,century1_norm,
                         century2_norm, daycent2_norm,LloydTaylor_norm,  Q10_norm,
                         standcarb_norm, demeter_norm, century1_norm, arrhenius_norm)


dim(temps)
names(temps)

png("temp_functions.png", width = 3000, height = 2000, res=300)
par(mfrow=c(2,2))
range_temp=range(scaled_temps[,1:10], na.rm = T)
plot(temp_test, scaled_temps[,1], type="l", ylim=range_temp, ylab="Standardized temp reduction (1 at 20\u00B0C)", xlab="Temperature (\u00B0C)", col=1, lty=1)
for(i in 2:3){
lines(temp_test, scaled_temps[,i], col=i, lty=i)
}
legend("topleft", names(temps)[1:3], lty=c(1:3), col=c(1:3), bty="n")

plot(temp_test, scaled_temps[,4], type="l", ylim=range_temp, ylab="Standardized temp reduction (1 at 20\u00B0C)", xlab="Temperature (\u00B0C)", col=4, lty=4)
lines(temp_test, scaled_temps[,1], col=1, lty=1)
for(i in 5:6){
  lines(scaled_temps[,i], col=i, lty=i)
}
legend("topleft", names(temps)[4:6], lty=c(4:6), col=c(4:6), bty="n")

plot(temp_test, scaled_temps[,7], type="l", ylim=range_temp, ylab="Standardized temp reduction (1 at 20\u00B0C)", xlab="Temperature (\u00B0C)", col=7,  lty=7)
lines(temp_test, scaled_temps[,1], col=1, lty=1)
for(i in 8:9){
  lines(scaled_temps[,i], col=i, lty=i)
}
legend("topleft", names(temps)[7:9], lty=c(7:9), col=c(7:9), bty="n")

plot(temp_test, scaled_temps[,10], type="l", ylim=range_temp, ylab="Standardized temp reduction (1 at 20\u00B0C)", xlab="Temperature (\u00B0C)", col=10,  lty=10)
lines(temp_test, scaled_temps[,1], col=1, lty=1)
for(i in 11:11){
  lines(scaled_temps[,i], col=i, lty=i)
}
legend("topleft", names(temps)[10:11], lty=c(10:11), col=c(10:11), bty="n")
dev.off()



png("moist_functions_zscore.png", width = 2000, height = 1500, res=300)
library(datawizard)
plot(gomperts_moist_ave$Group.1, standardize(gomperts_moist_ave$x), type="l", ylim=c(-4,4), ylab="Z-score standardized function value", xlab="year")
lines(moyano_moist_ave$Group.1, standardize(moyano_moist_ave$x), col=2, lty=2)
lines(candy_moist_ave$Group.1, standardize(candy_moist_ave$x), col=3, lty=3)
lines(century_moist_ave$Group.1, standardize(century_moist_ave$x),  col=4, lty=4)
lines(rothc_moist_ave$Group.1, standardize(rothc_moist_ave$x), col=5, lty=5)
legend("topleft", names(zscores_moists), lty=c(1:5), col=c(1:5), bty="n")

dev.off()


png("moist_functions_standard.png", width = 2500, height = 1500, res=300)
gomperts_moist_standard_long<-fW.Gompertz(seq(0, 1, by=0.01))
moyano_moist_standard_long<-fW.Moyano(seq(0, 1, by=0.01))
candy_moist_standard<-fW.Candy(seq(0, 1, by=0.01), 0.7)

plot(gomperts_moist_standard_long, type="l", ylim=c(0,1), ylab="Activity", xlab="Theta")
lines(moyano_moist_standard_long, col=2, lty=2)
lines(candy_moist_standard, col=3, lty=3)
legend("bottomright", names(zscores_moists)[1:3], lty=c(1:3), col=c(1:3), bty="n")
dev.off()


#working on the model(s)

param=c(k1 = 0.8,
           k2 = 0.00605,
           h = 0.13,
           r = 1.32)



#function to run the SOC model
runModel <- function(param, time,
                     c0 = c(Y0 = 0.3, O0 = 3.96), In=1){
  model =  ICBMModel(t=time, ks=c(k1=param[1], k2=param[2]),
                     h=param[3], r=param[4],
                     c0 = c0, In=In)
  predicted=rowSums(getC(model))
  return(predicted)
}



#likelihood function, calculated over all the experimental units
likelihood = function(param){

  predicted_long<-c()
  observed_long<-c()

  #running the model for one site
  for(j in 1:length(unique(SOC_full$site_ID))){ #j is the site ID
    for (k in 1:length(unique(SOC[SOC_full$site_ID==j,]$treat))){# k is the treatment ID
      #selecting the subset for the likelihood evaluation
      observed_ds=SOC[SOC_full$site_ID==j & SOC_full$treat==k,]

    observed=observed_ds$SOC
    time = observed_ds$year
    predicted = runModel(param, c0 = c(Y0 = 3000, O0 = 40000), In=1000, time = time)

    predicted_long<-c(predicted_long, predicted)
    observed_long<-c(observed_long, observed)
    } #end if k loop
  }#end of j loop

  #calculating the likelihood for one site
  ll = likelihoodIidNormal(predicted = predicted_long, observed = observed_long, sd=stdev)
  return(ll)
  # nse = NSE(sim = predicted_long, obs = observed_long)
  # return(nse)
}

# running the samples
setup = createBayesianSetup(likelihood,
                            lower = c(0.1,0.001,0.1,0.9), upper = c(1,0.01,0.5,1.1), names=c("k1", "k2", "h", "r"))
out = runMCMC(setup, sampler = "DEzs", settings = NULL)


chain<-out$chain


str(param)
str(chain[1,1:4][[1]])

time=seq(0,10)

runModel(chain[10,1:4][[1]], c0 = c(Y0 = 3000, O0 = 40000), In=1000, time = time)

sample<-getSample(out)
which(sample[,1]<=0)
which(sample[,2]<=0)
which(sample[,3]<=0)
which(sample[,4]<=0)
summary(out)
getVolume(out)

plot(out, start = 1000)

ML<-marginalLikelihood(out)

str(ML)

getPredictiveIntervals(
  sample,
  runModel,
  numSamples = 1000,
  quantiles = c(0.025, 0.975),
  error = NULL
)



#ICBM with Stan
#https://mc-stan.org/users/documentation/case-studies/soil-knit.html
library(rstan);
fit <- stan("ICBM.stan",
            data=c("totalC_t0", "t0", "N_t", "ts", "eCO2mean"),
            control=list(adapt_delta=0.90,
                         stepsize=0.005),
            chains=2, iter=200, seed=1234);






#Bare fallow example

library(FME)
library(SoilR)
library(Metrics)

LTBF_Ult_SOC<-SOC[SOC$treat==1,]




func <- function(pars) {
  t=LTBF_Ult_SOC$year-LTBF_Ult_SOC$year[1]
  C0=LTBF_Ult_SOC$SOC[1] # mass
  Model1=OnepModel(t,k=pars[1],C0,In=pars[2])
  Ct1=getC(Model1)
  model_cost=rmse(Ct1, LTBF_Ult_SOC$SOC)
  return(model_cost)
}

p = c(0.1, 2)
prior=function(p){p}
prior(p)

# simple Metropolis-Hastings
MCMC <- modMCMC(f = func,
                p = c(0.1, 2),
                niter = 500,
                outputlength = 1000, jump = 0.5)

MCMC <- modMCMC(f = func,
                p = c(0.1,2),
                niter = 500,
                outputlength = 1000, jump = 0.5, prior=log(p))


MCMC <- modMCMC(f = func,
                p = c(1),
                niter = 500,
                outputlength = 100, jump = 0.5)

str(MCMC)
summary(MCMC)
plot(MCMC) # check convergence
pairs(MCMC)


NL <- function(p) {
  mu <- 1:4
  -2*sum(log(dlnorm(p, mean = mu, sd = 1)))      #-2*log(probability)
}
MCMCl <- modMCMC(f = NL, p = log(1:2), niter = 3000,
                 outputlength = 1000, jump = 5)
plot(MCMCl)   # bad convergence
cumuplot(as.mcmc(MCMCl$pars))

