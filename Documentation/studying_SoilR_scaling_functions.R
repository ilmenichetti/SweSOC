

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
