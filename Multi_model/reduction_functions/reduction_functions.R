### Temperature

# Ratkovski temperature reduction function (ChatGPT, to check)
ratkovski_temp <- function(T, Topt, Tlow, Thigh, k) {
  if (T <= Tlow || T >= Thigh) {
    return(0)
  } else if (T > Tlow && T <= Topt) {
    return(exp(k * (T - Topt) / (Topt - Tlow)))
  } else {
    return(exp(k * (Topt - T) / (Thigh - Topt)))
  }
}

lloydtaylor_temp <- function(T, A) {
  # Temperature reduction function (Lloyd and Taylor 1994)
  # R: temperature reduction factor (-)

  # Parameters
  E_0 <- 308.56 #K
  T_0 <- 227.13 #K


  R <- A * exp(-E_0 / (T - T_0))


  return(R)
}




### Moisture
moyano_moist<-function(moisture, clay, SOC){

    # Moisture reduction function (Moyano et al. 2012)
    # moisture: soil moisture (m3/m3)
    # clay: clay content (g/kg)
    # SOC: soil organic carbon (g/kg)

    # Parameters
    beta1 <- -1.31
    beta2 <- 3.0
    beta3 <- -2.23
    beta4 <- 0.26
    beta5 <- -0.39
    beta6 <- 1.07

    # variance, not going to use it for now
    # beta1_var <- 0.05
    # beta2_var <- 0.2
    # beta3_var <- 0.2
    # beta4_var <- 0.02
    # beta5_var <- 0.05
    # beta6_var <- 0.07

    # Calculate moisture reduction function
    reduction <- beta1 * moisture + beta2 * moisture + beta3 * moisture + beta4 * clay + beta5 * clay + beta6 * SOC

    return(reduction)
    }


candy_moist<-function(moisture, clay, SOC, A=4){
  reduction<-c()

  for (i in 1:length(moisture)){
  # calculate porosity
  # according to Kätterer, T., O. Andrén, and P-E. Jansson. “Pedotransfer Functions for Estimating Plant Available
  # Water and Bulk Density in Swedish Agricultural Soils.” Acta Agriculturae Scandinavica, Section B - Plant Soil Science
  # 56, no. 4 (December 2006): 263–76. https://doi.org/10.1080/09064710500310170.
  PV = 0.3697 + clay + 0.0403 * SOC -0.0098 * clay * SOC
  # According to:Bauer, J., Herbst, M., Huisman, J.A., Weihermüller,
  # L., Vereecken, H., 2008. Sensitivity of simulated soil heterotrophic
  # respiration to temperature and moisture reduction functions.
  # Geoderma 145, 17–27. https://doi.org/10.1016/j.geoderma.2008.01.026
  if((moisture[i]/PV)<0.5){
    reduction[i]<-A*(moisture[i]/PV) * (1-(moisture[i]/PV))
  } else {
    reduction[i]<-1
  }

  }

  #reduction <- fW.Candy(moisture, PV)
  return(list(reduction=reduction, porosity=PV))

  }




gompertz_moist <- function(moisture, c) {
  a=1
  b=1
  exp_val <- exp(-b * exp(-c * moisture))
  y <- a * exp_val
  normalized_y <- (y - min(y)) / (max(y) - min(y))
  return(normalized_y)
}

