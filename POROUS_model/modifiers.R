
f_text_mic<-function(clay, phi_min){

  psi_w=-150 #meters, the wilting point pressure head
  theta_w=0.004+0.5*clay #psi_w (m) is the wilting point pressure head (= -150 m) and the corresponding water content theta_w (m3 m-3) is estimated from a pedotransfer function (Ostovari et al., 2015)
  lamda_mat=log(theta_w/phi_min)/log(-0.3/psi_w)

  f_text_mic=(-0.3/-6)^lamda_mat
  return(f_text_mic)
}

Delta_z<-function(f_agg=3,
                  Delta_z_min,
                  My_mic, Mo_mic, My_mes, Mo_mes,
                  phi_mat,
                  gamma_o){

  Mso=My_mic + Mo_mic + My_mes + Mo_mes
  Delta_z=(((1+f_agg)*(Mso/gamma_o))+Delta_z_min)/(1-phi_mac)

  return(Delta_z)
}



phi_mic<-function(My_mic, Mo_mic, My_mes, Mo_mes,
                  gamma_o, #density of organic matter
                  f_agg=3, #default suggested by Nick Jarvis
                  clay,
                  Delta_z_min, #soil layer thickess (m)
                  phi_min=1, #minimum matrix porosity, STILL MISSING
                  phi_mat
){

  f_text_mic_calc=f_text_mic(clay=clay, phi_min=phi_min)

  Delta_z_calc=Delta_z(f_agg=3,
                       Delta_z_min=Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)

  phi_mic=((f_agg*((My_mic+Mo_mic)/gamma_o))+f_text_mic_calc*Delta_z_min*phi_min)/Delta_z_calc
  return(phi_mic)
}





phi_mat<-function(phi_mic,
                  phi_mat,
                  f_agg=3,
                  gamma_o,
                  phi_min=1, #minimum matrix porosity, STILL MISSING
                  Delta_z_min,#soil layer thickess (m)
                  My_mic, Mo_mic, My_mes, Mo_mes){

  Mso=My_mic + Mo_mic + My_mes + Mo_mes

  Delta_z_calc=Delta_z(f_agg=3,
                       Delta_z_min = Delta_z_min,
                       My_mic, Mo_mic, My_mes, Mo_mes,
                       phi_mac,
                       gamma_o)


  phi_mat=(f_agg*(Mso/gamma_o)+Delta_z_min*phi_min)/Delta_z_calc
  return(phi_mat)
}


pore_frac<-function(phi_mac=0.2,
                    clay=0.2,
                    Delta_z_min=20,
                    gamma_o=1.2){
      phi_mat_calc<-phi_mat(phi_mic=phi_mic_calc, My_mic=1, Mo_mic=2, My_mes=1, Mo_mes=10, gamma_o=gamma_o, Delta_z_min = Delta_z_min, phi_mat=phi_mat)
      phi_mic_calc=phi_mic(My_mic=1, Mo_mic=2, My_mes=1, Mo_mes=10, gamma_o=gamma_o, clay=clay, Delta_z_min = Delta_z_min, phi_mat=phi_mat)
      phi_mes_calc=phi_mat_calc-phi_mic_calc

      mes_f=phi_mes_calc/(phi_mes_calc+phi_mic_calc)
      mic_f=phi_mic_calc/(phi_mes_calc+phi_mic_calc)
      return(c(mes_f, mic_f))
}
