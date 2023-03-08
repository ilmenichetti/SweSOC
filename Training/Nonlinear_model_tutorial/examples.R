library(SoilR)



Porous<-function(ky=0.8, ko=0.00605,
                 kmix=0.9,
                 e=0.13,
                 Im=1.1, Ir=0.5,
                 F_prot=0.0){
  time_symbol='t'

  ##### IN
  ifs=SoilR:::InFluxList_by_PoolName(
    c(
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mes',
        func=function(t){
          Im
        }
      ),
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mes',
        func=function(t){
          Ir*(0.5)
        }
      ),
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mic',
        func=function(t){
          Ir*(0.5)
        }
      )
    )
  )
  ##### OUT
  ofs=SoilR:::OutFluxList_by_PoolName(
    c(
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mes',
        func=function(My_mes){
          ky*My_mes
        }
      )
      ,
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mes',
        func=function(Mo_mes){
          (1-e)*ko*Mo_mes
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mic',
        func=function(My_mic){
          ky*F_prot*My_mic
        }
      ),
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mic',
        func=function(Mo_mic){
          (1-e)*ko*F_prot*Mo_mic
        }
      ),
      # Bioturbation, nonlinearity here, Ty and To are outfluxes from the micropores (and also internal fluxes below, into My_mes and Mo_mes)
      #Ty
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mic',
        func=function(My_mic,My_mes){
          kmix*((My_mic-My_mes)/2)
        }
      ),
      #To
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mic',
        func=function(Mo_mic,Mo_mes){
          kmix*((Mo_mic-Mo_mes)/2)
        }
      )
    )
  )

  ##### INT
  intfs=SoilR:::InternalFluxList_by_PoolName(
    list(
      SoilR:::InternalFlux_by_PoolName(
        sourceName='My_mes',
        destinationName='Mo_mes',
        func=function(My_mes){
          e*ky*My_mes
        }
      ),
      SoilR:::InternalFlux_by_PoolName(
        sourceName='My_mic',
        destinationName='Mo_mic',
        func=function(My_mic){
          e*ky*F_prot*My_mic
        }
      ),
      #bioturbation fluxes
      #Ty
      SoilR:::InternalFlux_by_PoolName(
        sourceName='My_mic',
        destinationName='My_mes',
        func=function(My_mic, My_mes){
          kmix*((My_mic-My_mes)/2)
        }
      ),
      #To
      SoilR:::InternalFlux_by_PoolName(
        sourceName='Mo_mic',
        destinationName='Mo_mes',
        func=function(Mo_mic, Mo_mes){
          kmix*((Mo_mic-Mo_mes)/2)
        }
      )
    )
  )

  smod <- SoilR:::SymbolicModel_by_PoolNames(
    in_fluxes=ifs,
    internal_fluxes=intfs,
    out_fluxes=ofs,
    timeSymbol=time_symbol
  )
  smod
}


##### building the functions for porosity

phi_mac=0.2
clay=0.2
Delta_z_min=20
gamma_o=1.2


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


phi_mic_calc=phi_mic(My_mic=1, Mo_mic=2, My_mes=1, Mo_mes=10, gamma_o=gamma_o, clay=clay, Delta_z_min = Delta_z_min, phi_mat=phi_mat)


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

phi_mat_calc<-phi_mat(phi_mic=phi_mic_calc, My_mic=1, Mo_mic=2, My_mes=1, Mo_mes=10, gamma_o=gamma_o, Delta_z_min = Delta_z_min, phi_mat=phi_mat)


phi_mat_calc+phi_mac


modelObject<-Porous()
plotPoolGraph(modelObject)

iv<-c(My_mes=1, Mo_mes=10,My_mic=0.6, Mo_mic=3)
times<-seq(0,20,by=0.1)
modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=iv)

Ct0<-getC(modrun0)
Rt0<-getReleaseFlux(modrun0)


#par(mfrow=c(2,1), mar=c(4,4,0,1))
matplot(times, Ct0, type="l", lty=1, col=1:4, xlab=" ", ylab="Pool contents", bty="n")
legend("topleft", c("My_mes", "Mo_mes", "My_mic", "Mo_mic"), lty=1, col=1:4, bty="n")
#matplot(times, Rt0,  type="l", lty=1, col=1:2, xlab="Time", ylab="Respiration", bty="n")