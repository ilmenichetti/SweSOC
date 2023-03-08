

require(SoilR)


Porous<-function(ky=0.8, ko=0.00605,
                 kmix=0.9,
                 e=0.13,
                 Im=1.1, Ir=0.5,
                 F_prot=0.0,
                 phi_mac=0.2,
                 clay=0.2,
                 Delta_z_min=20,
                 gamma_o=1.2){

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
          Ir*pore_frac(phi_mac, clay, Delta_z_min, gamma_o)[1]
        }
      ),
      SoilR:::InFlux_by_PoolName(
        destinationName='My_mic',
        func=function(t){
          Ir*pore_frac(phi_mac, clay, Delta_z_min, gamma_o)[2]
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




