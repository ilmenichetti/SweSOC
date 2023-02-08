library(SoilR)




FourpoolNonlinearInput<-function(ky=0.8, ko=0.00605,
                                 kmix=0.1,
                                 e=0.13,
                                 Iy=1.1, Io=0.5,
                                 F_prot=0.1){
  time_symbol='t'

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
      # Ty and To are outfluxes from the micropores (and also internal fluxes below, into My_mes and Mo_mes)
      #Ty
      SoilR:::OutFlux_by_PoolName(
        sourceName='My_mic',
        func=function(t,My_mic,My_mes){
          kmix*((My_mic-My_mes)/2)
        }
      ),
      #To
      SoilR:::OutFlux_by_PoolName(
        sourceName='Mo_mic',
        func=function(t,Mo_mic,Mo_mes){
          kmix*((Mo_mic-Mo_mes)/2)
        }
      )
    )
  )
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
        func=function(t, My_mic, My_mes){
          kmix*((My_mic-My_mes)/2)
        }
      ),
      #To
      SoilR:::InternalFlux_by_PoolName(
        sourceName='Mo_mic',
        destinationName='Mo_mes',
        func=function(t, Mo_mic, Mo_mes){
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



modelObject<-FourpoolNonlinearInput()
plotPoolGraph(modelObject)

iv<-c(My_mes=1, Mo_mes=10,My_mic=0.6, Mo_mic=3)
times<-seq(0,20,by=0.1)
modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=iv)

Ct0<-getC(modrun0)
Rt0<-getReleaseFlux(modrun0)


#par(mfrow=c(2,1), mar=c(4,4,0,1))
matplot(times, Ct0, type="l", lty=1, col=1:2, xlab=" ", ylab="Pool contents", bty="n")
legend("topleft", c("Cy", "Co"), lty=1, col=1:2, bty="n")
#matplot(times, Rt0,  type="l", lty=1, col=1:2, xlab="Time", ylab="Respiration", bty="n")