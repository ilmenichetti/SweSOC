library(SoilR)


TwopoolNonlinearInput<-function(ky, ko, h){
  time_symbol='t'

  ifs=SoilR:::InFluxList_by_PoolName(
    c(
      SoilR:::InFlux_by_PoolName(
        destinationName='Cy',
        func=function(t, Co){
          Iy
        }
      ),
      SoilR:::InFlux_by_PoolName(
        destinationName='Co',
        func=function(t){
          Io
        }
      )
    )
  )
  ofs=SoilR:::OutFluxList_by_PoolName(
    c(
      SoilR:::OutFlux_by_PoolName(
        sourceName='Cy',
        func=function(Cy){
          ky*Cy
        }
      )
      ,
      SoilR:::OutFlux_by_PoolName(
        sourceName='Co',
        func=function(Co){
          ko*Co
        }
      )
    )
  )
  intfs=SoilR:::InternalFluxList_by_PoolName(
    list(
      SoilR:::InternalFlux_by_PoolName(
        sourceName='Cy',
        destinationName='Co',
        func=function(Cy){
          h*ky*Cy
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

Iy=0.6
Io=12
iv<-c(Cy=1, Co=10)
duration=20
time_step=1

modelObject<-TwopoolNonlinearInput(ky=0.8, ko=0.00605, h=0.13, Iy, Io)

times<-seq(0,duration,by=time_step)
modrun0<-Model_by_PoolNames(smod=modelObject, times=times, initialValues=iv)
Ct0<-getC(modrun0)
Rt0<-getReleaseFlux(modrun0)

plot(Rt0[,1], type="l", ylim=c(0, max(Rt0)), col="blue")
lines(Rt0[,2], col="red")

### mass balance ###

#SOC difference
SOC_diff<-rowSums(Ct0)[length(times)]-rowSums(Ct0)[1]

#Respiration (cumulated)
RESP_tot<-sum(rowSums(Rt0*time_step))

#Inputs total
Input_tot<-sum(rep(Iy+Io, duration))

# check mass balance
Input_tot
(SOC_diff+RESP_tot)





ICBM_min<-function(ky, ko, h, Im){

  time_symbol='t'

  ##### IN
  ifs=SoilR:::InFluxList_by_PoolName(
    c(
      SoilR:::InFlux_by_PoolName(
        destinationName='Y',
        func=function(t){
          Im
        }
      )
    )
  )
  ##### OUT
  ofs=SoilR:::OutFluxList_by_PoolName(
    c(
      SoilR:::OutFlux_by_PoolName(
        sourceName='Y',
        func=function(Y){
          ky*Y
        }
      )
      ,
      SoilR:::OutFlux_by_PoolName(
        sourceName='O',
        func=function(O){
          ko*O
        }
      )
    )
  )

  ##### INT
  intfs=SoilR:::InternalFluxList_by_PoolName(
    list(
      SoilR:::InternalFlux_by_PoolName(
        sourceName='Y',
        destinationName='O',
        func=function(Y){
          h*ky*Y
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

init<-c(Y=1, O=10)

modelObject_minimal<-ICBM_min(ky=0.8, ko=0.00605, h=0.13, Im=Iy+Io)
plotPoolGraph(modelObject_minimal)

modrun0_minimal<-Model_by_PoolNames(smod=modelObject_minimal, times=times, initialValues=init)
Ct1<-getC(modrun0_minimal)
Rt1<-getReleaseFlux(modrun0_minimal)

plot(Rt1[,1], type="l", ylim=c(0, max(Rt1)), col="blue")
lines(Rt1[,2], col="red")

# checking the results (mass balance)
### mass balance

#SOC difference
SOC_diff<-rowSums(Ct1)[length(times)]-rowSums(Ct1)[1]
#Respiration (cumulated)
RESP_tot<-sum(rowSums(Rt1*time_step)[-1])

#Inputs total
Input_tot=(Iy+Io)*(duration)

# check mass balance
Input_tot
(SOC_diff+RESP_tot)


