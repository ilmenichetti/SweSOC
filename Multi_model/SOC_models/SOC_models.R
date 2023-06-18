
ICBM_desolve<-function(Is, Ir, Ia, re, parameters, init, time_sim, sim_steps=1){

  require(deSolve)

  Is_fun <- approxfun(Is, rule = 2)
  Ir_fun <- approxfun(Ir, rule = 2)
  Ia_fun <- approxfun(Ia, rule = 2)
  re_fun <- approxfun(re, rule = 2)

  ##ODE Porous for variable inputs
  ODE_ICBM <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      Is=Is_fun(t)
      Ir=Ir_fun(t)
      Ia=Ia_fun(t)
      re=re_fun(t)

      .Ys=Is-ky*re*Ys
      .Yr=Ir-ky*re*Yr
      .Ya=Ia-ky*re*Ya

      .O = hs*ky*re*Ys + hr*ky*re*Yr + ha*ky*re*Ya - ko*re*O

      return(list(c(.Ys, .Yr, .Ya, .O)))
    })
  }


constants=parameters

sim_time<-seq(0,time_sim,by=sim_steps)

out <- ode(y = init, time = sim_time, func = ODE_ICBM, parameters)
#out <- ode(y = c(Yr=1000,Ys=1000,Ya=1000,O=1000), time = sim_time, func = ODE_ICBM, constants)

SOC_tot=rowSums(out[,2:5])

return(as.data.frame(cbind(out, SOC_tot)))

}
