


allom_linear<-function(yields, A, R, E){

  aboveground=yields*A
  roots=yields*R
  exudates=yields*E

  return(list(aboveground, roots, exudates))
}



allom_intercept<-function(yields, A, a, R, r, E, e){

  aboveground=a+yields*A
  roots=r+yields*R
  exudates=e+yields*E

  return(list(aboveground, roots, exudates))

}