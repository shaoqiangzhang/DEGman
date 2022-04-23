DistrOfFit<-function(count,distribution){
## distribution is one of "poisson","nb","zinb",and "empty"
##output: a vector of probability for c(0:maxcount) 
  requiredPackages = c('pscl','fitdistrplus')
  for(p in requiredPackages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
  }
  if(sum(count) == 0){
    count <- append(count,1)## in order to fit NB for a zero-vector
  }
  data = data.frame(count)
  maxcount=2*max(round(count))
  if (distribution == "poisson") {
    model = glm("count ~ 1", data, family = "poisson")
    mu = exp(model$coefficients)
    x=c(0:maxcount)
    prob=dpois(x, lambda = mu)
  } else if (distribution == "nb") {
    #model = suppressWarnings(glm.nb("count ~ 1", data))
    model=suppressWarnings(glm.nb("count ~ 1",data,control = glm.control(maxit = 10, epsilon = 1e-8)))
    size = model$theta
    mu = exp(model$coefficients)
    x=c(0:maxcount)
    prob= dnbinom(x, size = size, mu = mu)
  } else if (distribution == "zinb") {
    fitResult = fitZINB(count)
    model = fitResult$modelZINB
    zeroProb = fitResult$zeroProb
    size = model$theta
    mu = exp(model$coefficients$count)
    x=c(1:maxcount)
    prob= dnbinom(x, size = size, mu = mu)## use NB to replace ZINB to eliminate the effect of false zeros 
    #prob=c(zeroProb +(1 - zeroProb) * pnbinom(0, size = size, mu = mu), (1 - zeroProb) * dnbinom(x, size = size, mu = mu))
  } else if(distribution == "empty"){
    maxexp=max(round(count))+1
    prob<-replicate(maxexp,0)
    for (i in 1: length(count)) {
    ##Frequency of expressionz values is counted
      prob[count[i]+1]=prob[count[i]+1]+1
    }
    prob=prob/sum(prob) 
  }
  return(prob)
}