pzinb = function(x, mu, size, zeroProb) {
  # calculate the probability of zero inflated negative binomial distribution

  # input
  # x: a vector to evaluate the probability
  # mu: the mu of the NB component
  # size: the size of the NB component, it is 1 / dispersion
  # zeroProb: the probability of the zero component

  prob = NA
  prob = zeroProb +
    (1 - zeroProb) * pnbinom(x, size = size, mu = mu)
  prob[x < 0] = 0
  return(prob)
}

fitZINB = function(count) {
  # use NB fit as an initial fit and give it to zeroinfl for further fitting
  # assume no offset for all counts

  # input
  # count: a vector of counts

  # output
  # a list of three components:
  #   fitted: the fitted probablity
  #   diffLogLikZINBNB: the difference of loglik between zinb and nb
  #   converge: if succeed, it is "converged", else other messages
  #   modelZINB: the zero inflated model

  fitted = NA
  diffLogLikZINBNB = NA
  converge = "failed"
  modelZINB = NULL
  zeroProb = NA
  # maker sure there is at least one 0
  if (min(count) != 0) {
    converge = "no zero count"
    return(list(fitted = fitted, diffLogLikZINBNB = diffLogLikZINBNB,
                converge = converge, modelZINB = modelZINB,
                zeroProb = zeroProb))
  }

  modelNB = list()
  modelZINB = list()

  modelNB$converged = NA
  modelZINB$converged = NA

  # fit with NB
  logLikNB = NA
  countPerCell = numeric(length(count))
  countPerCell[ ] = 1
  try({
    modelFormula = formula(paste("count ~ 1"))
    modelFrame = data.frame(rep(1, length(count)))
    result0 = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                        countPerCell, "pooled-ML")
    dispersionCR = result0$dispersion
    logLikNB = result0$logLik
    # browser()
    beta0CR = coefficients(result0$model)[1]
    modelNB = result0$model
  })
  thetaNB = 1 / dispersionCR

  try({
    if (!is.na(dispersionCR)) {
      EM = F
      start = list(count = coefficients(modelNB), zero = -3, theta = thetaNB)
    } else {
      EM = F
      start = NULL
    }
    modelZINB = suppressWarnings(
                   zeroinfl(formula(paste("count ~ 1 | 1")),
                         dist = "negbin",
                         EM = EM,
                         start = start)
    )
  })
  if (is.na(modelZINB$converged)) {
    try({
      modelZINB = suppressWarnings(
          zeroinfl(formula(paste("count ~1 | 1")),
                           dist = "negbin",
                           EM = T)
      )
    })
  }

  diffLogLikZINBNB = logLik(modelZINB) - logLikNB
  if (!is.na(modelZINB$converged) && !is.na(diffLogLikZINBNB) &&
      diffLogLikZINBNB > -0.5) {
    converge = "converged"
    thetaZINB = modelZINB$theta
    logthetaSE = modelZINB$SE.logtheta
    zeroProb = 1 / (1 + exp(-modelZINB$coefficients$zero))
    zeroPValue = summary(modelZINB)$coefficients$zero[4]
  }

  return(list(fitted = fitted, diffLogLikZINBNB = diffLogLikZINBNB,
              converge = converge, modelZINB = modelZINB,
              zeroProb = zeroProb))
}

combineTable = function(value, count, threshold) {
  # combine the table to make sure the count in each cell is above a threshold,
  # e.g., 5
  
  # input
  # value: the value 
  # count: the count of the observed value
  # threshold: a threshold value
  
  # output
  # combined: a 3 * n matrix with count >= threshold
  
  combined = NULL
  pointer = length(value)
  currentBucket = c(NA, Inf, 0) # (leftBound, rightBound, 0)
  while (pointer >= 1) {
    # add to the current bucket
    currentBucket[3] = currentBucket[3] + count[pointer]
    # check whether the bucket count >= threshold
    if (currentBucket[3] >= threshold || pointer == 1) {
      # set the left bound
      currentBucket[1] = value[pointer]
      # add to the combined
      combined = cbind(combined, currentBucket)
      # reset the currentBucket
      currentBucket = c(NA, value[pointer] - 1, 0)
    }
    pointer = pointer - 1
  }
  
  nBin = dim(combined)[2]
  
  # change the last bin to start from 0
  combined[1, nBin] = 0
  
  # make sure the last bin is >= threshold
  if (combined[3, nBin] < threshold) {
    # combine bin 1 and 2
    if (dim(combined)[2] > 1) {
      combined[1, nBin - 1] = combined[1, nBin]
      combined[3, nBin - 1] = combined[3, nBin] + combined[3, nBin - 1]
      combined = combined[, -nBin] # remove the last
    } else {
      print("not enough count with all combined")
    }
  }
  
  return(combined)
}



goodnessOfFit = function(count, distribution) {
  requiredPackages = c('pscl','fitdistrplus')
  for(p in requiredPackages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
  }
  if(sum(count) == 0){
   count <- append(count,1)## in order to fit NB for a zero-vector
  }
  # fit the count data with the distribution
  fitResult = NULL
  data = data.frame(count)
  size = NULL
  mu = NULL
  model = NULL
  pvalue = NA
  converge = "converged"
  if (distribution == "poisson") {
    model = glm("count ~ 1", data, family = "poisson")
    if (!model$converge) {
      converge = model$converge
    }
    mu = exp(model$coefficients)
  } else if (distribution == "nb") {
    #model = suppressWarnings(glm.nb("count ~ 1", data))
    model=suppressWarnings(glm.nb("count ~ 1",data,control = glm.control(maxit = 10, epsilon = 1e-5)))
    size = model$theta
    if (!model$converge) {
      converge = model$converge
    }
    mu = exp(model$coefficients)
  } else if (distribution == "zinb") {
    fitResult = fitZINB(count)
    if (fitResult$converge != "converged") {
      converge = fitResult$converge
      return(cbind(pvalue, converge))
    }
    model = fitResult$modelZINB
    zeroProb = fitResult$zeroProb
    size = model$theta
    mu = exp(model$coefficients$count)
  } else {
    stop("wrong distribution specified")
  }

  # group data
  tableResult = table(count)
  value = as.numeric(names(tableResult))
  countOfValue = as.vector(tableResult)
  groupBin = combineTable(value, countOfValue, 1)

  # Check the fitting with Chi-sq Goodnees of feet test
  testStat <- c()
  testDf <- c()
  TestPvalue <- c()

  #print(i)
  observed = groupBin[3,]
  expected = numeric(length(observed))
  expected[ ] = NA
  for (j in 1:ncol(groupBin)) {
    up = groupBin[2,j]
    low = groupBin[1,j] - 1
    if (distribution == "poisson") {
      upP = ppois(up, lambda = mu)
      lowP = ppois(low, lambda = mu)
    }
    if (distribution=="nb") {
      upP = pnbinom(up, size = size, mu = mu)
      lowP = pnbinom(low, size = size, mu = mu)
    } else if (distribution == "zinb") {
      upP = pzinb(up, size = size, mu = mu, zeroProb = zeroProb)
      lowP = pzinb(low, size = size, mu = mu, zeroProb = zeroProb)
    }
    expected[j] = upP-lowP
  }
  expected[abs(expected)<1e-15] = 1e-15
  if(length(observed)<2){
    testDfAdj = 0
  }else{
    chiqTest = suppressWarnings(chisq.test(x=observed, p=expected, correct = F))
    testStat = chiqTest$statistic
    testDf = chiqTest$parameter

    testDfAdj = NULL
    if (distribution=="poisson") {
      testDfAdj <- testDf-1
    } else if (distribution=="nb") {
      testDfAdj <- testDf-2
    } else if (distribution == "zinb") {
      testDfAdj <- testDf-3
    }
  }

  if (testDfAdj <= 0) {
    #warning("The number of bins is not enough for the test")
    converge = paste("dof is", testDfAdj)
    pvalue = NA
  } else {
    pvalue <- pchisq(testStat, df=testDfAdj, lower.tail = FALSE)
  }

  return(data.frame(pvalue, converge, stringsAsFactors = F))
}

bestFitModel = function(count) {
  # run goodness of fit on multiple distributions
  # input: count: a count vector
  # output: one of "poisson", "nb", "zinb", and "empty"
  dist=c("poisson", "nb", "zinb")
  goodnessOfFitResult <-c()
  for (j in 1:3) {
    result = goodnessOfFit(count, dist[j])
    if (result[, 2] == "converged" ) {
      goodnessOfFitResult<-c(goodnessOfFitResult,result[, 1]) 
    } else {
      goodnessOfFitResult<-c(goodnessOfFitResult,-1.0) ##if not fit, pvalue=-1
    }
  }
  if(max(goodnessOfFitResult)<0){
    bestmodel="empty"
  }else{
    bestmodel=dist[which.max(goodnessOfFitResult)]
  }
  return(bestmodel)
}