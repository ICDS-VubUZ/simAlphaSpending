alphaSpendingFunction <- function(tau, alpha, type = "OBF"){
  
  if(type == "OBF"){
    #alpha spending function approximating O'Brien-Fleming boundaries
    result <- 2 - 2*pnorm(qnorm(1-alpha/2)/sqrt(tau))
    
  }else if(type == "Pocock"){ 
    #alpha spending function approximating Pocock boundaries
    result <- alpha*log(1+(exp(1)-1)^tau)
    
  }else if(type == "compromise"){
    result <- alpha*tau
  }
  
  return(result)
}

simIntAtTestAB <- function(cohensD, Ns, alpha, beta, NrOfSimulations, sides = 1, alphaSpendingType = "OBF"){
  #number of (interim) analyses, total sample size, information quotient
  NoA <- length(Ns)
  N <- Ns[NoA]
  taus <- Ns/N
  
  # total alpha level at look-up (index corresponds to number of the look-up)
  alphaStar <- alphaSpendingFunction(taus, alpha, type = alphaSpendingType)
  betaStar  <- alphaSpendingFunction(taus, beta,  type = alphaSpendingType)
  
  # This vector will contain the alpha values at the interim analyses
  alphas <- matrix(nrow = 1, ncol = NoA)
  cumulativeAlpha <- matrix(nrow = 1, ncol = NoA)
  betas <- matrix(nrow = 1, ncol = NoA)
  cumulativeBeta  <- matrix(nrow = 1, ncol = NoA)
  cumulativePower <- matrix(nrow = 1, ncol = NoA)
  
  
  #data simulation
  x <- rnorm(n = N*NrOfSimulations, mean = 0)
  x2<- rnorm(n = N*NrOfSimulations, mean = 0)
  y <- rnorm(n = N*NrOfSimulations, mean = cohensD)
  
  #more convenient format
  x <- matrix(x , nrow = NrOfSimulations, ncol = N)
  x2<- matrix(x2, nrow = NrOfSimulations, ncol = N)
  y <- matrix(y , nrow = NrOfSimulations, ncol = N)
  
  # Vector which saves until which interim analysis you went
  stopAtHA <- rep(1, NrOfSimulations)
  stopAtH0 <- rep(1, NrOfSimulations)
  
  #index which will track which simulations haven't reached significance or rejection yet
  undecidedH0 <-  undecidedHA <-  undecidedTotal <- TRUE
  
  sseGroup1 <- sseGroup2H0 <- sseGroup2HA <- TH0 <- THA <- pH0 <- pHA <- matrix(nrow =1, ncol = NrOfSimulations)
  
  for(index in 1:NoA){
    if(index == 1){
      cN <- Ns[index]
      
      if(cN==1){
        groupSum1  <-  x[,1]
        groupSumH0 <- x2[,1]
        groupSumHA <-  y[,1]
      }else{
        groupSum1  <- apply( x[,1:cN],1,sum)
        groupSumH0 <- apply(x2[,1:cN],1,sum)
        groupSumHA <- apply( y[,1:cN],1,sum)
      }
      
      
    }else{
      cN <- Ns[index]
      oN <- Ns[index -1]
      
      if((cN-oN)==1){
        
        groupSum1[undecidedTotal] <- groupSum1[undecidedTotal] + x[undecidedTotal,(oN+1):cN]
        groupSumH0[  undecidedH0] <- groupSumH0[  undecidedH0] + x2[  undecidedH0,(oN+1):cN]
        groupSumHA[  undecidedHA] <- groupSumHA[  undecidedHA] +  y[  undecidedHA,(oN+1):cN]
      }else{
        
        groupSum1[undecidedTotal] <- groupSum1[undecidedTotal] + apply(x[undecidedTotal,(oN+1):cN],1,sum)
        groupSumH0[  undecidedH0] <- groupSumH0[  undecidedH0] + apply(x2[  undecidedH0,(oN+1):cN],1,sum)
        groupSumHA[  undecidedHA] <- groupSumHA[  undecidedHA] + apply( y[  undecidedHA,(oN+1):cN],1,sum)
      }
      
      # This is not an obvious place to put it, but it prevents the last analysis being count (as you cannot continu after the last analysis)
      stopAtH0 <- stopAtH0 + undecidedH0
      stopAtHA <- stopAtHA + undecidedHA
    }
    
    sseGroup1[undecidedTotal] <- apply(( x[undecidedTotal,1:cN]-groupSum1[undecidedTotal]/cN)^2, 1, sum)
    sseGroup2H0[undecidedH0]  <- apply((x2[undecidedH0,   1:cN]-groupSumH0[  undecidedH0]/cN)^2, 1, sum)
    sseGroup2HA[undecidedHA]  <- apply(( y[undecidedHA,   1:cN]-groupSumHA[  undecidedHA]/cN)^2, 1, sum)
    
    TH0[undecidedH0] <- sqrt((cN-1)/cN)*(groupSumH0[undecidedH0]-groupSum1[undecidedH0])/(sseGroup2H0[undecidedH0]+sseGroup1[undecidedH0])^0.5
    THA[undecidedHA] <- sqrt((cN-1)/cN)*(groupSumHA[undecidedHA]-groupSum1[undecidedHA])/(sseGroup2HA[undecidedHA]+sseGroup1[undecidedHA])^0.5
    
    if(sides == 2){
      TH0[undecidedH0] = abs(TH0[undecidedH0])
      THA[undecidedHA] = abs(THA[undecidedHA])
    }
    
    pH0[undecidedH0] <- pt(TH0[undecidedH0], df = 2*cN-2, lower.tail = FALSE)
    pHA[undecidedHA] <- pt(THA[undecidedHA], df = 2*cN-2, lower.tail = FALSE)
    
    
    #Determine the alpha at this interim analysis, such that our cumulative alpha matches the
    #desired cumulative alpha determined by the alpha spending function
    alphas[index] <- pt(quantile(TH0, probs = 1-alphaStar[index]), df = 2*cN-2, lower.tail = FALSE)
    betas[index]  <- pt(quantile(THA, probs =   betaStar[index]),  df = 2*cN-2, lower.tail = FALSE)
    
    if(alphas[index]>betas[index]){betas[index] <- alphas[index]}
    
    sigH0 <- pH0 < alphas[index]
    sigHA <- pHA < alphas[index]
    
    crossUpH0 <- pH0 > betas[index]
    crossUpHA <- pHA > betas[index]
    
    
    #previously significant requires you are now significant
    TH0[sigH0] <- Inf
    THA[sigHA] <- Inf
    
    pH0[sigH0] <- 0
    pHA[sigHA] <- 0
    
    TH0[crossUpH0] <- -Inf
    THA[crossUpHA] <- -Inf
    
    pH0[crossUpH0] <- 1
    pHA[crossUpHA] <- 1
    
    cumulativeAlpha[index] <- sum(sigH0)/NrOfSimulations
    cumulativePower[index] <- sum(sigHA)/NrOfSimulations
    cumulativeBeta[index]  <- sum(crossUpHA)/NrOfSimulations
    
    undecidedH0 <- !sigH0 & !crossUpH0
    undecidedHA <- !sigHA & !crossUpHA
    undecidedTotal <- undecidedH0 | undecidedHA
  }
  
  cumulativeBeta[NoA] <- 1-cumulativePower[NoA]
  
  #expected Stop and Number of Mice
  eStopHA <- sum(stopAtHA)/NrOfSimulations
  eNoMha <- sum(2*Ns[stopAtHA])/NrOfSimulations
  eStopH0 <- sum(stopAtH0)/NrOfSimulations
  eNoMh0 <- sum(2*Ns[stopAtH0])/NrOfSimulations
  
  #order of magnitude of the accuracy
  aOoM <- ceiling(log10(NrOfSimulations))-1
  
  results <- list(cumulativePower = round(cumulativePower, digits=aOoM), alphas = round(alphas*sides, digits=aOoM), betas = round(betas*sides, digits=aOoM),
                  expectedNumberofMiceHA = round(eNoMha, digits=aOoM), expectedStopHA = round(eStopHA, digits=aOoM) ,
                  expectedNumberofMiceH0 = round(eNoMh0, digits=aOoM), expectedStopH0 = round(eStopH0, digits=aOoM) ,
                  cumulativeAlpha = round(cumulativeAlpha, digits=aOoM), targetCumulativeAlpha = alphaStar,
                  cumulativeBeta = round(cumulativeBeta, digits=aOoM), targetCumulativeBeta = betaStar)
  return(results)
}

simIntAnFtestAB <- function(mu, sigma, Ns, alpha, beta, NrOfSimulations, alphaSpendingType = "OBF"){
  #number of (interim) analyses, total sample size, information quotient, number of groups
  NoA <- length(Ns)
  NoG <- length(mu)
  totalN <- Ns*NoG
  N <- totalN[NoA]
  taus <- totalN/N
  
  # total alpha level at look-up (index corresponds to number of the look-up)
  alphaStar <- alphaSpendingFunction(taus, alpha, type = alphaSpendingType)
  betaStar  <- alphaSpendingFunction(taus, beta,  type = alphaSpendingType)
  
  # This vector will contain the alpha values at the interim analyses
  alphas <- matrix(nrow = 1, ncol = NoA)
  cumulativeAlpha <- matrix(nrow = 1, ncol = NoA)
  betas <- matrix(nrow = 1, ncol = NoA)
  cumulativeBeta  <- matrix(nrow = 1, ncol = NoA)
  cumulativePower <- matrix(nrow = 1, ncol = NoA)
  
  #data simulation (no rescaling for df yet)
  # Denominator
  H0 <- rnorm(n = N*NrOfSimulations, mean = 0)
  # Numerator H0
  HA <- rnorm(n = N*NrOfSimulations, mean = mu, sd = sigma)
  
  #More convenient format
  H0 <- array(H0, dim = c(NoG, NrOfSimulations, N/NoG))
  HA <- array(HA, dim = c(NoG, NrOfSimulations, N/NoG))

  # Vector which saves until which interim analysis you went
  stopAtHA <- rep(1, NrOfSimulations)
  stopAtH0 <- rep(1, NrOfSimulations)
  
  #index which will track which simulations haven't reached significance or rejection yet
  undecidedH0 <-  undecidedHA <- TRUE
  
  totalSumH0 <- totalSumHA <- betweenGroupH0 <- betweenGroupHA <- matrix(nrow =1, ncol = NrOfSimulations)
  withinGroupH0 <- withinGroupHA <- FH0 <- FHA <- pH0 <- pHA <- matrix(nrow =1, ncol = NrOfSimulations)
  
  for(index in 1:NoA){
    if(index == 1){
      cN <- Ns[index]
      
      groupSumH0 <- apply(H0[,,1:cN],c(1,2),sum)
      groupSumHA <- apply(HA[,,1:cN],c(1,2),sum)
      
      CNoSH0 <- CNoSHA <- NrOfSimulations
      
    }else{
      cN <- Ns[index]
      oN <- Ns[index -1]
      
      groupSumH0[,undecidedH0] <- groupSumH0[,undecidedH0] + apply(H0[,undecidedH0,(oN+1):cN],c(1,2),sum)
      groupSumHA[,undecidedHA] <- groupSumHA[,undecidedHA] + apply(HA[,undecidedHA,(oN+1):cN],c(1,2),sum)
      
      CNoSH0 <- sum(undecidedH0)
      CNoSHA <- sum(undecidedHA)
      
      # This is not an obvious place to put it, but it prevents the last analysis being count (as you cannot continu after the last analysis)
      stopAtH0 <- stopAtH0 + undecidedH0
      stopAtHA <- stopAtHA + undecidedHA
    }
    
    totalSumH0[undecidedH0] <- apply(groupSumH0[,undecidedH0], 2, sum)
    totalSumHA[undecidedHA] <- apply(groupSumHA[,undecidedHA], 2, sum)
    
    betweenGroupH0[undecidedH0] <- apply(cN/(NoG-1)*(groupSumH0[,undecidedH0]/cN- matrix(data = totalSumH0[undecidedH0]/(cN*NoG), 
                                                                                   nrow = NoG, ncol = CNoSH0, byrow = TRUE))^2, 2, sum)
    
    betweenGroupHA[undecidedHA] <- apply(cN/(NoG-1)*(groupSumHA[,undecidedHA]/cN- matrix(data = totalSumHA[undecidedHA]/(cN*NoG), 
                                                                                   nrow = NoG, ncol = CNoSHA, byrow = TRUE))^2, 2, sum)
    
    withinGroupH0[undecidedH0] <- apply(1/((cN-1)*NoG)*(H0[,undecidedH0,1:cN]-array(groupSumH0[,undecidedH0]/cN, 
                                                                              dim = c(NoG,CNoSH0,cN)))^2, 2, sum)
    withinGroupHA[undecidedHA] <- apply(1/((cN-1)*NoG)*(HA[,undecidedHA,1:cN]-array(groupSumHA[,undecidedHA]/cN, 
                                                                              dim = c(NoG,CNoSHA,cN)))^2, 2, sum)
    
    FH0[undecidedH0] <- betweenGroupH0[undecidedH0]/withinGroupH0[undecidedH0]
    FHA[undecidedHA] <- betweenGroupHA[undecidedHA]/withinGroupHA[undecidedHA]
    
    pH0[undecidedH0] <- pf(FH0[undecidedH0], df1 = NoG-1, df2 = totalN[index]-NoG, lower.tail = FALSE)
    pHA[undecidedHA] <- pf(FHA[undecidedHA], df1 = NoG-1, df2 = totalN[index]-NoG, lower.tail = FALSE)
    
    
    #Determine the alpha at this interim analysis, such that our cumulative alpha matches the
    #desired cumulative alpha determined by the alpha spending function
    alphas[index] <- pf(quantile(FH0, probs = 1-alphaStar[index]), df1 = NoG-1,
                        df2 = totalN[index]-NoG, lower.tail = FALSE)
    betas[index]  <- pf(quantile(FHA, probs = betaStar[index]), df1 = NoG-1,
                        df2 = totalN[index]-NoG, lower.tail = FALSE)
    
    if(alphas[index]>betas[index]){betas[index] <- alphas[index]}
    
    sigH0 <- pH0 < alphas[index]
    sigHA <- pHA < alphas[index]
    
    crossUpH0 <- pH0 > betas[index]
    crossUpHA <- pHA > betas[index]
    
    
    #previously significant requires you are now significant
    FH0[sigH0] <- Inf
    FHA[sigHA] <- Inf
    
    pH0[sigH0] <- 0
    pHA[sigHA] <- 0
    
    FH0[crossUpH0] <- -Inf
    FHA[crossUpHA] <- -Inf
    
    pH0[crossUpH0] <- 1
    pHA[crossUpHA] <- 1
    
    cumulativeAlpha[index] <- sum(sigH0)/NrOfSimulations
    cumulativePower[index] <- sum(sigHA)/NrOfSimulations
    cumulativeBeta[index]  <- sum(crossUpHA)/NrOfSimulations
    
    undecidedH0 <- !sigH0 & !crossUpH0
    undecidedHA <- !sigHA & !crossUpHA
  }
  
  cumulativeBeta[NoA] <- 1-cumulativePower[NoA]
  
  #expected Stop and Number of Mice
  eStopHA <- sum(stopAtHA)/NrOfSimulations
  eNoMha <- sum(2*Ns[stopAtHA])/NrOfSimulations
  eStopH0 <- sum(stopAtH0)/NrOfSimulations
  eNoMh0 <- sum(2*Ns[stopAtH0])/NrOfSimulations
  
  #order of magnitude of the accuracy
  aOoM <- ceiling(log10(NrOfSimulations))-1
  
  results <- list(cumulativePower = round(cumulativePower, digits=aOoM), alphas = round(alphas*sides, digits=aOoM), betas = round(betas*sides, digits=aOoM),
                  expectedNumberofMiceHA = round(eNoMha, digits=aOoM), expectedStopHA = round(eStopHA, digits=aOoM) ,
                  expectedNumberofMiceH0 = round(eNoMh0, digits=aOoM), expectedStopH0 = round(eStopH0, digits=aOoM) ,
                  cumulativeAlpha = round(cumulativeAlpha, digits=aOoM), targetCumulativeAlpha = alphaStar,
                  cumulativeBeta = round(cumulativeBeta, digits=aOoM), targetCumulativeBeta = betaStar)
  
  return(results)
}
