alphaSpendingFunction <- function(tau, alpha, type = "OBF"){
  
  if(type == "OBF"){
    #alpha spending function approximating O'Brien-Fleming boundaries
    result <- 2 - 2*pnorm(qnorm(1-alpha/2)/sqrt(tau))
    
  }else if(type == "Pocock"){ 
    #alpha spending function approximating Pocock boundaries
    result <- alpha*log(1+(exp(1)-1)^tau)
    
  }else if(type == "compromis"){
    result <- alpha*tau
  }
  
  return(result)
}

simIntAtTest <- function(cohensD, Ns, alpha, NrOfSimulations, sides = 1, alphaSpendingType = "OBF"){
  #number of (interim) analyses, total sample size, information quotient
  NoA <- length(Ns)
  N <- Ns[NoA]
  taus <- Ns/N
  
  # total alpha level at look-up (index corresponds to number of the look-up)
  alphaStar <- alphaSpendingFunction(taus, alpha, type = alphaSpendingType)
  
  # This vector will contain the alpha values at the interim analyses
  alphas <- matrix(nrow = 1, ncol = NoA)
  cumulativeAlpha <- matrix(nrow = 1, ncol = NoA)
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
  stopAt <- rep(1, NrOfSimulations)
  
  #index which will track which simulations haven't reached significance yet
  notSigH0 <- notSigHA <- notSigTotal <- TRUE
  
  sseGroup1 <- sseGroup2H0 <- sseGroup2HA <- TH0 <- THA <- pH0 <- pHA <- matrix(nrow =1, ncol = NrOfSimulations)
  
  for(index in 1:NoA){
    if(index == 1){
      cN <- Ns[index]
      
      groupSum1  <- apply( x[,1:cN],1,sum)
      groupSumH0 <- apply(x2[,1:cN],1,sum)
      groupSumHA <- apply( y[,1:cN],1,sum)
      
    }else{
      cN <- Ns[index]
      oN <- Ns[index -1]
      
      groupSum1[notSigTotal] <- groupSum1[notSigTotal] + apply(x[notSigTotal,(oN+1):cN],1,sum)
      groupSumH0[  notSigH0] <- groupSumH0[  notSigH0] + apply(x2[  notSigH0,(oN+1):cN],1,sum)
      groupSumHA[  notSigHA] <- groupSumHA[  notSigHA] + apply( y[  notSigHA,(oN+1):cN],1,sum)
      
      # This is not an obvious place to put it, but it prevents the last analysis being count (as you cannot continu after the last analysis)
      stopAt <- stopAt + notSigHA
    }
    
    sseGroup1[notSigTotal] <- apply(( x[notSigTotal,1:cN]-groupSum1[notSigTotal]/cN)^2, 1, sum)
    sseGroup2H0[notSigH0]  <- apply((x2[notSigH0,   1:cN]-groupSumH0[  notSigH0]/cN)^2, 1, sum)
    sseGroup2HA[notSigHA]  <- apply(( y[notSigHA,   1:cN]-groupSumHA[  notSigHA]/cN)^2, 1, sum)
    
    TH0[notSigH0] <- sqrt((cN-1)/cN)*(groupSumH0[notSigH0]-groupSum1[notSigH0])/(sseGroup2H0[notSigH0]+sseGroup1[notSigH0])^0.5
    THA[notSigHA] <- sqrt((cN-1)/cN)*(groupSumHA[notSigHA]-groupSum1[notSigHA])/(sseGroup2HA[notSigHA]+sseGroup1[notSigHA])^0.5
    
    if(sides == 2){
      TH0[notSigH0] = abs(TH0[notSigH0])
      THA[notSigHA] = abs(THA[notSigHA])
    }
    
    pH0[notSigH0] <- pt(TH0[notSigH0], df = 2*cN-2, lower.tail = FALSE)
    pHA[notSigHA] <- pt(THA[notSigHA], df = 2*cN-2, lower.tail = FALSE)
    
    
    #Determine the alpha at this interim analysis, such that our cumulative alpha matches the
    #desired cumulative alpha determined by the alpha spending function
    alphas[index] <- pt(quantile(TH0, probs = 1-alphaStar[index]), df = 2*cN-2, lower.tail = FALSE)
    
    notSigH0 <- pH0 > alphas[index]
    notSigHA <- pHA > alphas[index]
    
    #previously significant requires you are now significant
    TH0[!notSigH0] <- Inf
    THA[!notSigHA] <- Inf
    
    pH0[!notSigH0] <- 0
    pHA[!notSigHA] <- 0
    
    cumulativeAlpha[index] <- 1-sum(notSigH0)/NrOfSimulations
    cumulativePower[index] <- 1-sum(notSigHA)/NrOfSimulations
  }
  
  #expected Stop and Number of Mice
  eStop <- sum(stopAt)/NrOfSimulations
  eNoM <- sum(2*Ns[stopAt])/NrOfSimulations
  
  #order of magnitude of the accuracy
  aOoM <- ceiling(log10(NrOfSimulations))-1
  
  results <- list(cumulativePower = round(cumulativePower, digits=aOoM), alphas = round(alphas*sides, digits=aOoM), 
                  expectedNumberofMice = round(eNoM, digits=aOoM), expectedStop = round(eStop, digits=aOoM) ,
                  cumulativeAlpha = round(cumulativeAlpha, digits=aOoM), 
                  targetCumulativeAlpha = round(alphaStar, digits=aOoM))
  return(results)
}

simIntAnFtest <- function(mu, sigma, Ns, alpha, NrOfSimulations, alphaSpendingType = "OBF"){
  #number of (interim) analyses, total sample size, information quotient, number of groups
  NoA <- length(Ns)
  NoG <- length(mu)
  totalN <- Ns*NoG
  N <- totalN[NoA]
  taus <- totalN/N
  
  # total alpha level at look-up (index corresponds to number of the look-up)
  alphaStar <- alphaSpendingFunction(taus, alpha, type = alphaSpendingType)
  
  #data simulation (no rescaling for df yet)
  # Denominator
  H0 <- rnorm(n = N*NrOfSimulations, mean = 0)
  # Numerator H0
  HA <- rnorm(n = N*NrOfSimulations, mean = mu, sd = sigma)
  
  #More convenient format
  H0 <- array(H0, dim = c(NoG, NrOfSimulations, N/NoG))
  HA <- array(HA, dim = c(NoG, NrOfSimulations, N/NoG))
  
  # This vector will contain the alpha values at the interim analyses
  alphas <- matrix(nrow = 1, ncol = NoA)
  cumulativeAlpha <- matrix(nrow = 1, ncol = NoA)
  cumulativePower <- matrix(nrow = 1, ncol = NoA)
  
  # Vector which saves until which interim analysis you went
  stopAt <- rep(1, NrOfSimulations)
  
  #index which will track which simulations haven't reached significance yet
  notSigH0 <- notSigHA <- TRUE
  
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
      
      groupSumH0[,notSigH0] <- groupSumH0[,notSigH0] + apply(H0[,notSigH0,(oN+1):cN],c(1,2),sum)
      groupSumHA[,notSigHA] <- groupSumHA[,notSigHA] + apply(HA[,notSigHA,(oN+1):cN],c(1,2),sum)
      
      CNoSH0 <- sum(notSigH0)
      CNoSHA <- sum(notSigHA)
      
      # This is not an obvious place to put it, but it prevents the last analysis being count (as you cannot continu after the last analysis)
      stopAt <- stopAt + notSigHA
    }
    
    totalSumH0[notSigH0] <- apply(groupSumH0[,notSigH0], 2, sum)
    totalSumHA[notSigHA] <- apply(groupSumHA[,notSigHA], 2, sum)
    
    betweenGroupH0[notSigH0] <- apply(cN/(NoG-1)*(groupSumH0[,notSigH0]/cN- matrix(data = totalSumH0[notSigH0]/(cN*NoG), 
                                                                                   nrow = NoG, ncol = CNoSH0, byrow = TRUE))^2, 2, sum)
    
    betweenGroupHA[notSigHA] <- apply(cN/(NoG-1)*(groupSumHA[,notSigHA]/cN- matrix(data = totalSumHA[notSigHA]/(cN*NoG), 
                                                                                   nrow = NoG, ncol = CNoSHA, byrow = TRUE))^2, 2, sum)
    
    withinGroupH0[notSigH0] <- apply(1/((cN-1)*NoG)*(H0[,notSigH0,1:cN]-array(groupSumH0[,notSigH0]/cN, 
                                                                              dim = c(NoG,CNoSH0,cN)))^2, 2, sum)
    withinGroupHA[notSigHA] <- apply(1/((cN-1)*NoG)*(HA[,notSigHA,1:cN]-array(groupSumHA[,notSigHA]/cN, 
                                                                              dim = c(NoG,CNoSHA,cN)))^2, 2, sum)
    
    FH0[notSigH0] <- betweenGroupH0[notSigH0]/withinGroupH0[notSigH0]
    FHA[notSigHA] <- betweenGroupHA[notSigHA]/withinGroupHA[notSigHA]
    
    pH0[notSigH0] <- pf(FH0[notSigH0], df1 = NoG-1, df2 = totalN[index]-NoG, lower.tail = FALSE)
    pHA[notSigHA] <- pf(FHA[notSigHA], df1 = NoG-1, df2 = totalN[index]-NoG, lower.tail = FALSE)
    
    #Determine the alpha at this interim analysis, such that our cumulative alpha matches the
    #desired cumulative alpha determined by the alpha spending function
    alphas[index] <- pf(quantile(FH0, probs = 1-alphaStar[index]), df1 = NoG-1,
                        df2 = totalN[index]-NoG, lower.tail = FALSE)
    
    notSigH0 <- pH0 > alphas[index]
    notSigHA <- pHA > alphas[index]
    
    #previously significant requires you are now significant
    FH0[!notSigH0] <- Inf
    FHA[!notSigHA] <- Inf
    
    pH0[!notSigH0] <- 0
    pHA[!notSigHA] <- 0
    
    cumulativeAlpha[index] <- 1-sum(notSigH0)/NrOfSimulations
    cumulativePower[index] <- 1-sum(notSigHA)/NrOfSimulations
  }
  
  
  #expected Stop and Number of Mice
  eStop <- sum(stopAt)/NrOfSimulations
  eNoM <- sum(totalN[stopAt])/NrOfSimulations
  
  #order of magnitude of the accuracy
  aOoM <- ceiling(log10(NrOfSimulations))-1
  
  results <- list(cumulativePower = round(cumulativePower, digits=aOoM), alphas = round(alphas, digits=aOoM), 
                  expectedNumberofMice = round(eNoM, digits=aOoM), expectedStop = round(eStop, digits=aOoM) ,
                  cumulativeAlpha = round(cumulativeAlpha, digits=aOoM), 
                  targetCumulativeAlpha = round(alphaStar, digits=aOoM))
  return(results)
}
