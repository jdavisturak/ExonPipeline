## Function to simulate genes/model

# First just need to go from Dist2End to splicing
modelIntron = function(Dist2End=NA, k_splice, ElongRate=NA,Time=NA){
  if(is.na(Time[1])){
    Time = Dist2End/ElongRate
  }
  1-exp(-k_splice*Time)  
}

modelGeneEfficiency = function(Dists=NA,ks,Elongs=NA,Time=NA){
  prod(modelIntron(Dists,ks,Elongs,Time))
}


### limit the absolute value of a distribution to the MAX of the distro
limitToMax = function(X, myMax=max(X,na.rm=T)){
  X[X < -myMax] <- -myMax
  X
}


### Calculate time from Distance, Elong, and modified Elong over a short distance
timeFromModifiedElong = function(Dist,Elong, Elong2, Dist2){
  max(c(Dist-Dist2),0)/Elong + Dist2/Elong2
}


### Go from Acceptor Strength to k:
Acceptor2k = function(Acceptors, K=1, div_factor=3){
  Acceptors[is.na(Acceptors)] <- 0
  #K*(1 + Acceptors/div_factor)
  K*exp(Acceptors/div_factor)
}

### Go from Stability Score to Elong:
Stability2Elong = function(Stability, Elong=3000, div_factor=4, exponent=.1){
  Stability[is.na(Stability)] <- 0
  # Raise stability to an exponent (< 1 generally)
  polarity = ifelse(Stability<0,-1,1)
  Stability = polarity*(Stability*polarity)^exponent
  Elong*(1 - Stability/div_factor^exponent)
}

Stability2Elong2 = function(Stability, Elong=3000, MinOffset=5,sdlog=0.5){
   Y = dlnorm(Stability + MinOffset, sdlog=sdlog)
   Y*Elong/median(Y)
}
  

### Given an ORDERED series of Dist2End and modified Elong times, calculate the total Elong time for each segment
getTimes = function(Dists, ModElongs, Elong=3000, extraDist=rep(80,length(Dists))){
  
  Times = rep(NA,length(Dists))
  oldDist = 0
  for(i in length(Dists):1){
    Times[i] = timeFromModifiedElong(Dists[i]-oldDist, Elong, ModElongs[i],extraDist[i])
    oldDist = Dists[i]    
  }

  rev(cumsum(rev(Times)))
}


simulateGene = function(myData, Elong=3000, K_splice=1, Stability_div=4, Acceptor_div=3, extraDist=80, maxAcceptor=2.678306, readThrough=1000,Stability_exp=.25){
  # This should be rows from mergedData

  #ORDER by feature
  #myData=myData[order(myData$FeatureCount),]
  #b=proc.time(); 
   
  # Get all the Modified Elongation Times (Energy)
  modElongs = Stability2Elong(myData$Stability, Elong, Stability_div,exponent=Stability_exp)
  #c=proc.time(); elongTime <<- elongTime+c-b
  
  # Get all the total elongation times
  Times = getTimes(myData$Dist2End, modElongs, Elong,extraDist) + readThrough/Elong
  #d=proc.time(); timeTime <<- timeTime+d-c
  
  # Get all the K_splice (Acceptor)
  #Ks = Acceptor2k(limitToMax(myData$Acceptor,maxAcceptor), K_splice, Acceptor_div)
  Ks = rep(K_splice,nrow(myData))
  Ks[length(Ks)] = K_splice * (1 + 1/Acceptor_div)
  #e=proc.time(); kTime <<- kTime+e-d
  
  # Simulate!
  eff=modelGeneEfficiency(ks=Ks,Time=Times)    
  #f=proc.time(); simulateTime <<- simulateTime+f-e
    
  eff
 } 

 simulateGene2 = function(myData, Elong=3000, K_splice=1, Stability_MinOffset=5,Stability_sdlog=0.5, Acceptor_div=3, extraDist=80, maxAcceptor=2.678306, readThrough=1000){
  # Get all the Modified Elongation Times (Energy)
  modElongs = Stability2Elong2(myData$Stability, Elong, Stability_MinOffset,Stability_sdlog)
  
  # Get all the total elongation times
  Times = getTimes(myData$Dist2End, modElongs, Elong,extraDist) + readThrough/Elong
  
  # Get all the K_splice (Acceptor)
  #Ks = Acceptor2k(limitToMax(myData$Acceptor,maxAcceptor), K_splice, Acceptor_div)
  Ks = rep(K_splice,nrow(myData))
  Ks[length(Ks)] = K_splice * (1 + 1/Acceptor_div)
  
  # Simulate!
  eff=modelGeneEfficiency(ks=Ks,Time=Times)    
    
  eff
 } 
