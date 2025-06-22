#load packages
library(tidyverse)

#Pedigree derived abundance estimate 
## Creel-Rosenblatt
# 2022
# summary of 10 pedigrees
results_22 <- data.frame(iteration = seq(1,10), 
                         N.match.alive = c(23,24,24,24,22,24,22,24,22,22), 
                         N.inferred.alive = c(15,17,16,18,15,15,16,15,15,14), 
                         N.Sampled.alive = rep(57, 10))

# function
nHat_est_coyote_22 <- function(dataset, iter) {
  require(tidyverse)
  
  #First, we create storage for posterior estimates of probability of detection (pdet.para.blank) and probability of matching individuals (pmatch.para.blank), using sample data from the population-based approach. These storage objects store beta distribution shape parameters alpha and beta.
  pdet.para.blank <- list(data.frame("pedigree"  =  NA,  
                                     "Year 1"  =  NA),  
                          data.frame("pedigree"  =  NA,  
                                     "Year 1"  =  NA))
  
  names(pdet.para.blank) <- c("alpha", "beta") #Name the detection shape parameter storage
  
  pmatch.para.blank <- list(data.frame("pedigree"  =  NA, 
                                       "Year 1"  =  NA),  
                            data.frame("pedigree"  =  NA, 
                                       "Year 1"  =  NA))
  
  names(pmatch.para.blank) <- c("alpha", "beta")
  
  #Create storage for beta distribution parameters - these lists will include copies of pdet.para and pmatch.para, each corresponding to a simulation.
  pdet.list <- list(NA)
  pmatch.list <- list(NA)
  
  uninformative.prior <- 1 #Set uninformative prior for all parameter estimation (results in a uniform beta distribution)
  
  for(i in 1:iter){ #For every simulation...
    
    pedigree <- sample(seq(1,nrow(dataset)), 1, replace = TRUE)
    
    pdet.list[[i]] <- pdet.para.blank #Create a blank copy of probability of detection parameters in running list...
    names(pdet.list)[i] <- i #... and name the copy the simulation number
    
    pmatch.list[[i]] <- pmatch.para.blank #Create a blank copy of probability of matching parameters in running list...
    names(pmatch.list)[i] <- i #... and name the copy the simulation number
    
    temp <- dataset[pedigree, ] #Isolate the summary values for this particular simulation
    pdet.list[[i]]$alpha[1] <- pedigree
    pdet.list[[i]]$alpha[2] <- uninformative.prior+temp$N.match.alive #Calculate pdetection alpha shape parameter
    
    pdet.list[[i]]$beta[1] <- pedigree
    pdet.list[[i]]$beta[2] <- uninformative.prior+round(runif(n=1, min = 2, max = temp$N.inferred.alive),0) #Calculate pdetection beta shape parameter, drawing from a uniform distribution with min and max number of possible inferred individuals
    
    pmatch.list[[i]]$alpha[1] <- pedigree
    pmatch.list[[i]]$alpha[2] <- uninformative.prior+temp$N.match.alive  #Calculate pmatch alpha shape parameter
    
    pmatch.list[[i]]$beta[1] <- pedigree
    pmatch.list[[i]]$beta[2] <- (uninformative.prior+temp$N.Sampled.alive - temp$N.match.alive) #Calculate pmatch beta shape parameter
  } #close iteration loop
  
  pop_est <- data.frame(det = NULL, match = NULL, N.match.alive = NULL, N.inferred = NULL)
  
  for(i in 1:length(pdet.list)){ #Now, for each iteration
    
    #...Isolate shape parameters for p detection and p match
    det.alpha <- pdet.list[[i]]$alpha[,2]
    det.beta <- pdet.list[[i]]$beta[,2]
    match.alpha <- pmatch.list[[i]]$alpha[,2]
    match.beta <- pmatch.list[[i]]$beta[,2]
    
    #Using these shape parameters, randomly sample from beta distributions to estimate detection and matching probabilities
    det <- rbeta(n  =  1000, shape1 = det.alpha,  shape2  =  det.beta)
    match  <-rbeta(n  =  1000, shape1 = match.alpha,  shape2   =  match.beta)
    
    pop_est_i <- data.frame(det = det, match = match, N.sampled= dataset$N.Sampled.alive[pdet.list[[i]]$alpha[,1]], N.matched = dataset$N.match.alive[pdet.list[[i]]$alpha[,1]], N.inferred = pdet.list[[i]]$beta[,2])
    pop_est <- rbind(pop_est, pop_est_i)
    #Calculate adult abundance by dividing the number of matched individuals by the joint probability of being detected(sampled) and matched with another individual.
  } #End det and match prob draws
  
  #estimate pop size
  pop_est %>% 
    mutate(., nhat.1 = N.sampled/(det), #ignores information from matched/unmatched individuals
           nhat.2 = N.matched/(det*match), #Because of the random draws from two probabilities, this can result in Nhat < nsampled...
           nhat.3 = N.sampled + N.inferred + ((N.sampled-N.matched)/det), #Again, ignores information about pMatch
           nhat.4 = (N.sampled + N.inferred)/((det*match)+(det*(1-match))+((1-det)*match)))-> pop_est #Go with this one, incorporates both probabilities and ensures that nhat is greater than nsampled+ninferred
  
  print(paste0("Average population size: ", round(mean(pop_est$nhat.4),0)))
  print(paste0("95% confidence interval: ", round(quantile(pop_est$nhat.4, probs = 0.025),0), " - ", round(quantile(pop_est$nhat.4, probs = 0.975),0)))
} #End function

# Estimate
nHat_est_coyote_22(results_22, iter = 1000)

# 2023
# summary of 10 pedigrees 
results_23 <- data.frame(iteration = seq(1,10), 
                         N.match.alive = c(25,25,25,27,27,27,26,27,28,25), 
                         N.inferred.alive = c(21,21,21,21,21,21,21,22,20,21), 
                         N.Sampled.alive = rep(55, 10))
# function
nHat_est_coyote_23 <- function(dataset, iter) {
  require(tidyverse)
  
  #First, we create storage for posterior estimates of probability of detection (pdet.para.blank) and probability of matching individuals (pmatch.para.blank), using sample data from the population-based approach. These storage objects store beta distribution shape parameters alpha and beta.
  pdet.para.blank <- list(data.frame("pedigree"  =  NA,  
                                     "Year 1"  =  NA),  
                          data.frame("pedigree"  =  NA,  
                                     "Year 1"  =  NA))
  
  names(pdet.para.blank) <- c("alpha", "beta") #Name the detection shape parameter storage
  
  pmatch.para.blank <- list(data.frame("pedigree"  =  NA, 
                                       "Year 1"  =  NA),  
                            data.frame("pedigree"  =  NA, 
                                       "Year 1"  =  NA))
  
  names(pmatch.para.blank) <- c("alpha", "beta")
  
  #Create storage for beta distribution parameters - these lists will include copies of pdet.para and pmatch.para, each corresponding to a simulation.
  pdet.list <- list(NA)
  pmatch.list <- list(NA)
  
  uninformative.prior <- 1 #Set uninformative prior for all parameter estimation (results in a uniform beta distribution)
  
  for(i in 1:iter){ #For every simulation...
    
    pedigree <- sample(seq(1,nrow(dataset)), 1, replace = TRUE)
    
    pdet.list[[i]] <- pdet.para.blank #Create a blank copy of probability of detection parameters in running list...
    names(pdet.list)[i] <- i #... and name the copy the simulation number
    
    pmatch.list[[i]] <- pmatch.para.blank #Create a blank copy of probability of matching parameters in running list...
    names(pmatch.list)[i] <- i #... and name the copy the simulation number
    
    temp <- dataset[pedigree, ] #Isolate the summary values for this particular simulation
    pdet.list[[i]]$alpha[1] <- pedigree
    pdet.list[[i]]$alpha[2] <- uninformative.prior+temp$N.match.alive #Calculate pdetection alpha shape parameter
    
    pdet.list[[i]]$beta[1] <- pedigree
    pdet.list[[i]]$beta[2] <- uninformative.prior+round(runif(n=1, min = 2, max = temp$N.inferred.alive),0) #Calculate pdetection beta shape parameter, drawing from a uniform distribution with min and max number of possible inferred individuals
    
    pmatch.list[[i]]$alpha[1] <- pedigree
    pmatch.list[[i]]$alpha[2] <- uninformative.prior+temp$N.match.alive  #Calculate pmatch alpha shape parameter
    
    pmatch.list[[i]]$beta[1] <- pedigree
    pmatch.list[[i]]$beta[2] <- (uninformative.prior+temp$N.Sampled.alive - temp$N.match.alive) #Calculate pmatch beta shape parameter
  } #close iteration loop
  
  pop_est <- data.frame(det = NULL, match = NULL, N.match.alive = NULL, N.inferred = NULL)
  
  for(i in 1:length(pdet.list)){ #Now, for each iteration
    
    #...Isolate shape parameters for p detection and p match
    det.alpha <- pdet.list[[i]]$alpha[,2]
    det.beta <- pdet.list[[i]]$beta[,2]
    match.alpha <- pmatch.list[[i]]$alpha[,2]
    match.beta <- pmatch.list[[i]]$beta[,2]
    
    #Using these shape parameters, randomly sample from beta distributions to estimate detection and matching probabilities
    det <- rbeta(n  =  1000, shape1 = det.alpha,  shape2  =  det.beta)
    match  <-rbeta(n  =  1000, shape1 = match.alpha,  shape2   =  match.beta)
    
    pop_est_i <- data.frame(det = det, match = match, N.sampled= dataset$N.Sampled.alive[pdet.list[[i]]$alpha[,1]], N.matched = dataset$N.match.alive[pdet.list[[i]]$alpha[,1]], N.inferred = pdet.list[[i]]$beta[,2])
    pop_est <- rbind(pop_est, pop_est_i)
    #Calculate adult abundance by dividing the number of matched individuals by the joint probability of being detected(sampled) and matched with another individual.
  } #End det and match prob draws
  
  #estimate pop size
  pop_est %>% 
    mutate(., nhat.1 = N.sampled/(det), #ignores information from matched/unmatched individuals
           nhat.2 = N.matched/(det*match), #Because of the random draws from two probabilities, this can result in Nhwhat at < nsampled...
           nhat.3 = N.sampled + N.inferred + ((N.sampled-N.matched)/det), #Again, ignores information about pMatch
           nhat.4 = (N.sampled + N.inferred)/((det*match)+(det*(1-match))+((1-det)*match)))-> pop_est #Go with this one, incorporates both probabilities and ensures that nhat is greater than nsampled+ninferred
  
  print(paste0("Average population size: ", round(mean(pop_est$nhat.4),0)))
  print(paste0("95% confidence interval: ", round(quantile(pop_est$nhat.4, probs = 0.025),0), " - ", round(quantile(pop_est$nhat.4, probs = 0.975),0)))
} #End function

# Estimate
nHat_est_coyote_23(results_23, iter = 1000)
