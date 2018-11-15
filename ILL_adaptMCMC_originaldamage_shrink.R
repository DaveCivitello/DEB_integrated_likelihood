### General adaptive MCMC ###
library("adaptMCMC")
library("deSolve")

rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

setwd("C:/RData")

# compile my model from C definition
dyn.unload("IndividualModel_Shrink_P.dll") # unload dll
system("R CMD SHLIB IndividualModel_Shrink_P.c")
dyn.load("IndividualModel_Shrink_P.dll") # Load dll


# sim.data = rnorm(1000, mean=4, sd=1)
# mean(sim.data)
# sd(sim.data)
# 
# density_func = function(x){
#   sum(dnorm(sim.data, mean=x[1], sd=x[2], log=T))
# }
# 
# test = MCMC(density_func, init=c("mean"=2, "sd"=2), scale=c(0.1, 0.1), adapt=T, acc.rate = 0.35, n=1000)
# 
# testc = convert.to.coda(test)
# plot(testc)



#### Fixed information ####

# 1 - Data
setwd("C:/RData")

data = read.csv("AllSnailsLTSurv.csv")
data2 = read.csv("schisto_periodic_starve_LT2.csv")

data = list(t = data$Date, L = data$Length, Negg = data$C_Eggs, Nworms = data$C_Worms, Alive=data$Last_Alive,
            L2 = data2$Length, Negg2 = data2$C_Eggs, Nworms2 = data2$C_Worms, Alive2=data2$Last_Alive)


# 2 - Initial conditions
setinits.Food<-function(F0 = 16.5, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Starve<-function(F0 = 10.74, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}




# 3 - Functions for feeding events
periodic.starvation.events = function(initial.food, starve.period, infection.date, weeks.duration){
  infection.value = 2.85e-5
  all.fed.dates = infection.date + c(1,4) # Snails were infected on a Thursday, All fed on Friday and Monday, then treatments
  potential.feedings = infection.date + 5 + sort(c((1:weeks.duration)*7 - 6,((1:weeks.duration)*7 - 4),((1:weeks.duration)*7 - 2)))
  #print(0:(weeks.duration/starve.period - 1)*6 + 1)
  potential.feedings.vals = rep(initial.food, times = length(potential.feedings))
  # work out starvation
  if( starve.period == 2){
    potential.feedings.vals[(potential.feedings - (infection.date + 5))  %% 14 <= 7] = 0
  }
  if (starve.period == 3){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 21 <= 14] = 0
  }
  if (starve.period == 4){
    potential.feedings.vals[(potential.feedings - (infection.date + 5)) %% 28 <= 21] = 0
  }
  event.dates = c(infection.date, all.fed.dates, potential.feedings, infection.date)
  event.values = c(infection.value, initial.food, initial.food, potential.feedings.vals, 0)
  methods = rep("replace", times=length(event.dates))
  vars = c("P", rep("F", times = length(event.dates)-2), "HAZ")
  data.frame(var=vars, time= event.dates, value = event.values, method= methods)
}

feeding.events <- function(dates, var="F", value, method="replace", Infected=0){
  # Assemble all data
  events = length(dates) # number of events
  vars = rep(var, times = events)
  values = rep(value, times = events)
  methods = rep(method, times = events)
  
  #build data.frame
  result = data.frame(var = vars, time = dates, value = values, method = methods)
  if(Infected > 0){
    result = rbind(result, data.frame(var="P", time=Infected, value=2.85e-5, method="replace"))
  }
  result = rbind(result, data.frame(var="HAZ", time=63, value=0, method="replace"))
  result
}


dur.F = 245
dur.S = 140

in.F = setinits.Food()
in.S = setinits.Starve()
Feed.F = feeding.events(dates = sort(c((1:35)*7,((1:35)*7 - 3))), var="F", in.F[1], method="replace", Infected=28)


# 4 - Functions to solve DEBs
solve.DEB.food<-function(params, inits, duration, feeding.events){
  feed.juv <- feeding.events
  feed.juv.U <- subset(feed.juv, var != "P")
  feed.juv.U[71, 2] = 28
  parms = as.numeric(params[1:24])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], kapc=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], kapnew=parms[15], yEF=parms[16], LM=parms[17],kd=parms[18], 
             z=parms[19], kk=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
  
  Juv.R <- as.data.frame(lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                               params[1:25],  rtol=1e-6, atol=1e-6,   
                               events = list(data = feed.juv)))
  
  Juv.RU <- as.data.frame(lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.juv.U)))
  
  feed.juv[1:70,3] <- 11
  feed.juv.U[1:70,3] <- 11
  Juv.G <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.juv))
  Juv.GU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 5.5
  feed.juv.U[1:70,3] <- 5.5
  Juv.U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.juv))  
  Juv.UU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 2.75
  feed.juv.U[1:70,3] <- 2.75
  Juv.P <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.juv))
  Juv.PU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 1.375
  feed.juv.U[1:70,3] <- 1.375
  Juv.BR <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv))
  Juv.BRU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                   initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                   params[1:25],  rtol=1e-6, atol=1e-6,   
                   events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 0.6875
  feed.juv.U[1:70,3] <- 0.6875
  Juv.BL <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv))
  Juv.BLU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                   initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                   params[1:25],  rtol=1e-6, atol=1e-6,   
                   events = list(data = feed.juv.U))
  # We need to overwrite the survival result for the day of conditioning, because lsoda records the
  # value before overwriting the state variable during the "event"
  if(length(Juv.R[,1]) >= 64 & length(Juv.RU[,1]) >= 29){Juv.R[64, "Survival"] = 1; Juv.RU[29, "Survival"] = 1}
  if(length(Juv.G[,1]) >= 64 & length(Juv.GU[,1]) >= 29){Juv.G[64, "Survival"] = 1; Juv.GU[29, "Survival"] = 1}
  if(length(Juv.U[,1]) >= 64 & length(Juv.UU[,1]) >= 29){Juv.U[64, "Survival"] = 1; Juv.UU[29, "Survival"] = 1}
  if(length(Juv.P[,1]) >= 64 & length(Juv.PU[,1]) >= 29){Juv.P[64, "Survival"] = 1; Juv.PU[29, "Survival"] = 1}
  if(length(Juv.BR[,1]) >= 64 & length(Juv.BRU[,1]) >= 29){Juv.BR[64, "Survival"] = 1; Juv.BRU[29, "Survival"] = 1}
  if(length(Juv.BL[,1]) >= 64 & length(Juv.BLU[,1]) >= 29){Juv.BL[64, "Survival"] = 1; Juv.BLU[29, "Survival"] = 1}
  
  result <- rbind(#Infecteds (n=66)
    Juv.BL, Juv.BL, Juv.BL, Juv.BL, Juv.BL, Juv.BL, Juv.BL, Juv.BL, Juv.BL, Juv.BL, 
    Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, Juv.BR, 
    Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, Juv.P, 
    Juv.U, Juv.U, Juv.U, Juv.U, Juv.U, Juv.U, Juv.U, Juv.U, Juv.U, Juv.U, 
    Juv.G, Juv.G, Juv.G, Juv.G, Juv.G, Juv.G, Juv.G, Juv.G, Juv.G, Juv.G, 
    Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, Juv.R, 
    #Uninfecteds (n=30)
    Juv.BLU, Juv.BLU, Juv.BLU, Juv.BLU, Juv.BLU, 
    Juv.BRU, Juv.BRU, Juv.BRU, Juv.BRU, Juv.BRU, 
    Juv.PU, Juv.PU, Juv.PU, Juv.PU, Juv.PU, 
    Juv.UU, Juv.UU, Juv.UU, Juv.UU, Juv.UU, 
    Juv.GU, Juv.GU, Juv.GU, Juv.GU, Juv.GU,
    Juv.RU, Juv.RU, Juv.RU, Juv.RU, Juv.RU)
  
  result
  
}

# This function is customized to my model and data
solve.DEB.starve<-function(params, inits, duration){
  # Collect the params the way C likes them
  parms = as.numeric(params[1:25])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], kapc=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], kapnew=parms[15], yEF=parms[25], LM=parms[17],kd=parms[18], 
             z=parms[19], kk=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=14)
  
  # Simulate dynamics for 1-0, set up feeding regime and fix conditional survival
  food_1_0 = periodic.starvation.events(2.69, 0, 14, 18)
  food_1_0.U <- subset(food_1_0, var != "P")
  food_1_0[which(food_1_0$var=="HAZ"),2] = 35
  
  out_1_0 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_1_0)))
  out_1_0U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_1_0.U)))
  
  # Simulate dynamics for 2-2
  food_2_2 = periodic.starvation.events(5.37, 2, 14, 18)
  food_2_2.U <- subset(food_2_2, var != "P")
  food_2_2[which(food_2_2$var=="HAZ"),2] = 35
  
  out_2_2 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_2_2)))
  out_2_2U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_2_2.U)))
  
  # Simulate dynamics for 3-3
  food_3_3 = periodic.starvation.events(8.06, 3, 14, 18)
  food_3_3.U <- subset(food_3_3, var != "P")
  food_3_3[which(food_3_3$var=="HAZ"),2] = 42
  
  out_3_3 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_3_3)))
  out_3_3U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_3_3.U)))
  
  # Simulate dynamics for 4-4
  food_4_4 = periodic.starvation.events(10.74, 4, 14, 18)
  food_4_4.U <- subset(food_4_4, var != "P")
  food_4_4[which(food_4_4$var=="HAZ"),2] = 49
  
  out_4_4 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_4_4)))
  out_4_4U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_4_4.U)))
  # No need to overwrite one day of survival for correct conditional probabilities, b/c nobody died in first conditional week
  result <- rbind(
    out_1_0, out_1_0, out_1_0, out_1_0, out_1_0, # 5
    out_2_2, out_2_2, out_2_2, # 3
    out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, out_3_3, # 8    
    out_4_4, out_4_4, out_4_4, out_4_4, # 4
    
    
    out_1_0U, out_1_0U, out_1_0U, out_1_0U, out_1_0U,
    out_2_2U, out_2_2U, out_2_2U, out_2_2U, out_2_2U,
    out_3_3U, out_3_3U, out_3_3U, out_3_3U, out_3_3U,
    out_4_4U, out_4_4U, out_4_4U, out_4_4U, out_4_4U)
  
  result
  
}

extract.data<-function(data, w.t=7, start.age=0){
  ww<- which(data$time%%w.t==0 & data$time>=start.age)
  data[ww,]
}

## This function takes the simulator, params, inits, etc, and runs a
## forward simulation, and then extracts a subset of the points from
## the simulator (determined by w.t), without noise
make.states<-function(params, inits.F, inits.S, duration.F, duration.S, feeding.events.F, w.t=7){
  result.F = solve.DEB.food(params, inits.F, duration.F, feeding.events.F)
  result.F = extract.data(result.F, start.age=28)
  Survival.F = result.F$Survival
  
  # This fix for survival data is currently specific to the sample size (96) and duration (32) of the experiment
  Survival.F[-((1:96)*32)] =   Survival.F[-((1:96)*32)] - Survival.F[-(1+(0:95)*32)]
  
  result.S = solve.DEB.starve(params, inits.S, duration.S)
  result.S = extract.data(result.S, start.age=14)
  Survival.S = result.S$Survival
  # This fix for survival data is currently specific to the sample size (40) and duration (19) of the experiment
  Survival.S[-((1:40)*19)] =   Survival.S[-((1:40)*19)] - Survival.S[-(1+(0:39)*19)]
  
  return(list(time=result.F$time, L=result.F$LG, RH=result.F$RH, RP=result.F$RP,
              L2=result.S$LG, W2=result.S$RP, E2=result.S$RH, SurvF=Survival.F, SurvS=Survival.S))
}

# 5 - Prior likelihood
prior.likelihood = function(x){
  prior.lik = with(as.list(x),
                   sum(dbeta(c(yPE, yEF, yEF2, yRP, yVE, mP, eh, k, fe), 1, 1, log=T)) + 
                     sum(dunif(c(sd.LI1, sd.LU1, sd.EI1, sd.EU1, sd.W1, sd.LI2, sd.LU2, sd.EI2, sd.EU2, sd.W2), min=0, max=10, log=T)) +
                     sum(dunif(c(ph, alpha, iPM, EM, DR, Fh, muD, kd, z, kk, theta, mR, hb), min=0, max=1000000, log=T)) +
                     dnorm(iM, mean=0.0183, sd=0.0016, log=T) + dnorm(M, mean=0.004, sd=0.00047, log=T) + dnorm(LM, mean=30, sd=1, log=T)
  )
  return(prior.lik)
}


integrated.likelihood = function(x){
  n=10
  # observation model
  gammaH<-0.015 # C content of eggs
  gammaP<-4e-5 # C content of cercs
  
  sd.L <- x["sd.L"] # on average w/in 0.2 mm measurement error
  sd.E <- x["sd.E"]# on average, counting w/in 10% eggs
  sd.W <- x["sd.W"] # on average, counting w/in 10% worms
  
  # random effect in parameter
  sd.LM = x["sd.LM"]
  LM.mean = x["LM"]
  qs = seq(from= 0.001, to = 0.999, length.out = n) #select quantiles
  LMs = qnorm(qs, mean=LM.mean, sd=sd.LM) # obtain random effect values
  ds = dnorm(LMs, mean=LM.mean, sd=sd.LM) # get the density of these values
  probs = ds/sum(ds)
  LLs = numeric()
  sim.data = list()
  for(i in 1:length(LMs)){
    x["LM"] = LMs[i]
    sim.data[[i]] = make.states(x, in.F, in.S, dur.F, dur.S, Feed.F, w.t=7)
  }
  NLL = 0
  for(snail in 1:96){ # subsetting by snail
    LL.ind = numeric()
    rows = 1:32 + (snail-1)*32
    snail.L = data$L[rows]
    snail.E = data$Negg[rows]
    snail.W = data$Nworms[rows]
    snail.death = which(data$Alive == 1)
    for(i in 1:length(n)){ # scrolling over different values of random effect parameter
      if(length(sim.data[[i]]$L[rows]) != length(snail.L) | anyNA(sim.data[[i]]$RH[rows])){
        LL.ind[i] = -1e6
        next
      }
      LL.ind[i] = sum(dnorm(x=log(snail.L), mean=log(sim.data[[i]]$L[rows]), sd=sd.L, log=T), na.rm=T) + 
                  sum(dnorm(x=log(snail.E+1), mean=log(1+sim.data[[i]]$RH[rows]/0.015), sd=sd.E, log=T), na.rm=T) + 
                  sum(dnorm(x=log(snail.W+1), mean=log(1+sim.data[[i]]$RP[rows]/4e-5), sd=sd.W, log=T), na.rm=T) + 
                  sum(log(sim.data[[i]]$SurvF[snail.death[snail]])) + #indexing to snail to get death date of focal individual
                  log(probs[i])
      }
    LL.ind_max = max(LL.ind)
    relative_probs = sum(exp(LL.ind - LL.ind_max))
    NLL = NLL - (LL.ind_max + log(relative_probs))
  }
  NLL
}


pars2 = pars[1:25]
pars2["sd.LM"] = 0.1
pars2["sd.L"] = pars["sd.LI1"]
pars2["sd.E"] = pars["sd.EI1"]
pars2["sd.W"] = pars["sd.W1"]
integrated.likelihood(pars2)

sd.LMs = 1:10/4
NLLs = numeric()
for(i in 1:length(sd.LMs)){
  pars2["sd.LM"] = sd.LMs[i]
  NLLs[i] = integrated.likelihood(pars2)
  print(i)
}
plot(sd.LMs, NLLs)

# 6 - data likelihood
full.likelihood<-function(x){
  
  
  
  # simulate data
  sim.data = make.states(x, in.F, in.S, dur.F, dur.S, Feed.F, w.t=7)
  
  # data likelihood
  e.c<-1
  
  ## observation model

  
  ## convert predictions into correct count units
  l.temp<-sim.data$L
  n.temp<-sim.data$RH/gammaH
  w.temp<-sim.data$RP/gammaP
  
  
  l2.temp<-sim.data$L2
  n2.temp<-sim.data$E2/gammaH
  w2.temp<-sim.data$W2/gammaP
  
  SF.temp<-sim.data$SurvF
  SS.temp<-sim.data$SurvS
  
  sd.LI1<-as.numeric(x["sd.LI1"])
  sd.LU1<-as.numeric(x["sd.LU1"])
  sd.EI1<-as.numeric(x["sd.EI1"])
  sd.EU1<-as.numeric(x["sd.EU1"])
  sd.W1<-as.numeric(x["sd.W1"])
  
  sd.LI2<-as.numeric(x["sd.LI2"])
  sd.LU2<-as.numeric(x["sd.LU2"])
  sd.EI2<-as.numeric(x["sd.EI2"])
  sd.EU2<-as.numeric(x["sd.EU2"])
  sd.W2<-as.numeric(x["sd.W2"])
  
  # Avoids simulations that fell short
  NObs = length(data$L)
  NObs2 = length(data$L2)
  if(length(n.temp) != NObs){print("Simulation too short"); return(-1e6)}
  if(length(n2.temp) != NObs2){print("Simulation too short"); return(-1e6)}    
  
  if(anyNA(l.temp)){print("NaNs in l.temp"); return(-1e6)}
  if(anyNA(l2.temp)){print("NaNs in l2.temp"); return(-1e6)}
  if(anyNA(SF.temp)){print("NaNs in SF.temp");return(-1e6)}
  if(anyNA(SS.temp)){print("NaNs in SS.temp"); return(-1e6)}
  
  sds = c(rep(sd.LI1, times = 2112), rep(sd.LU1, times=960))
  sd.Eggs = c(rep(sd.EI1, times = 2112), rep(sd.EU1, times=960))
  
  ## likelihoods from food gradient
  llik.L<- sum(dnorm(log(data$L), mean=log(l.temp), sd=sds, log=TRUE), na.rm=T)
  llik.Negg<- sum(dnorm(log(data$Negg+e.c), mean=log(n.temp+e.c), sd=sd.Eggs, log=TRUE), na.rm=T)
  llik.Nworms<- sum(dnorm(log(data$Nworms+e.c), mean=log(w.temp+e.c), sd=sd.W1, log=TRUE), na.rm=T)
  SF =SF.temp[which(data$Alive == 1)]
  llik.Survival <- sum(log(SF))
  
  sds2 = c(rep(sd.LI2, times = 380), rep(sd.LU2, times=380))
  sd.Eggs2 = c(rep(sd.EI2, times = 380), rep(sd.EU2, times=380))
  
  ## likelihoods from density gradient
  llik.L2<- sum(dnorm(log(data$L2), mean=log(l2.temp), sd=sds2, log=TRUE), na.rm=T)
  llik.Negg2<- sum(dnorm(log(data$Negg2+e.c), mean=log(n2.temp+e.c), sd=sd.Eggs2, log=TRUE), na.rm=T)
  llik.Nworms2<- sum(dnorm(log(data$Nworms2+e.c), mean=log(w2.temp+e.c), sd=sd.W2, log=TRUE), na.rm=T)
  
  
  SS =SS.temp[which(data$Alive2 == 1)]
  llik.Survival2 <- sum(log(SS)) 
  
  llik<-(llik.L + llik.Negg + llik.Nworms + llik.Survival + llik.L2 + llik.Negg2 + llik.Nworms2 + llik.Survival2)
  
  if(is.na(llik)|!is.finite(llik)){
    print("Infinite NLL")
    return(-1e6)}
  
  lprior = prior.likelihood(x)
  
  return( (llik + lprior)/1)
}




### Tuning ###

setwd("C:/RData")
samps = readRDS("FullStarve_shrink_adaptMCMC_original_DAM.Rda")
pars = samps$samples[which.max(samps$log.p),]
pars = as.numeric(pars[1:35])
names(pars) = c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2",
                "sd.LI1", "sd.LU1", "sd.EI1", "sd.EU1", "sd.W1", 
                "sd.LI2", "sd.LU2", "sd.EI2", "sd.EU2", "sd.W2")

full.likelihood(pars)
round(1 - rejectionRate(mcmc(samps$samples)), 2)

variances = samps$cov.jump
# 
# ### running the mcmc ###
# start.time = proc.time()
# #test = MCMC(full.likelihood, init=pars, scale=as.matrix(variances), adapt=50000, acc.rate = 0.3, n=250000)
# test = readRDS("FullStarve_shrink_adaptMCMC_original_DAM_run1.Rda")
# test = MCMC.add.samples(test, n.update=50000)
# 
# ### converting to coda
# testc = convert.to.coda(test)
# testc = cbind(testc, "lpost" = test$log.p)
# round(1 - rejectionRate(mcmc(testc)), 3)
# plot(1:length(testc[,"lpost"]), testc[,"lpost"])
# saveRDS(test, file="FullStarve_shrink_adaptMCMC_original_DAM_run1.Rda")
# 
# #test = MCMC.add.samples(test, n.update=10000)
