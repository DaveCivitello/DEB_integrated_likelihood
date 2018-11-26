#### General MCMC with BayesianTools package ####

library(Matrix)
library(deSolve)
library("mvtnorm")
library(LaplacesDemon)
library(coda)

rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

setwd("C:/RData")

# compile my model from C definition
dyn.unload("IndividualModel_Shrink_P_dilute.dll") # unload dll
system("R CMD SHLIB IndividualModel_Shrink_P_dilute.c")
dyn.load("IndividualModel_Shrink_P_dilute.dll") # Load dll


# # These results are really insensitive to futzing with initial energy density,
# # so I'm just going to set it to 'high'


# Set initial conditions for food gradient
setinits.Food<-function(F0 = 16.5, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

setinits.Starve<-function(F0 = 10.15, L0=0.7, e0=0.9, D0 = 0, RH0 = 0, P0 = 0, RP0 = 0, DAM0=0, HAZ0=0){
  inits<-c(F = F0, L = L0, e = e0, D = D0, RH = RH0, P = P0, RP = RP0, DAM=DAM0, HAZ=HAZ0)
  return(inits)
}

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


setup.mcmc.params<-function(hyper=NULL, w.p=NULL, samp.type=NULL,
                            prop.var=NULL, prop.mean=NULL, start=NULL,
                            joint=FALSE, n.joint=1, scale.var=NULL){
  
  
  ## first load the parameters for the priors
  if(is.null(hyper)) hyper<-make.hypers()
  
  ## Need to change these bits if parameters are only being proposed
  ## individually. Additionally, most of these can be changed by
  ## passing in a new vector with the appropriate name and format as
  ## an arguement to the function. However, most bits for any sets of
  ## parameters that are being proposed jointly must be changed
  ## manually below.
  
  if(is.null(w.p)) w.p <- c("k","EM", "Fh", "muD", "DR", 
                            "fe", "yRP", "ph", "yPE", "iPM", "eh", "mP",
                            "alpha", "yEF", "LM", "hb", "yVE", "yEF2",
                            "iM", "M", "kd", "z", "kk", "theta", "mR", "alpha.k",
                            "sd.L", "sd.E", "sd.W")
  if(is.null(samp.type)) samp.type<-list(iM="rw", k="rw", M="rw", EM="rw", Fh="rw", muD="rw", DR="rw", fe="rw",
                                         yRP="rw", ph="rw", yPE="rw", iPM = "rw", eh="rw", yEF="rw", LM="rw",
                                         kd="rw", z="rw", kk="rw", hb="rw", theta="rw", mR="rw", yVE="rw", yEF2="rw",
                                         mP="rw", alpha="rw", alpha.k="rw", sd.L="rw", sd.E="rw", sd.W="rw")
  if(is.null(prop.var)) prop.var<-list()                              
  if(is.null(prop.mean)) prop.mean<-list()
  if(is.null(start)) start<-c()
  
  ## The following lines must be changed if parameters are being
  ## proposed jointly
  
  if(length(n.joint)==1) n.joint<-seq(1, n.joint, by=1)
  
  # w1<-c()
  # cor1<-prev.matrix1
  # hyp1<-c(hyper$k, hyper$EM, hyper$Fh, hyper$muD, hyper$DR,
  #         hyper$fe, hyper$yRP, hyper$ph, hyper$yPE, hyper$iPM, hyper$eh, hyper$mP,
  #         hyper$alpha, hyper$yEF, hyper$LM,hyper$hb, hyper$yVE, hyper$yEF2)
  # m1<-c(1,1)
  # t1<-"rw"
  # 
  # w2<-c()
  # cor2<-prev.matrix2
  # hyp2<-c(hyper$iM, hyper$M)
  # m2<-c(1,1)
  # t2<-"rw"
  # 
  # w3<-c()
  # cor3<-prev.matrix3
  # hyp3<-c( hyper$kd, hyper$z, hyper$kk, hyper$theta, hyper$mR)
  # m3<-c(1,1)
  # t3<-"rw"
  # 
  # # 
  # # w2<-c()
  # # cor2<-prev.matrix2
  # # hyp2<-c(hyper$kk, hyper$kd, hyper$z, hyper$theta, hyper$hb)
  # # m2<-c(1,1)
  # # t2<-"rw"
  # # 
  # # 
  # ws<-list(w1, w2, w3)
  # ms<-list(m1, m2, m3)
  # cors<-list(cor1, cor2, cor3)
  # hyps<-list(hyp1, hyp2, hyp3)
  # types<-list(t1, t2, t3)
  ## The code below here builds the appropriate structure with all the
  ## info needed by the mcmc to propose samples, etc.
  
  all<-list()
  
  if(joint){
    n.joint<-n.joint[!is.na(n.joint)]
    j<-0
    for(i in n.joint){
      j<-j+1
      all[[j]] <- list(params = ws[[i]],
                       ##var = diag(x=c(0.00055, 0.0005), nrow=2),
                       mean   = ms[[i]],
                       var    = make.sigNN(var=prop.var[ws[[i]]],
                                           cor=cors[[i]],
                                           scale=scale.var),
                       hyp    = matrix(hyps[[i]], nrow=length(ws[[i]]),
                                       byrow=TRUE),
                       start  = start[ws[[i]]],
                       type   = types[[i]]
      )
    }
  }
  
  if(!joint) n.joint<-NULL
  if(!is.null(w.p)){
    n.params <- length(n.joint)+length(w.p)
    
    for(i in (length(n.joint)+1):n.params){
      all[[i]]<-list(params=w.p[i-length(n.joint)])
    }
    
    for(i in 1:n.params){
      pp<-all[[i]]$params
      if(length(pp)==1){
        if(is.element(pp, w.p)){
          all[[i]]$type<-samp.type[[pp]]
          all[[i]]$mean<-prop.mean[[pp]]
          all[[i]]$var<-prop.var[[pp]]
          all[[i]]$hyp<-hyper[[pp]]
          all[[i]]$start<-start[[pp]]
        }
        else{
          all[[i]]$var<-0
          ##all[[i]]$hyp<-0
          ##all[[i]]$
        }
      }
    }
    
  }
  return(all)
}

make.sigNN<-function(var=list(g=0.0009, eg=0.0009), cor, scale=NULL){
  
  if(!is.null(scale)) var.scale<-(2.38^2)/scale
  else var.scale<-1
  
  
  sd<-sqrt(as.numeric(unlist(var)))
  l<-length(var)
  
  sigma<-matrix(0, ncol=l, nrow=l)
  
  for(i in 1:l){
    for(j in 1:l){
      sigma[i,j]<-cor[i,j]*sd[i]*sd[j]*var.scale
    }
  }
  sigma<-as.matrix(nearPD(sigma)$mat)
  
  return(sigma)
}


#################### Model specification 2 ################
## This file holds many of the functions needed to simulate the DEB
## model forward in time. It needs the package deSolve to function.


# Discrete feeding events
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



solve.DEB.food<-function(params, inits, duration, feeding.events){
  feed.juv <- feeding.events
  feed.juv.U <- subset(feed.juv, var != "P")
  feed.juv.U[71, 2] = 28
  parms = as.numeric(params[1:24])
  params = c(iM=parms[1], k=parms[2], M=parms[3], EM=parms[4], Fh=parms[5], muD=parms[6],
             DR=parms[7], fe=parms[8], yRP=parms[9], kapc=parms[10], yPE=parms[11], iPM=parms[12],
             eh=parms[13], mP=parms[14], kapnew=parms[15], yEF=parms[16], LM=parms[17],kd=parms[18], 
             z=parms[19], kk=parms[20], hb=parms[21], theta=parms[22], mR=parms[23], yVE=parms[24], startage=28)
  
  Juv.R <- as.data.frame(lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                               params[1:25],  rtol=1e-6, atol=1e-6,   
                               events = list(data = feed.juv)))
  
  Juv.RU <- as.data.frame(lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                                params[1:25],  rtol=1e-6, atol=1e-6,   
                                events = list(data = feed.juv.U)))
  
  feed.juv[1:70,3] <- 11
  feed.juv.U[1:70,3] <- 11
  Juv.G <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.juv))
  Juv.GU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 5.5
  feed.juv.U[1:70,3] <- 5.5
  Juv.U <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.juv))  
  Juv.UU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 2.75
  feed.juv.U[1:70,3] <- 2.75
  Juv.P <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                 initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                 params[1:25],  rtol=1e-6, atol=1e-6,   
                 events = list(data = feed.juv))
  Juv.PU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 1.375
  feed.juv.U[1:70,3] <- 1.375
  Juv.BR <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv))
  Juv.BRU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                   initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                   params[1:25],  rtol=1e-6, atol=1e-6,   
                   events = list(data = feed.juv.U))
  
  feed.juv[1:70,3] <- 0.6875
  feed.juv.U[1:70,3] <- 0.6875
  Juv.BL <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                  initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"),maxsteps=5e5,
                  params[1:25],  rtol=1e-6, atol=1e-6,   
                  events = list(data = feed.juv))
  Juv.BLU <- lsoda(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
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
  
  out_1_0 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_1_0)))
  out_1_0U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_1_0.U)))
  
  # Simulate dynamics for 2-2
  food_2_2 = periodic.starvation.events(5.37, 2, 14, 18)
  food_2_2.U <- subset(food_2_2, var != "P")
  food_2_2[which(food_2_2$var=="HAZ"),2] = 35
  
  out_2_2 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_2_2)))
  out_2_2U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_2_2.U)))
  
  # Simulate dynamics for 3-3
  food_3_3 = periodic.starvation.events(8.06, 3, 14, 18)
  food_3_3.U <- subset(food_3_3, var != "P")
  food_3_3[which(food_3_3$var=="HAZ"),2] = 42
  
  out_3_3 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_3_3)))
  out_3_3U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                                initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                                params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                                events = list(data = food_3_3.U)))
  
  # Simulate dynamics for 4-4
  food_4_4 = periodic.starvation.events(10.74, 4, 14, 18)
  food_4_4.U <- subset(food_4_4, var != "P")
  food_4_4[which(food_4_4$var=="HAZ"),2] = 49
  
  out_4_4 <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
                               initfunc = "initmod",  nout=2, outnames=c("Survival", "LG"), method="lsoda",
                               params[1:25],  rtol=1e-6, atol=1e-6, maxsteps=5e5,
                               events = list(data = food_4_4)))
  out_4_4U <- as.data.frame(ode(inits, 0:duration, func = "derivs", dllname = "IndividualModel_Shrink_P_dilute", 
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

make.states.ILL<-function(params, inits.F, inits.S, duration.F, duration.S, feeding.events.F, w.t=7){
  n = 10
  alpha.k = as.numeric(params["alpha.k"])
  k.mean = as.numeric(params["k"])
  beta.k = alpha.k*(1 - k.mean)/k.mean
  qs = seq(from= 0.001, to = 0.999, length.out = n) #select quantiles
  ks = qbeta(qs, shape1=alpha.k, shape2=beta.k) # obtain random effect values
  ds = dbeta(ks, shape1=alpha.k, shape2=beta.k) # get the density of these values
  probs = ds/sum(ds)
  LLs = numeric()
  sim.data = list()
  for(i in 1:length(ks)){
    params2 = params
    params2["k"] = ks[i]
    sim.data[[i]] = make.states(params2, in.F, in.S, dur.F, dur.S, Feed.F, w.t=7)
  }
  sim.data
}


#################### Meat of the MCMC ##################################

## Bayesian inference for a deterministic DEB model (with models
## solved via an ode solver in the function specified with argument
## "sim") with an observation model, and with observation error
## variances specified in sds.

deb.mcmc<-function(N, data, all.params, inits.S, inits.F, samp.p,
                   cnt=10, burnin=0.1, duration.F, duration.S, feeding.events.F, 
                   plot=TRUE, w.t=7, test=TRUE, my.par=c(1,1))
{
  
  ## first we calculate a few lengths, etc, that we use for the for
  ## loops later. 
  
  l<-length(samp.p)
  np<-length(all.params)
  p<-NULL
  
  ltot<-0
  
  for(i in 1:l) ltot<-ltot+length(samp.p[[i]]$params)
  
  ## for testing, here is code that calculates (and prints out) the
  ## posterior prob of the "real" parameters, which can be passed in
  ## through all.params
  
  if(test){
    sim.old<-make.states.ILL(all.params, inits.F, inits.S, duration.F, duration.S, feeding.events.F, w.t)
    prob.old<-log.post.params(all.params, data, samp.p, sim.old)
    print(paste("posterior probability of the real parameters= ", prob.old, sep=""))
  }
  
  ## build a data frame to hold the posterior samples. These include
  ## "samples" even for the parameters that are being held fixed, to
  ## make the code more straightforward and generic, as well as being
  ## useful for debugging. The samps structure will also keep track of 
  ## the log posterior probability of that particular sample.
  
  samps<-data.frame(matrix(0, ncol=np+1, nrow=N))
  names(samps)<-c(names(all.params), "lpost")
  
  ## initialize the samps structure and p
  samps[1,1:np]<-p<-all.params  
  
  ## initialize the data frame with the starting values passed into
  ## the function in the structure samp.p
  for(i in 1:l){
    samps[1, samp.p[[i]]$params]<-samp.p[[i]]$start
    p[samp.p[[i]]$params]<-samp.p[[i]]$start
  }
  
  
  ## run the data simulation to make the underlying states (for
  ## determining the likelihoods) using the parameters stored in p.
  
  sim.old<-make.states.ILL(p, inits.F, inits.S, duration.F, duration.S, feeding.events.F, w.t)
  
  ## check the posterior probability to make sure you have reasonable
  ## starting values, and to initialize prob.old
  samps$lpost[1]<-log.post.params(samps[1,], data, samp.p, sim.old)
  
  print(paste("initial posterior probability = ", samps$lpost[1], sep=""))
  
  if(!is.finite(samps$lpost[1])) stop("bad starting values")
  
  ## now we begin the MCMC
  
  out<-list(s=samps[1,], p=p, sim.old=sim.old)
  
  for(i in 2:N){
    ## printing and plotting output so we can watch the progress
    if(i%%cnt == 0){
      print(paste("sample number", i, sep=" "))
      if(plot) plot.output(i, samps, samp.p, l, ltot, my.par)
    }
    
    ## the meat of the MCMC is found in the function update.samps (see below)
    
    out<-update.sample(samps[i-1,], samp.p, data, inits.F, inits.S, out,
                       duration.F, duration.S, w.t, l, i, cnt, feeding.events.F)
    samps[i,]<-out$s
    if(test){
      if(-samps$lpost[i-1]+samps$lpost[i-1]<=-10){
        stop("we've had a really large swing in the posterior prob")
      }
    }
    
  }
  
  ##plot(b[100:N],type="l")
  lim<-min(1, round(burnin*N))
  samps <- samps[lim:N,]
  
  return(list(samps=samps))
  
}



update.sample<-function(samps, samp.p, data, inits.F, inits.S, out, duration.F, duration.S, w.t = 7,
                        l, i, cnt, feeding.events.F, test=TRUE)
{
  ## read in some bits
  s<-samps
  sim.old<-out$sim.old
  p<-out$p
  
  x<-1:l
  s.x<-sample(x)
  
  for(k in s.x){   
    
    s.new<-s
    p.new<-p
    
    q<-propose.params(s, samp.p[[k]])
    
    ## automatically reject if the proposed parameter value is
    ## outside of the reasonable limits (i.e. < 0 )
    zeros<-0
    zeros<-check.zeros(samp.p[[k]], q$b)
    if(zeros){
      problem = samp.p[[k]]$params[which(q$b < 0 | q$b > 1)]
      print(paste("proposed", toString(problem), "out of range. moving on"), sep="")
      next
    }
    ## write the proposed params into the p.new and s.new.
    
    for(j in 1:length(samp.p[[k]]$params)){
      ww<-samp.p[[k]]$params[j]
      p.new[ww]<-s.new[ww]<-q$b[j]
    }
    
    ## simulate the dynamics forward with the new parameters
    sds = c("sd.L", "sd.E", "sd.W") # Allows us not to simulate dynamics for observation parameters
    #sim.new<-make.states(sim, p.new, inits, duration, feeding.events, w.t)
    ifelse(ww %in% sds, sim.new<-sim.old, sim.new<-make.states.ILL(p.new, inits.F, inits.S, duration.F, duration.S, feeding.events.F, w.t))
    ## The posteriorprob of the previous sample is saved as
    ## s$lpost. If we accept a draw, we will set s$lpost<-s.new$lpost
    
    s.new$lpost<-log.post.params(s.new, data, samp.p, sim.new)
    
    if(is.finite(s.new$lpost) && is.finite(s$lpost)){
      A<-exp( s.new$lpost + q$lbak - s$lpost - q$lfwd )
    }
    else{
      A<-0
      print("whoops! must have proposed outside the correct region")
    }
    
    ## print some output so we can follow the progress
    if(i%%cnt==0){
      print(paste("proposing " , samp.p[[k]]$params, ": prob.old = ",
                  signif(s$lpost, digits=5),
                  "; prob.new = ", signif(s.new$lpost, digits=5), "; A = ",
                  signif(A, digits=5),
                  sep=""))
    }
    
    ## take a draw from a unif distrib, and use it to accept/reject
    u<-runif(1)    
    if( u < A ){ ## accept
      sim.old<-sim.new
      p<-p.new
      s<-s.new
    }
  }
  
  return(list(s=s, p=p, sim.old=sim.old))
  
}


check.zeros<-function(s.p, q.b){
  z<-0
  for(j in 1:length(s.p$params)){
    if(s.p$params[j]=="l.M.HP") next
    if(q.b[j]<0)z<-1
    if(s.p$params[j] %in% c("yVE", "yEF", "fe", "k", "yPE", "yRP", "mP")){
      if(q.b[j]>1) z<-1
    }
  }
  return(z)
}




propose.params<-function(samps, s.p)
{
  if(length(s.p$params)==1){
    ##print(paste(s.p$params, " proposing single ", sep=" "))
    q<-propose.single(samps, s.p)
  }
  else{
    ##print(paste(s.p$params, " proposing jointly ", sep=" "))
    q<-propose.joint(samps, s.p)
  }         
  return(q)
  
}

## I'm feeding in the variance, so I need to take the square root....
propose.single<-function(samps, s.p)##, i, freq=50, size=50 )##l=5, h=6)
{
  
  b<-as.numeric(samps[s.p$params])
  var<-s.p$var
  type<-s.p$type
  hyps<-s.p$hyp
  
  
  if(type=="rw"){
    if(length(var)>1){
      l<-var[1]
      h<-var[2]
      b.new<-runif(1, l/h*b, h/l*b)
      lfwd<-dunif(b.new, l/h*b, h/l*b, log=TRUE)
      lbak<-dunif(b, l/h*b.new, h/l*b.new, log=TRUE)
    }
    else{
      sd<-sqrt(var)
      b.new<-rnorm(1, b, sd=sd)
      lfwd<-dnorm(b.new, b, sd=sd, log=TRUE)
      lbak<-dnorm(b, b.new, sd=sd, log=TRUE)
      #if(s.p$params == "M"){print(b-b.new)}
    }
    return(list(b=b.new, lfwd=lfwd, lbak=lbak))
  }
  else if(type=="ind"){
    out<-prior.draw(b, hyps, s.p$params)
    return(out)
  }
  
}  

# UP = function(proportion){ # Easy uniform proposal function
#   c(1 - proportion/2, 1 + proportion/2)
# }

propose.joint<-function(samp, samp.p){
  
  
  b<-NULL
  if(samp.p$type=="rw"){
    b<-as.numeric(samp[samp.p$params])
    sigma<-samp.p$var
    
    b.new<-rmvnorm(1, mean=b, sigma=sigma)
    lfwd<-dmvnorm(b.new, b, sigma, log=TRUE)
    lbak<-dmvnorm(b, b.new, sigma, log=TRUE)
  }
  else if(samp.p$type=="ind"){
    if(is.null(samp.p$mean)) stop("not enough info for the independence sampler")
    mean<-as.numeric(samp.p$mean)
    b<-as.numeric(samp[samp.p$params])
    sigma<-samp.p$var
    
    b.new<-rmvnorm(1, mean=mean, sigma=sigma)
    lfwd<-dmvnorm(b.new, mean, sigma, log=TRUE)
    lbak<-dmvnorm(b, mean, sigma, log=TRUE)
  }
  
  ##print(c(b, b.new, lfwd, lbak))
  
  
  ##samp[s]<-b.new
  
  ##stop()
  return(list(b=b.new, lfwd=lfwd, lbak=lbak))
  
}



plot.output<-function(i, samps, samp.p, l, ltot, my.par=c(2,4), plot.post=TRUE){
  
  if( ltot > 1 ) par(mfrow=my.par, bty="n")
  for( j in 1:29 ){ # Plots only the 25 process parameters
    ww<-samp.p[[j]]$params
    for(k in 1:length(ww)){
      plot(samps[1:(i-1),ww[k]], type="l", xlab="sample", ylab=ww[k])
    }
  }
  if(plot.post){
    my.min<-100
    if(i>my.min){
      x<-seq(my.min, (i-1), by=1)
      plot(x, samps$lpost[x], type="l", xlab="sample", ylab="log posterior prob")
    }
  }
}


prior.draw<-function(b, hyp, p){
  
  param1<-hyp[1]
  param2<-hyp[2]
  
  gams<-c()
  norms<-c("LM", "iM", "M")
  betas<-c("yVE", "yEF", "yEF2", "yRP", "fe", "k", "yPE", "mP", "eh")
  unifs<-c("EM", "DR", "Fh", "iPM", "alpha", "muD", "kd", "z", "hb", "kk", "theta", "mR", "ph",
           "alpha.k", "sd.L", "sd.E", "sd.W")
  
  if( p %in% gams ){
    b.new<-rgamma(1, shape=param1, rate=param2)
    lfwd<-dgamma(b.new, shape=param1, rate=param2, log=TRUE)
    lbak<-dgamma(b, shape=param1, rate=param2, log=TRUE)
  }
  else if( p %in% norms ){
    b.new<-rnorm(1, mean=param1, sd=param2)
    lfwd<-dnorm(b.new, mean=param1, sd=param2, log=TRUE)
    lbak<-dnorm(b, mean=param1, sd=param2, log=TRUE)
  }
  else if( p %in% betas ){
    b.new<-rbeta(1,  shape1=param1, shape2=param2)
    lfwd<-dbeta(b.new,  shape1=param1, shape2=param2, log=TRUE)
    lbak<-dbeta(b, shape1=param1, shape2=param2, log=TRUE)
  }
  else if( p %in% unifs ){
    b.new<-runif(1,  param1, param2)
    lfwd<-dunif(b.new,  param1, param2, log=TRUE)
    lbak<-dunif(b, param1, param2, log=TRUE)
  }
  
  return(list(b=b.new, lfwd=lfwd, lbak=lbak)) 
}

## In this file we put the functions that calculate prior and
## posterior probabilities of the parameters. The function that
## returns the posterior probabilities must be called
## "log.post.params", and the prior functions should in
## "log.prior.params" (which will only be called from
## log.post.params). The user can update these functions so they are
## appropriate for their data and parameters. Additionally the
## function "make.hypers" makes an appropriately structured set of
## hyper parameters (which determine the prior) for the deb.mcmc
## function.


log.post.params<-function(samps, data, samp.p, sim.data){
  
  e.c<-1
  
  ## observation model
  gammaH<-0.015
  gammaP<-4e-5
  sd.L <- as.numeric(samps["sd.L"])
  sd.E <- as.numeric(samps["sd.E"])
  sd.W <- as.numeric(samps["sd.W"])
  
  # probabilities for integrating over random effect in parameter
  n = 10
  alpha.k = as.numeric(samps["alpha.k"])
  k.mean = as.numeric(samps["k"])
  beta.k = alpha.k*(1 - k.mean)/k.mean
  qs = seq(from= 0.001, to = 0.999, length.out = n) #select quantiles
  ks = qbeta(qs, shape1=alpha.k, shape2=beta.k) # obtain random effect values
  ds = dbeta(ks, shape1=alpha.k, shape2=beta.k) # get the density of these values
  probs = ds/sum(ds)

  LL = 0
  # Scrolling over snails in first experiment
  for(snail in 1:96){ # subsetting by snails in supply gradient experiment
    LL.ind = numeric()
    rows = 1:32 + (snail-1)*32
    snail.L = data$L[rows]
    snail.E = data$Negg[rows]
    snail.W = data$Nworms[rows]
    snail.death = which(data$Alive == 1)
    for(i in 1:n){ # scrolling over different values of random effect parameter
      if(length(sim.data[[i]]$L[rows]) != length(snail.L) | ks[i] <= 0 | anyNA(sim.data[[i]]$RH[rows])){
        LL.ind[i] = -1e6; #print("failed in 1st exp")
        next
      }
      LL.ind[i] = (sum(dnorm(x=log(snail.L), mean=log(sim.data[[i]]$L[rows]), sd=sd.L, log=T), na.rm=T) + 
                     sum(dnorm(x=log(snail.E+1), mean=log(1+sim.data[[i]]$RH[rows]/0.015), sd=sd.E, log=T), na.rm=T) + 
                     sum(dnorm(x=log(snail.W+1), mean=log(1+sim.data[[i]]$RP[rows]/4e-5), sd=sd.W, log=T), na.rm=T) + 
                     #indexing to snail to get death date of focal individual
                     sum(log(sim.data[[i]]$SurvF[snail.death[snail]])) + log(probs[i]))
    }
    LL.ind_max = max(LL.ind)
    relative_probs = sum(exp(LL.ind - LL.ind_max))
    LL = LL + (LL.ind_max + log(relative_probs))
  }
  # Scrolling over snails in second experiment
  for(snail in 1:40){ # subsetting by snails in supply gradient experiment
    LL.ind = numeric()
    rows = 1:19 + (snail-1)*19
    snail.L = data$L2[rows]
    snail.E = data$Negg2[rows]
    snail.W = data$Nworms2[rows]
    snail.death = which(data$Alive2 == 1)
    for(i in 1:n){ # scrolling over different values of random effect parameter
      if(length(sim.data[[i]]$L2[rows]) != length(snail.L) | ks[i] <= 0  | anyNA(sim.data[[i]]$E2[rows])){
        LL.ind[i] = -1e6; #print("failed in 2nd exp")
        next
      }
      LL.ind[i] = (sum(dnorm(x=log(snail.L), mean=log(sim.data[[i]]$L2[rows]), sd=sd.L, log=T), na.rm=T) +
                     sum(dnorm(x=log(snail.E+1), mean=log(1+sim.data[[i]]$E2[rows]/0.015), sd=sd.E, log=T), na.rm=T) +
                     sum(dnorm(x=log(snail.W+1), mean=log(1+sim.data[[i]]$W2[rows]/4e-5), sd=sd.W, log=T), na.rm=T) +
                     #indexing to snail to get death date of focal individual
                     sum(log(sim.data[[i]]$SurvS[snail.death[snail]])) + log(probs[i]))
      
    }
    LL.ind_max = max(LL.ind)
    relative_probs = sum(exp(LL.ind - LL.ind_max))
    LL = LL + (LL.ind_max + log(relative_probs))
  }

  lprior<- sum(log.prior.params(samps, samp.p))
  
  if(is.na(LL)|!is.finite(LL))return(-1e6)#stop("something is wrong in the likelihood bit")
  if(is.na(lprior)|!is.finite(lprior)) return(-1e6)#stop("something is wrong in the prior bit")
  
  return(LL + lprior) # trying this for a second
}


## eventually change this to make it more efficient
log.prior.params<-function(samps, samp.p){
  
  lp<-NULL
  cnt<-0
  len<-length(samp.p)
  
  ##print(c(as.numeric(samp), w.p, hyper[1]))
  ## this won't work as it is for joint proposals
  for(i in 1:len){
    l<-length(samp.p[[i]]$params)
    for(j in 1:l){
      cnt<-cnt+1
      p<-samp.p[[i]]$params[j]
      s<-as.numeric(samps[p])
      
      if(l==1){
        param1<-samp.p[[i]]$hyp[1]
        param2<-samp.p[[i]]$hyp[2]
      }
      else{
        param1<-samp.p[[i]]$hyp[j,1]
        param2<-samp.p[[i]]$hyp[j,2]
      }
      
      if( p == "yPE" ) lp<-c(lp, yPE=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "yEF" ) lp<-c(lp, yEF=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "yEF2" ) lp<-c(lp, yEF2=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "yRP" ) lp<-c(lp, yRP=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "yVE" ) lp<-c(lp, yEF=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      
      if( p == "sd.L" )  lp<-c(lp, sd.LI1=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "sd.E" )  lp<-c(lp, sd.EU1=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "sd.W" )  lp<-c(lp, sd.W1=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "alpha.k" )  lp<-c(lp, sd.W1=dunif(s, min=param1, max=param2, log=TRUE))

      if( p == "ph" ) lp<-c(lp, ph=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "alpha" ) lp<-c(lp, alpha=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "mP" ) lp<-c(lp, mP=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "eh" ) lp<-c(lp, eh=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "iPM" )  lp<-c(lp, iPM=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "iM" ) lp<-c(lp, iM=dnorm(s, mean=param1, sd=param2,  log=TRUE))
      if( p == "k" ) lp<-c(lp, k=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))      
      if( p == "EM" ) lp<-c(lp, EM=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "M" ) lp<-c(lp, M=dnorm(s, mean=param1, sd=param2,  log=TRUE))
      if( p == "LM" ) lp<-c(lp, LM=dnorm(s, mean=param1, sd=param2,  log=TRUE))
      if( p == "DR" ) lp<-c(lp, DR=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "Fh" ) lp<-c(lp, Fh=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "fe" ) lp<-c(lp, fe=dbeta(s, shape1=param1, shape2=param2,  log=TRUE)) 
      if( p == "muD" ) lp<-c(lp, muD=dunif(s, min=param1, max=param2, log=TRUE)) 
      if( p == "kd" ) lp<-c(lp, kd=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "z" ) lp<-c(lp, z=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "kk" ) lp<-c(lp, kk=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "hb" ) lp<-c(lp, hb=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "theta" ) lp<-c(lp, theta=dunif(s, min=param1, max=param2, log=TRUE))
      if( p == "mR" ) lp<-c(lp, mR=dunif(s, min=param1, max=param2, log=TRUE))
      
      if(!is.finite(lp[cnt])){lp[cnt]<- -1e6} #stop("one of the prior probs isn't finite") 
    }
  }
  
  return(lp)
  
}


## here's a function to make the hyper params in the correct format
make.hypers<-function(yPE=c(1,1), yEF=c(1,1), yEF2=c(1,1), LM=c(35,2), yRP=c(1,1), 
                      kap.b=c(0,10), ph=c(0, 100), eh=c(1,1), mP=c(1,1), iPM=c(0, 10), alpha=c(0,10), yVE=c(1,1),
                      iM=c(0.0183, 0.0016), k=c(1,1),  EM=c(0,1000), M=c(0.004, 0.00047), DR=c(0,10), kd=c(0,10), z=c(0,10), kk=c(0,10), hb=c(0,10), theta=c(0,10000),
                      Fh=c(0, 10), fe=c(1,1), muD=c(0,10), mR=c(0,10), alpha.k=c(0,10000),
                      sd.L=c(0,10), sd.E=c(0,10), sd.W=c(0,10)){
  
  hyper<-list(yPE=yPE, yEF=yEF, yEF2=yEF2, LM=LM, yRP=yRP, kap.b=kap.b, ph=ph, 
              eh=eh, mP=mP, iPM=iPM, alpha=alpha, iM=iM, k=k, EM=EM, M=M, DR=DR, kd=kd, z=z, kk=kk, hb=hb, theta=theta, mR=mR, yVE=yVE,
              Fh=Fh, fe=fe, muD=muD, alpha.k=alpha.k,
              sd.L=sd.L, sd.E=sd.E, sd.W=sd.W)
  
  return(hyper)
  
}


######## Tuning #########
setwd("C:/RData")
samps = readRDS("ILL_shrink_damage.Rda")

samps <- as.mcmc(samps[,c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                          "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                          "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2", "alpha.k")])

d = as.dist(1 - abs(cor(samps)))
par(mfrow=c(1,1))
plot(hclust(d))

samps = readRDS("ILL_shrink_damageA.Rda")
samps <- as.mcmc(samps[,c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                          "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                          "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2", "alpha.k", "sd.L", "sd.E", "sd.W",
                          "lpost")])

round(1 - rejectionRate(mcmc(samps)), 2)

variances = apply(samps, 2, var)

de_vars = c(iM=1e-17, k=1e-14, M=4e-19, EM=1e-15, Fh=2e-20, muD=3e-13,
            DR=3e-12, fe=2e-14, yRP=5e-13, ph=4e-16, yPE=3e-15, iPM=1e-15,
            eh=4e-16, mP=2e-16, alpha=3e-13, yEF=3e-17, LM=6e-14, kd=1e-13, 
            z=1e-14, kk=5e-11, hb=2e-12, theta=5e-9, mR=1e-16, yVE=4e-15, 
            yEF2=3e-9)

obs_vars = c(alpha.k=4e-5, sd.L=0.0001, sd.E=0.001, sd.W=0.01)
propvar = c(de_vars, obs_vars)




pars = as.vector(data.frame(samps)[max(which(data.frame(samps)$lpost >= max(data.frame(samps)$lpost) -0.001)),])
names(pars) = c("iM", "k", "M", "EM", "Fh", "muD", "DR", "fe", "yRP",
                "ph", "yPE", "iPM", "eh", "mP", "alpha", "yEF", "LM",
                "kd", "z", "kk", "hb", "theta", "mR", "yVE", "yEF2", "alpha.k", "sd.L", "sd.E", "sd.W")

start.time = proc.time()
## here's code to run a sample of our nifty mcmc
## run with: nohup R CMD batch run_mcmc.R &

#set.seed(15)
plt=T
my.par<-c(5,6)

N<-100
cnt<-10


# To avoid having to press RETURN
par(ask=F)
par(mar = c(2, 4, 2, 2))



setwd("C:/RData")
data = read.csv("AllSnailsLTSurv.csv")
data2 = read.csv("schisto_periodic_starve_LT2.csv")
data = list(t = data$Date, L = data$Length, Negg = data$C_Eggs, Nworms = data$C_Worms, Alive=data$Last_Alive,
            L2 = data2$Length, Negg2 = data2$C_Eggs, Nworms2 = data2$C_Worms, Alive2=data2$Last_Alive)


dur.F = 245
dur.S = 140



params = c(pars, gammaH=0.015, gammaP=4e-5)
hyps = make.hypers(yPE=c(1,1), yEF=c(1,1), yEF2=c(1,1), LM=c(35,2), yRP=c(1,1), sd.L=c(0,10), sd.E=c(0,10), sd.W=c(0,10),
                   kap.b=c(0,10), ph=c(0,100), eh=c(1,1), mP=c(1,1), iPM=c(0, 10), alpha=c(0,10), mR=c(0,10), yVE=c(1,1), alpha.k=c(0,10000),
                   iM=c(0.0183, 0.0016), k=c(1,1),  EM=c(0,1000), M=c(0.004, 0.00047), DR=c(0,10), kd=c(0,100), z=c(0,10), kk=c(0,10), hb=c(0,10), theta=c(0,100000),
                   Fh=c(0, 10), fe=c(1,1), muD=c(0,10))

mcmc.p<-setup.mcmc.params(start=pars, prop.var = propvar, hyper=hyps, joint=F)


in.F = setinits.Food()
in.S = setinits.Starve()
Feed.F = feeding.events(dates = sort(c((1:35)*7,((1:35)*7 - 3))), var="F", in.F[1], method="replace", Infected=28)

samps.2<-deb.mcmc(N, data, params, inits.F=in.F, inits.S=in.S, mcmc.p, cnt, 
                  burnin=0.1, duration.F=dur.F, duration.S=dur.S, feeding.events.F=Feed.F, 
                  plot=T, w.t=7, test=F, my.par=c(5,6))

samps<-samps.2$samps
saveRDS(samps, file="ILL_shrink_damageA.Rda")

print(proc.time() - start.time)