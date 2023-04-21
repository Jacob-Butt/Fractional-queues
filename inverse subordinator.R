######################################################
# This code is used to generate plots of our functional limits for the GI/Gi/1 queue model, using the formulation in terms
# of inverse subordinators. We assume that alpha_1 = alpha_2 = alpha for our two fractional Poisson Processes in this case, with arrival "rates" 
# lambda and departures "rates" mu

tstep = 0.01 #time partition size
xstep = 0.01 # space partition size
xlim = 10^3 # spatial limit
nsteps = 10^5 #number of steps in the algorithm

######################################
# Initialisation of vectors for our subordinators and inverse subordinators. XA and XD represent vectors for the subordinators of the arrival and departure processes,
# while YA and YD are vectors of their inverses (that is, YA will be the inverse of XA and so on).
# which will be used to gen
XA <- c(0, 0) 
XD <- c(0,0)
YA <- c(0,0)
YD <- c(0,0)

Q <- c(0,0) # vector initialisation for the functional limit given in the final line of (3.20)

alpha = 0.9 #shape parameter
lambda = 1
mu = 1
t = 0 # time


k = 0 # dummy variables keeping track of number of steps in the algorithm
xa = 0 # initial values for the subordinators at time 0
xd = 0

#######################################
# algorithm to generate subordinators XA and XD. This is done by using the CMS formula for Levy subordinators.
########################################

while (k < nsteps){
  k = k + 1
  uA <- runif(1)
  vA <- runif(1)
  uD <- runif(1)
  vD <- runif(1)
  phiA = pi * (vA - 1/2)
  phiD = pi * (vD - 1/2)
  L =  tstep^(1/alpha) * (1 + (tan((pi * alpha)/2))^2)^(1/(2 * alpha))
  xStepA = L * (sin(alpha * pi * vA)/(cos(phiA))^(1/alpha)) * (cos(phiA - alpha * pi * vA)/(-log(uA)))^((1-alpha)/alpha)
  xa = xa + xStepA
  xStepD = L * (sin(alpha * pi * vD)/(cos(phiD))^(1/alpha)) * (cos(phiD - alpha * pi * vD)/(-log(uD)))^((1-alpha)/alpha)
  xd = xd + xStepD
  t = t + tstep
  XA = rbind(XA, c(t,xa))
  XD = rbind(XD, c(t,xd))
}

j = 1 # dummy variables to keep track of algorithm steps
k = 1

fA = 0 # These are used to denote the total number of steps taken by the two algorithms below
fD = 0


################################################################
# two algorithms to generate the inverse subordinators YA and YD from XA, XD respectively. These are generated directly using the definition of inverse functions
#############################################################

while (j < nsteps){
  kStep = k * xstep
  if (kStep > xlim){
    break
  }
  z = XA[j,2]
  if (z > kStep){
    i = XA[j,1]
    YA = rbind(YA, c(kStep, i2))
    fA = k
    k = k + 1
  }
  else{
    j = j + 1
  }
}

k = 1
j = 1


while (j < nsteps){
  kStep = k * xstep
  if (kStep > xlim){
    break
  }
  z = XD[j,2]
  if (z > kStep){
    i = XD[j,1]
    YD = rbind(YD, c(kStep, i))
    fD = k
    k = k + 1
  }
  else{
    j = j + 1
  }
}

#######################################################################
# Final section that fills in vector Q using the values of lambda, mu and the inverse suboridnators YA, YD
#######################################################################

k = 1
f = min(fA, fD)
R = 0
Rdiff = 0
infR = 0

while (k < f){
  kStep = k * xstep
  R = (lambda^(alpha)) * YA[k, 2] - (mu^(alpha)) * YD[k, 2]
  if(R < infR){
    infR = R
  }
  Rdiff = R - infR
  Q = rbind(Q, c(kStep, Rdiff))
  k = k + 1
}

plot(Q, type="l", xlab = "time", ylab = "X") # plots Q as a function of time