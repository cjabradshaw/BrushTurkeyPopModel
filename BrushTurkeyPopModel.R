##########################################################################################################################################
## brush turkey (Alectura lathami) demographic model
## 
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, September 2021
##########################################################################################################################################

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


##############################
## ALECTURA (lathami) (AL)

# mass
AL.mass <- 2.2 # Alectura lathami (brush turkey) adult female (Jones, R.W.R.J. Dekker, C.S. Roselaar. The Megapodes: Megapodiidae, Oxford University Press, Oxford, New York, Tokyo (1995))
## for comparison, female malleefowl range from 1.68-2.235 kg (Leipoa ocellata) Priddel & Wheeler 2003 Wildl Res 30:451-464

## Omnivorous birds (Juanes 1986 Am Nat 128: 921-929)
AL.D.pred <- (10^(1.63 - 0.23*(log10(AL.mass*1000))))/2

# non-passerines (tmax=7.02*M^0.174) (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
AL.age.max <- round(7.02*(AL.mass*1000)^0.174,0) # reproductive senescence at 25 years for malleefowl (Bode & Brennan 2011-Oryx 45:513-520) 

## age vector
AL.age.vec <- 0:AL.age.max

## fertility
AL.eggs <- 16.6 / 2 # (for females only)
AL.hatch <- 0.866
AL.F.pred <- AL.eggs*AL.hatch

## age at primiparity
## alpha (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
AL.alpha <- ceiling(0.214*(AL.mass*1000)^0.303) # for birds
AL.alpha <- 2 # Bode & Brennan 2011 (Frith 1959)

## define m function with age
AL.m.vec <- c(rep(0, AL.alpha-1), rep(0.5*AL.F.pred, round(AL.alpha/2,0)), 0.75*AL.F.pred, rep(AL.F.pred, (AL.age.max+1-((AL.alpha+round(AL.alpha/2,0))))))
AL.m.sd.vec <- 0.05*AL.m.vec
plot(AL.age.vec, AL.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# Bleasdale yield density function y = x(a + bx^θ)^(-1/θ)
AL.m.dat <- data.frame(AL.age.vec, AL.m.vec)
param.init <- c(8.8077e-03, 4.571e-04, 3.89)
AL.fit.byd <- nls(AL.m.vec ~ AL.age.vec * ((a + b*(AL.age.vec)^theta)^(-1/theta)), 
                   data = AL.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], theta = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
plot(AL.age.vec, AL.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
AL.age.vec.cont <- seq(0,max(AL.age.vec),1)
AL.pred.p.mm <- AL.age.vec.cont * ((coef(AL.fit.byd)[1] + coef(AL.fit.byd)[2]*(AL.age.vec.cont)^coef(AL.fit.byd)[3])^(-1/coef(AL.fit.byd)[3]))
lines(AL.age.vec.cont, AL.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -1.78; b.s <- -0.21 # for birds
AL.s.tran <- ln.a.s + b.s*log(AL.mass*1000) + log(1)
AL.s.ad.yr <- exp(-exp(AL.s.tran))

## calculate lmax from Dillingham et al. 2016 Ecol Appl 26:322-333
# lmax_DIM occurs when lmax.fun returns 0
lmax.fun <- function(lmax, alpha, s, ar.aT) {
  abs(lmax - exp(ar.aT*(alpha + s/(lmax-s) )^(-1) ))
}
# get.lmax calculates lmax_DIM for scalar alpha and scalar/vector s
get.lmax <- function(alpha, s, ar.aT) {
  # Fixed alpha, allows a single value or vector for s
  lmax <- rep(NA, length(s))
  for (i in 1:length(s)) {
    lmax[i] <- optimize(lmax.fun, c(1, 5), tol = 0.0001, alpha = alpha, s=s[i], ar.aT=ar.aT)$minimum
  }
  return(list(lmax=lmax, alpha=alpha, s=s, ar.aT=ar.aT))
}

AL.lm.pred <- (get.lmax(alpha=AL.alpha, s=AL.s.ad.yr, ar.aT=1.107))$lmax # for birds
AL.rm.pred <- log(AL.lm.pred)

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.56*AL.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 0.3 # rate of mortality decline (also known as bt)
a2 <- 1 - AL.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 2.5e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.31 # rate of mortality increase
longev <- AL.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
AL.Sx <- c(0.1*AL.s.ad.yr, 1 - qx)
plot(x, AL.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
AL.s.sd.vec <- 0.05*AL.Sx

## create matrix
AL.popmat <- matrix(data = 0, nrow=AL.age.max+1, ncol=AL.age.max+1)
diag(AL.popmat[2:(AL.age.max+1),]) <- AL.Sx[-(AL.age.max+1)]
AL.popmat[AL.age.max+1,AL.age.max+1] <- AL.Sx[AL.age.max+1]
AL.popmat[1,] <- AL.pred.p.mm
colnames(AL.popmat) <- c(0:AL.age.max)
rownames(AL.popmat) <- c(0:AL.age.max)
AL.popmat.orig <- AL.popmat ## save original matrix

## matrix properties
max.lambda(AL.popmat.orig) ## 1-yr lambda
AL.lm.pred
max.r(AL.popmat.orig) # rate of population change, 1-yr
AL.ssd <- stable.stage.dist(AL.popmat.orig) ## stable stage distribution
plot(AL.age.vec, AL.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(AL.popmat.orig, AL.age.max) # reproductive value
AL.gen.l <- G.val(AL.popmat.orig, AL.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
AL.pop.found <- round(area*AL.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
AL.init.vec <- AL.ssd * AL.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*AL.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

AL.tot.F <- sum(AL.popmat.orig[1,])
AL.popmat <- AL.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
AL.n.mat <- matrix(0, nrow=AL.age.max+1,ncol=(t+1))
AL.n.mat[,1] <- AL.init.vec

## set up projection loop
for (i in 1:t) {
  AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
}

AL.n.pred <- colSums(AL.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(AL.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
AL.K.max <- 1*AL.pop.found
AL.K.vec <- c(1, AL.K.max/2, 0.75*AL.K.max, AL.K.max) 
AL.red.vec <- c(1,0.96,0.88,0.805)
plot(AL.K.vec, AL.red.vec,pch=19,type="b")
AL.Kred.dat <- data.frame(AL.K.vec, AL.red.vec)

# logistic power function a/(1+(x/b)^c)
AL.param.init <- c(1, 2*AL.K.max, 2)
AL.fit.lp <- nls(AL.red.vec ~ a/(1+(AL.K.vec/b)^c), 
                 data = AL.Kred.dat,
                 algorithm = "port",
                 start = c(a = AL.param.init[1], b = AL.param.init[2], c = AL.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
AL.fit.lp.summ <- summary(AL.fit.lp)
plot(AL.K.vec, AL.red.vec, pch=19,xlab="N",ylab="reduction factor")
AL.K.vec.cont <- seq(1,2*AL.pop.found,1)
AL.pred.lp.fx <- coef(AL.fit.lp)[1]/(1+(AL.K.vec.cont/coef(AL.fit.lp)[2])^coef(AL.fit.lp)[3])
lines(AL.K.vec.cont, AL.pred.lp.fx, lty=3,lwd=3,col="red")

AL.a.lp <- coef(AL.fit.lp)[1]
AL.b.lp <- coef(AL.fit.lp)[2]
AL.c.lp <- coef(AL.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
AL.n.mat <- matrix(0, nrow=AL.age.max+1, ncol=(t+1))
AL.n.mat[,1] <- AL.init.vec
AL.popmat <- AL.popmat.orig

## set up projection loop
for (i in 1:t) {
  AL.totN.i <- sum(AL.n.mat[,i])
  AL.pred.red <- as.numeric(AL.a.lp/(1+(AL.totN.i/AL.b.lp)^AL.c.lp))
  diag(AL.popmat[2:(AL.age.max+1),]) <- (AL.Sx[-(AL.age.max+1)])*AL.pred.red
  AL.popmat[AL.age.max+1,AL.age.max+1] <- (AL.Sx[AL.age.max+1])*AL.pred.red
  AL.popmat[1,] <- AL.pred.p.mm
  AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
}

AL.n.pred <- colSums(AL.n.mat)
plot(yrs, AL.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=AL.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 100
itdiv <- iter/10

AL.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
AL.s.arr <- AL.m.arr <- array(data=NA, dim=c(t+1, AL.age.max+1, iter))

for (e in 1:iter) {
  AL.popmat <- AL.popmat.orig
  
  AL.n.mat <- matrix(0, nrow=AL.age.max+1,ncol=(t+1))
  AL.n.mat[,1] <- AL.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    AL.s.alpha <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$alpha
    AL.s.beta <- estBetaParams(AL.Sx, AL.s.sd.vec^2)$beta
    AL.s.stoch <- rbeta(length(AL.s.alpha), AL.s.alpha, AL.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    AL.fert.stch <- rnorm(length(AL.popmat[,1]), AL.pred.p.mm, AL.m.sd.vec)
    AL.m.arr[i,,e] <- ifelse(AL.fert.stch < 0, 0, AL.fert.stch)
    
    AL.totN.i <- sum(AL.n.mat[,i], na.rm=T)
    AL.pred.red <- AL.a.lp/(1+(AL.totN.i/AL.b.lp)^AL.c.lp)
    
    diag(AL.popmat[2:(AL.age.max+1),]) <- (AL.s.stoch[-(AL.age.max+1)])*AL.pred.red
    AL.popmat[AL.age.max+1,AL.age.max+1] <- (AL.s.stoch[AL.age.max+1])*AL.pred.red
    AL.popmat[1,] <- AL.m.arr[i,,e]
    AL.n.mat[,i+1] <- AL.popmat %*% AL.n.mat[,i]
    
    AL.s.arr[i,,e] <- AL.s.stoch * AL.pred.red
    
  } # end i loop
  
  AL.n.sums.mat[e,] <- ((as.vector(colSums(AL.n.mat))/AL.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

AL.n.md <- apply(AL.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
AL.n.up <- apply(AL.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
AL.n.lo <- apply(AL.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,AL.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(AL.n.lo),1.05*max(AL.n.up)))
lines(yrs,AL.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,AL.n.up,lty=2,col="red",lwd=1.5)

AL.s.add <- AL.m.add  <- rep(0, AL.age.max+1)
for (m in 1:iter) {
  AL.s.add <- rbind(AL.s.add, AL.s.arr[ceiling(AL.gen.l):(t+1),,m])
  AL.m.add <- rbind(AL.m.add, AL.m.arr[ceiling(AL.gen.l):(t+1),,m])
}
AL.s.add <- AL.s.add[-1,]
AL.m.add <- AL.m.add[-1,]

AL.s.md <- apply(AL.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
AL.s.up <- apply(AL.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
AL.s.lo <- apply(AL.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(AL.age.vec,AL.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(AL.s.lo),1.05*max(AL.s.up)))
lines(AL.age.vec,AL.s.lo,lty=2,col="red",lwd=1.5)
lines(AL.age.vec,AL.s.up,lty=2,col="red",lwd=1.5)

AL.m.md <- apply(AL.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
AL.m.up <- apply(AL.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
AL.m.lo <- apply(AL.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(AL.age.vec,AL.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(AL.m.lo),1.05*max(AL.m.up)))
lines(AL.age.vec,AL.m.lo,lty=2,col="red",lwd=1.5)
lines(AL.age.vec,AL.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))

