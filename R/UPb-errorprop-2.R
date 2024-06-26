setwd('~/Documents/geotopes/R/')

blankcorr <- function(dat,bi,si){
    return(dat[si] - mean(dat[bi]))
}

meanUPbRatio <- function(D,P,bi,si){
    ratiovec <- blankcorr(D,bi,si)/blankcorr(P,bi,si)
    out <- c(0,0)
    out[1] <- mean(ratiovec)
    out[2] <- sd(ratiovec)/sqrt(length(ratiovec))
    return(out)
}

calibration <- function(D,P,tstd,lambda,bi,si){
    atomic <- exp(lambda*tstd)-1
    measured <- meanUPbRatio(D,P,bi,si)
    out <- c(0,0)
    out[1] <- atomic/measured[1]
    out[2] <- out[1]*measured[2]/measured[1]
    return(out)
}

atomic <- function(samp,stand,num='Pb206',den='U238',tstd,lambda,tb,ts){
    bi <- samp$Time<tb
    si <- samp$Time>ts
    measured <- meanUPbRatio(D=samp[[num]],P=samp[[den]],bi=bi,si=si)
    cal <- calibration(D=std[[num]],P=std[[den]],
                       tstd=tstd,lambda=lambda,bi=bi,si=si)
    out <- c(0,0)
    out[1] <- cal[1]*measured[1]
    out[2] <- out[1]*sqrt( (cal[2]/cal[1])^2 + (measured[2]/measured[1])^2 )
    return(out)
}

tUPb <- function(samp,stand,num='Pb206',den='U238',tstd,lambda,tb,ts){
    a <- atomic(samp,stand,num,den,tstd,lambda,tb,ts)
    out <- c(0,0)
    out[1] <- log(1+a[1])/lambda
    out[2] <- (a[2]/(1+a[1]))*(1/lambda)
    return(out)
}

tb <- 20 # end of blank (in seconds)
ts <- 25 # start of signal (in seconds)
tstd <- 600.4 # Ma
l38 <- 1.55125e-4 # Myr-1
l35 <- 9.8485e-4 # Myr-1
l32 <- 0.495e-4 # Myr-1

smp <- read.csv("91500.csv",header=TRUE)
std <- read.csv("GJ1.csv",header=TRUE)
t0638 <- tUPb(samp=smp,stand=std,num='Pb206',den='U238',
              tstd=tstd,lambda=l38,tb=tb,ts=ts)
t0735 <- tUPb(samp=smp,stand=std,num='Pb207',den='U235',
              tstd=tstd,lambda=l35,tb=tb,ts=ts)
t0832 <- tUPb(samp=smp,stand=std,num='Pb208',den='Th232',
              tstd=tstd,lambda=l32,tb=tb,ts=ts)
