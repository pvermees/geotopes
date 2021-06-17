setwd('~/Documents/IsoplotR/')

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

tUPb <- function(samp,stand,num='Pb206',den='U238',tstd,lambda,tb,ts){
    measured <- meanUPbRatio(D=samp[[num]],P=samp[[den]],
                             bi=(samp$Time<tb),si=(samp$Time>ts))
    cal <- calibration(D=std[[num]],P=std[[den]],
                       tstd=tstd,lambda=lambda,
                       bi=(stand$Time<tb),si=(stand$Time>ts))
    atomic <- cal[1]*measured[1]
    atomic_err <- atomic * sqrt( (cal[2]/cal[1])^2 + (measured[2]/measured[1])^2 )
    out <- c(0,0)
    out[1] <- log(1+atomic)/lambda
    out[2] <- (atomic_err/(1+atomic))*(1/lambda)
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
