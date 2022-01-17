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

# samp: a data frame of time-resolved mass spectrometer sample data
# stand: a data frame of time-resolved mass spectrometer standard data
# num: the column header of the numerator isotope
# den: the column header of the denominator isotope
# tstd: the age of the standard (in Ma)
# lambda: the decay constant of the denumerator isotope (in Myr-1)
# tb: the end of the blank (in seconds)
# ts: the start of the signal (in seconds)
# returns the atomic num/den ratio and its standard error
atomic <- function(samp,stand,num='Pb206',den='U238',tstd,lambda,tb,ts){
    bi <- samp$Time<tb
    si <- samp$Time>ts
    measured <- meanUPbRatio(D=samp[[num]],P=samp[[den]],bi=bi,si=si)
    cal <- calibration(D=std[[num]],P=std[[den]],tstd=tstd,
                       lambda=lambda,bi=bi,si=si)
    out <- c(0,0)
    out[1] <- cal[1]*measured[1]
    out[2] <- out[1] * sqrt( (cal[2]/cal[1])^2 + (measured[2]/measured[1])^2 )
    return(out)
}
