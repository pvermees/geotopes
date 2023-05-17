# replace with your own working directory:
setwd('~/Documents/geotopes/R/')

intercepts <- function(dat){
    out <- rep(0,3)
    names(out) <- c('Ar36','Ar39','Ar40')
    out['Ar36'] <- lm(Ar36 ~ t, data=dat)$coefficients[1]
    out['Ar39'] <- lm(Ar39 ~ t, data=dat)$coefficients[1]
    out['Ar40'] <- lm(Ar40 ~ t, data=dat)$coefficients[1]
    return(out)
}

blankcorr <- function(sig,blk){
    out <- sig
    for (Ar in names(sig)){
        out[Ar] <- sig[Ar] - blk[Ar]
    }
    return(out)
}

process <- function(smpfile,stdfile,blkfile){
    std <- read.csv(stdfile,header=TRUE)
    smp <- read.csv(smpfile,header=TRUE)
    blk <- read.csv(blkfile,header=TRUE)

    stdint <- intercepts(std)
    smpint <- intercepts(smp)
    blkint <- intercepts(blk)

    stdblk <- blankcorr(stdint,blkint)
    smpblk <- blankcorr(smpint,blkint)

    # atmospheric correction
    a4036 <- 298.5
    std40r <- stdblk['Ar40'] - stdblk['Ar36']*a4036
    smp40r <- smpblk['Ar40'] - smpblk['Ar36']*a4036
    
    # calculate J
    tstd <- 27.8    # standard age, in Ma
    l40 <- 5.543e-4 # decay constant, in Ma-1
    J <- (exp(l40*tstd)-1)*(stdblk['Ar39']/std40r)

    # calculate age
    tx = log(1 + J*smp40r/smpblk['Ar39'])/l40
    out <- c(smpblk['Ar39'],tx)
    return(out)
}

setwd('~/Documents/geotopes/R/')
age <- process('smpl.csv','stnd.csv','blnk.csv')
