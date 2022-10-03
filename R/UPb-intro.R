setwd('~/Documents/geotopes/R/')

# 1. load the data
smp <- read.csv("91500.csv",header=TRUE)
std <- read.csv("GJ1.csv",header=TRUE)

# 2. plot the 238U signal against time
plot(smp$Time,smp$U238,type='b')

# 3. blank correction
bi <- (smp$Time<20) # blank indices
si <- (smp$Time>25) # signal indices
smp06b <- smp$Pb206[si] - mean(smp$Pb206[bi])
smp07b <- smp$Pb207[si] - mean(smp$Pb207[bi])
smp38b <- smp$U238[si] - mean(smp$U238[bi])
smp35b <- smp$U235[si] - mean(smp$U235[bi])
std06b <- std$Pb206[si] - mean(std$Pb206[bi])
std07b <- std$Pb207[si] - mean(std$Pb207[bi])
std38b <- std$U238[si] - mean(std$U238[bi])
std35b <- std$U235[si] - mean(std$U235[bi])

# 4. ratios
smp0638 <- smp06b/smp38b
smp0735 <- smp07b/smp35b
std0638 <- std06b/std38b
std0735 <- std07b/std35b

# 5. ratio means
smp0638m <- mean(smp0638)
smp0735m <- mean(smp0735)
std0638m <- mean(std0638)
std0735m <- mean(std0735)

# 6. expected ratios for the standard
tstd <- 600.4    # age of the standard in Ma
l38 <- 1.55125e-4 # 238U decay constant in Ma-1
l35 <- 9.8485e-4 # 235U decay constant in Ma-1
std0638a <- exp(l38*tstd)-1 # 206Pb/238U ingrowth equation
std0735a <- exp(l35*tstd)-1 # 207Pb/235U ingrowth equation

# 7. correction factors
c0638 <- std0638a/std0638m
c0735 <- std0735a/std0735m

# 8. sample-standard bracketing of the sample
smp0638a <- c0638*smp0638m
smp0735a <- c0735*smp0735m

# 9. calculate the age
t0638 <- log(1+smp0638a)/l38
t0735 <- log(1+smp0735a)/l35

# 10. Wetherill concordia
# A. basic solution
tt <- seq(from=0,to=2000,length.out=50)
x <- exp(l35*tt)-1
y <- exp(l38*tt)-1
plot(x,y,type='l')
points(smp0735a,smp0638a)
# B. fancy solution
tt <- seq(from=0,to=2000,length.out=50)
x <- exp(l35*tt)-1
y <- exp(l38*tt)-1
plot(x,y,type='l',
     xlab=expression(''^207*'Pb/'^235*'U'),
     ylab=expression(''^206*'Pb/'^238*'U'))
points(smp0735a,smp0638a,pch=21,bg='white')
# add age ticks to the concordia plot:
ttick <- seq(from=0,to=2000,by=250)
xt <- exp(l35*ttick)-1
yt <- exp(l38*ttick)-1
points(xt,yt,pch=16)
text(xt,yt,labels=ttick,pos=2)
