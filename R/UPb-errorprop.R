setwd('~/Documents/IsoplotR/')

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

# 10. propagate the standard errors of the ratio means from step 5:
se_smp0638m <- sd(smp0638)/sqrt(length(smp0638))
se_smp0735m <- sd(smp0735)/sqrt(length(smp0735))
se_std0638m <- sd(std0638)/sqrt(length(std0638))
se_std0735m <- sd(std0735)/sqrt(length(std0735))

# 11. standard errors of the atomic ratios from steps 7 & 8:
se_smp0638a <- smp0638a*sqrt( (se_smp0638m/smp0638m)^2 + (se_std0638m/std0638m)^2 )
se_smp0735a <- smp0735a*sqrt( (se_smp0735m/smp0735m)^2 + (se_std0735m/std0735m)^2 )

# 12. standard errors of the U-Pb ages:
x <- 1 + smp0638a
se_x <- se_smp0638a
se_t0638 <- (se_x/x)*(1/l38)
# repeating for Pb207/U238
y <- 1 + smp0735a
se_y <- se_smp0735a
se_t0735 <- (se_y/y)*(1/l35)