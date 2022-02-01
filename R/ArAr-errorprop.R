# load the data
std <- read.csv("stnd.csv",header=TRUE)
smp <- read.csv("smpl.csv",header=TRUE)
blk <- read.csv("blnk.csv",header=TRUE)

# plot the sample's Ar40 signal against time
plot(smp$t,smp$Ar40,type="p")
# plot the blank's Ar36 signal against time
plot(blk$t,blk$Ar36,type="p")

# linear regression of the data
std36fit <- lm(Ar36 ~ t, data=std)
std39fit <- lm(Ar39 ~ t, data=std)
std40fit <- lm(Ar40 ~ t, data=std)
smp36fit <- lm(Ar36 ~ t, data=smp)
smp39fit <- lm(Ar39 ~ t, data=smp)
smp40fit <- lm(Ar40 ~ t, data=smp)
blk36fit <- lm(Ar36 ~ t, data=blk)
blk39fit <- lm(Ar39 ~ t, data=blk)
blk40fit <- lm(Ar40 ~ t, data=blk)

# blank correction
std36b <- std36fit$coefficients[1] - blk36fit$coefficients[1]
std39b <- std39fit$coefficients[1] - blk39fit$coefficients[1]
std40b <- std40fit$coefficients[1] - blk40fit$coefficients[1]
smp36b <- smp36fit$coefficients[1] - blk36fit$coefficients[1]
smp39b <- smp39fit$coefficients[1] - blk39fit$coefficients[1]
smp40b <- smp40fit$coefficients[1] - blk40fit$coefficients[1]

# atmospheric correction
a4036 <- 298.5
std40r <- std40b - std36b*a4036
smp40r <- smp40b - smp36b*a4036

# calculate J
tstd <- 27.8    # standard age, in Ma
l40 <- 5.543e-4 # decay constant, in Ma-1
J <- (exp(l40*tstd)-1)*(std39b/std40r)

# calculate age
tx = log(1 + J*smp40r/smp39b)/l40

# calculate the variances of the intercepts:
std36merr2 <- vcov(std36fit)[1,1]
std39merr2 <- vcov(std39fit)[1,1]
std40merr2 <- vcov(std40fit)[1,1]
smp36merr2 <- vcov(smp36fit)[1,1]
smp39merr2 <- vcov(smp39fit)[1,1]
smp40merr2 <- vcov(smp40fit)[1,1]
blk36merr2 <- vcov(blk36fit)[1,1]
blk39merr2 <- vcov(blk39fit)[1,1]
blk40merr2 <- vcov(blk40fit)[1,1]

# variances of the blank corrections:
std36berr2 <- std36merr2 + blk36merr2
std39berr2 <- std39merr2 + blk39merr2
std40berr2 <- std40merr2 + blk40merr2
smp36berr2 <- smp36merr2 + blk36merr2
smp39berr2 <- smp39merr2 + blk39merr2
smp40berr2 <- smp40merr2 + blk40merr2

# variance of the atmospheric correction
std40rerr2 <- std40berr2 + std36berr2*a4036^2
smp40rerr2 <- smp40berr2 + smp36berr2*a4036^2

# J error
Jerr <- (exp(l40*tstd)-1) * (std39b/std40r) * 
  sqrt((std39berr2/std39b^2) + (std40rerr2/std40r^2))

# age error
txerr = (J*smp40r/smp39b)*
  sqrt((Jerr/J)^2 + (smp40rerr2/smp40r^2)+(smp39berr2/smp39b^2))/l40