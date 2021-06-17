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

# calculate the standard errors of the intercepts:
std36merr <- summary(std36fit)$coefficients[1,2]
std39merr <- summary(std39fit)$coefficients[1,2]
std40merr <- summary(std40fit)$coefficients[1,2]
smp36merr <- summary(smp36fit)$coefficients[1,2]
smp39merr <- summary(smp39fit)$coefficients[1,2]
smp40merr <- summary(smp40fit)$coefficients[1,2]
blk36merr <- summary(blk36fit)$coefficients[1,2]
blk39merr <- summary(blk39fit)$coefficients[1,2]
blk40merr <- summary(blk40fit)$coefficients[1,2]

# standard errors of the blank corrections:
std36berr <- std36merr + blk36merr
std39berr <- std39merr + blk39merr
std40berr <- std40merr + blk40merr
smp36berr <- smp36merr + blk36merr
smp39berr <- smp39merr + blk39merr
smp40berr <- smp40merr + blk40merr

# standard error of the atmospheric correction
std40rerr <- std40berr + std36berr*a4036
smp40rerr <- smp40berr + smp36berr*a4036

# J error
Jerr <- (exp(l40*tstd)-1) * (std39b/std40r) * 
  sqrt((std39berr/std39b)^2+(std40rerr/std40r)^2)

# age error
txerr = (J*smp40r/smp39b)*
  sqrt((Jerr/J)^2 + (smp40rerr/smp40r)^2+(smp39berr/smp39b)^2)/l40
