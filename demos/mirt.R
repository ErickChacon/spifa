library(mirt)
#make some data
set.seed(1234)
N <- 750
a <- matrix(rlnorm(10,.3,1),10,1)
d <- matrix(rnorm(10), 10)
Theta <- matrix(sort(rnorm(N)))
pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
pseudoIQ <- (pseudoIQ - mean(pseudoIQ))/10  #rescale variable for numerical stability
group <- factor(rep(c('G1','G2','G3'), each = N/3))
data <- simdata(a,d,N, itemtype = rep('2PL',10), Theta=Theta)

covdata <- data.frame(group, pseudoIQ)
#use parallel computing
mirtCluster()

#specify IRT model
model <- 'Theta = 1-10'

#model with no person predictors
mod0 <- mirt(data, model, itemtype = 'Rasch')

#group as a fixed effect predictor (aka, uniform dif)
mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items)
mod1 <- mixedmirt(data, covdata, model = 1, fixed = ~ 0 + group + items)
anova(mod0, mod1)
summary(mod1)
coef(mod1)

#
