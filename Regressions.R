### Variable switching in regression

# simulate data
Age <- rep(seq(10,50,by=10),each=10)
Exersize <- 2* Age + rnorm(length(Age),0,10)
Cholesterol <- 100 * Age + -15 * Exersize + rnorm(length(Age),0,1)
x3 <- rnorm(length(Age),10,1)  # unrelated x3 is needed to have enough d.f. to test the SEM model below
d <- data.frame(Age,Exersize,Cholesterol,x3)

# plotting
par(mfrow=c(2,2))
plot(Cholesterol~Exersize,d)
plot(Cholesterol~Age,d)
plot(Exersize~Age,d)
par(mfrow=c(1,1))

## ordinary regressions
summary(lm(Cholesterol~Exersize+Age))   # correct model and parameter estimates
summary(lm(Cholesterol~Exersize))       # Exersize increases cholesterol!


## SEM approach
library(lavaan)
mod.good <- "
Cholesterol ~ Exersize+Age
Exersize ~ Age
x3 ~1   # unrelated x3 is needed to have enough d.f. to test the model
"
fit.good <- sem(mod.good,data=d)
summary(fit.good)   # Correct results

mod.bad <- "
Cholesterol ~ Exersize
Exersize ~ Age
x3 ~1   # unrelated x3 is needed to have enough d.f. to test the model
"
fit.bad <- sem(mod.bad,data=d)
summary(fit.bad)   # SEM can lead to the same wrong parameter estimates if the model is misspecified, luckily it does tell you the model does not fit the data well
modificationIndices(fit.bad)  # Inspection of the model will tell you to include Age as a predictor of Cholesterol (high mi and high epc)

