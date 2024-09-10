## Mock soil multifunctionality analysis using latent variable models

library(lavaan)
library(faux)



### Simulating correlated data
n = 100   # N soils
mu = 0     # mean of errors
sd = 10    # sd of errors
r = 0.4    # correlations among errors

# setup variables for soil functions
LVprod <- c(1:n)
LVbio <- 10 + -2 * 0.01 * LVprod
LVclim <- 0.5 * 0.01 * LVprod
LVwater <- 10 + 2*0.01* LVclim

# generate data
set.seed(08029185)
df1 <- 0.01*LVprod + rnorm_multi(n = n, 
                  mu = mu,
                  sd = sd,
                  r = r, 
                  varnames = c("wheat","festuca","trifolium","arabidopsis"
                             ),
                  empirical = FALSE)
df2 <- LVclim + rnorm_multi(n = n, 
                  mu = mu,
                  sd = sd,
                  r = r, 
                  varnames = c("control.substrate.co2","bean.co2","FYM.co2","sawdust.co2",
                               "control.substrate.n2o","bean.n2o","FYM.n2o","sawdust.n2o",
                               "control.substrate.ch4","bean.ch4","FYM.ch4","sawdust.ch4"
                               ),
                  empirical = FALSE)
df3 <- LVwater+rnorm_multi(n = n, 
                  mu = mu,
                  sd = sd,
                  r = r, 
                  varnames = c("wr","infil","field.cap","wilt",
                               "no3","nh4","po4","pb","cd","gly","fluo"
                               ),
                  empirical = FALSE)
df4 <- LVbio + rnorm_multi(n = n, 
                  mu = mu,
                  sd = sd,
                  r = r, 
                  varnames = c("psf1","psf2","psf3","psf4"),
                  empirical = FALSE)
df <-cbind(df1,df2,df3,df4)

pairs(df1)
pairs(df2[1:4])
pairs(df2[,5:12])
pairs(df3[,1:4])
pairs(df3[,5:11])
pairs(df4)

# separate measurement model per soil function
modLVpp <-"
LVpp =~ LVwheat + LVfestuca + LVtrifolium + LVarabidopsis
LVwheat =~ wheat
LVfestuca =~ festuca
LVtrifolium =~ trifolium
LVarabidopsis =~ arabidopsis

## variances
wheat~~108.0799*wheat              
festuca~~109.8652*festuca           
trifolium~~111.9223*trifolium      
arabidopsis~~107.0755*arabidopsis   
"

modLVclim <-"
LVcarbon=~ 1*LVcControl.co2 + LVcBean.co2 + LVcFYM.co2 + LVcSawdust.co2
LVcControl.co2 =~ control.substrate.co2
LVcBean.co2 =~ bean.co2
LVcFYM.co2 =~ FYM.co2
LVcSawdust.co2 =~ sawdust.co2

LVghg =~ 1*LVghgControl.n2o + LVghgBean.n2o + LVghgFYM.n2o + LVghgSawdust.n2o + LVghgControl.ch4 + LVghgBean.ch4 + LVghgFYM.ch4 + LVghgSawdust.ch4
LVghgControl.n2o =~ control.substrate.n2o
LVghgBean.n2o =~ bean.n2o
LVghgFYM.n2o =~ FYM.n2o
LVghgSawdust.n2o =~ sawdust.n2o
LVghgControl.ch4 =~ control.substrate.ch4
LVghgBean.ch4 =~ bean.ch4
LVghgFYM.ch4 =~ FYM.ch4
LVghgSawdust.ch4 =~ sawdust.ch4

## variances
control.substrate.co2 ~~ 107.0976*control.substrate.co2
bean.co2 ~~ 106.226*bean.co2
FYM.co2 ~~ 102.4883*FYM.co2 
sawdust.co2 ~~ 106.9629*sawdust.co2

control.substrate.n2o ~~ 107.698* control.substrate.n2o
bean.n2o ~~ 112.5011*bean.n2o
FYM.n2o ~~ 104.3462*FYM.n2o
sawdust.n2o ~~ 96.64822*sawdust.n2o
control.substrate.ch4 ~~ 106.7509*control.substrate.ch4
bean.ch4 ~~  107.9279*bean.ch4
FYM.ch4 ~~ 101.601*FYM.ch4
sawdust.ch4 ~~ 99.88639*sawdust.ch4
"

modLVwater <-"
LVwhc =~ 1*LVwr + LVinfil + LVfield.cap + LVwilt
LVwr =~ wr
LVinfil =~ infil
LVfield.cap =~ field.cap
LVwilt =~ wilt

LVpurif =~ 1*LVno3 + LVnh4 + LVpo4 + LVpb + LVcd + LVgly + LVfluo
LVno3 =~ no3
LVnh4 =~ nh4
LVpo4 =~ po4
LVpb =~ pb
LVcd =~ cd
LVgly =~ gly
LVfluo =~ fluo

## variances
wr ~~ 105.1202*wr
infil ~~ 101.6277*infil
field.cap ~~ 101.3991*field.cap
wilt ~~ 101.0982*wilt

no3 ~~ 99.05046*no3
nh4 ~~ 101.2444*nh4
po4 ~~ 103.8008*po4
pb ~~ 99.72934*pb
cd ~~ 110.2981*cd
gly ~~ 105.9757*gly
fluo ~~ 106.3622*fluo
"

modLVbio <-"
LVbio =~ 1*LVpsf1 + LVpsf2 + LVpsf3 + LVpsf4
LVpsf1 =~ psf1
LVpsf2 =~ psf2
LVpsf3 =~ psf3
LVpsf4 =~ psf4

## variances
psf1 ~~ 133.1146*psf1
psf2 ~~ 136.5257*psf2
psf3 ~~ 139.8648*psf3
psf4 ~~ 140.3208*psf4

"
mod.LVpp.fit1<-sem(model=modLVpp, data = df)
mod.LVclim.fit1<-sem(model=modLVclim, data = df)
mod.LVwater.fit1<-sem(model=modLVwater, data = df)
mod.LVbio.fit1<-sem(model=modLVbio, data = df)

summary(mod.LVpp.fit1,rsquare = F)
summary(mod.LVclim.fit1,rsquare = F)
summary(mod.LVwater.fit1,rsquare = F)
summary(mod.LVbio.fit1,rsquare = F)






