########## Data analyse Judith

########## Data preparation ###########

### set working directory (where your datafiles are)
setwd("~/My publications/Multifunctionality viewpoint/SOIL")

#### libraries
library(lavaan)
library(lme4)
library(openxlsx)
library(lattice)
library(arm)
library(vegan)

#### load data
d1 <-read.xlsx(sheet=4,"Supplementary information 3.xlsx")   # Exp 1 on 30 soils
d1 <- data.frame(Exp="Exp1",d1,Fertilized="Control")
d2 <-read.xlsx(sheet=5,"Supplementary information 3.xlsx")   # Exp 2 fertilization experiment
d2 <- data.frame(Exp="Exp2",d2)
names(d1)==names(d2)

### data clean up
d2$Fertilized <- factor(d2$Fertilized); levels(d2$Fertilized) <- c("Control","Control","High","Low","Low")
d1$Soil2 <- d1$Soil; d2$Soil2 <- factor(paste(d2$Soil,d2$Fertilized,sep="-"))
dat <- rbind(d1,d2)  ### adding the data from both experiments together
names(dat) <- c("exp","soil","Rep","pot","species","individu","mass_in_g","length_in_cm","age_day","fertilized","soil2")
dat$species <- factor(dat$species); levels(dat$species) <- c("arabidopsis","festuca","trifolium","trifolium","wheat")

dat$mass_in_g[dat$mass_in_g=="-"] <- NA       # record not observed values as NA (= not available)
dat$length_in_cm[dat$length_in_cm=="-"] <- NA
dat$mass_in_g <- as.numeric(dat$mass_in_g) #+ 0.0001
dat$soil3 <- factor(paste(dat$soil2,dat$Rep,sep="-"))

table(dat$species,is.na(dat$mass_in_g))  # not so many missing values


######### Basic information ############
table(dat$species,dat$age_day)  # harvest days
boxplot(mass_in_g~species,dat,log="y",ylab="Biomass (g/plant)")
bwplot(mass_in_g~factor(age_day)|species,dat,log=TRUE,ylab="Biomass (g/plant)")


######### RGR models ############
## After Daou & Shipley 2019

options(na.action="na.exclude")
#RGR were calculated as ln(biomass, mg) ~ age (fixed effect, days) + soil (random effect) 
### factor *1000 (see below) is to express the RGRs ihead() mg g-1 d-1 instead of g g-1 d-1
fit.wheat<-lmer(log(mass_in_g)~age_day+(age_day-1|soil3),data=dat,subset=(species=="wheat"))
RGR.wheat.by.soil<-coef(fit.wheat)$soil[,2]

fit.festuca<-lmer(log(mass_in_g)~age_day+(age_day-1|soil3),data=dat,subset=(species=="festuca"))
RGR.festuca.by.soil<-coef(fit.festuca)$soil[,2]

fit.trifolium<-lmer (log(mass_in_g)~age_day+(age_day-1|soil3),data=dat,subset= (species=="trifolium"))
RGR.trifolium.by.soil<-coef(fit.trifolium)$soil[,2]

fit.arabidopsis<-lmer(log(mass_in_g)~age_day+(age_day-1|soil3),data=dat,subset= (species=="arabidopsis"))
RGR.arabidopsis.by.soil<-coef(fit.arabidopsis)$soil[,2]


soils<-unique(dat$soil3)
RGR.all.species<-data.frame(soil=soils,wheat=RGR.wheat.by.soil*1000,festuca=RGR.festuca.by.soil*1000,trifolium=RGR.trifolium.by.soil*1000,arabidopsis=RGR.arabidopsis.by.soil*1000)
any(is.na(RGR.all.species))
cor(RGR.all.species[,2:5])  # correlations among variables
pairs(RGR.all.species[,2:5])

################# plot RGRs ###############
n.soils<-length(soils)
wheattransition <- data.frame(soils, RGR =  RGR.all.species$wheat, species = rep ("T.aestivum",n.soils))
trifoliumtransition <- data.frame(soils, RGR =  RGR.all.species$trifolium,species = rep ("T.pratense",n.soils))
arabidopsistransition <- data.frame(soils, RGR =  RGR.all.species$arabidopsis,species = rep ("A.thaliana",n.soils))
festucatransition <- data.frame(soils, RGR =  RGR.all.species$festuca,species = rep ("F.rubra",n.soils))

RGRall <- rbind (wheattransition,trifoliumtransition,arabidopsistransition,festucatransition)
RGRall$species <- factor(RGRall$species)
RGRall$species <- factor(RGRall$species,levels=c("F.rubra","T.pratense","T.aestivum","A.thaliana"))

# Boxplot of RGR of the four species
#RGRall$species=factor(RGRall$species , levels=levels(RGRall$species)[c(4,2,1,3)])
par(mfrow = (c(1,1)),xpd = F, mar=c(5,5.5,4,3))
boxplot(RGR~species,data=RGRall, main="RGR vs species", xlab="", ylab= bquote('RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'))


##### extract mean measurement error variances ###########

#Variances
SE.random.slopes <-function(x){
  # x is the output of a mixed model for ln(mass)~day
  library(arm)
  fixed.slope.se <- se.fixef(x)[2]
  random.slope.se <- se.ranef(x)$soil3
  se.sum <- sqrt(fixed.slope.se^2+random.slope.se^2)
  se.sum}
VarMeasErrorwheat <-  1000000*mean(SE.random.slopes(fit.wheat)^2)
VarMeasErrorfestuca <- 1000000*mean(SE.random.slopes(fit.festuca)^2)
VarMeasErrortrifolium <- 1000000*mean(SE.random.slopes(fit.trifolium)^2)
VarMeasErrorarabidopsis <-1000000*mean(SE.random.slopes(fit.arabidopsis)^2)

######## latent variable model for fertility ##############

# part 1 of the model. It uses the latent variable measurement model applied to the four RGRs.
# This part allows to predict the values of FG of the soils without considering the soil variables. 
# Thus, those predicted values would be comparable from a study to another.
# L is the latent variable we want to predict (generalized soil fertility; FG)
# Lwheat, Lfestuca, Ltrifolium and Larabidopsis represent the true responses of plants (RGRij' in Fig.1).
# wheat, festuca, trifolium and arabidopsis represent the measures of RGR (RGRij in Fig.1) considering the measurement errors (see below for calculation).
Modelall<-"
FG=~ 1*Lwheat+Lfestuca+Ltrifolium+Larabidopsis

Lwheat =~ wheat
Lfestuca =~ festuca
Ltrifolium =~ trifolium
Larabidopsis =~ arabidopsis

FG ~~ VarFG*FG
Lwheat ~~ Var1*Lwheat
Lfestuca ~~ Var2*Lfestuca
Ltrifolium ~~ Var3*Ltrifolium
Larabidopsis ~~ Var4*Larabidopsis

# error variances
wheat~~47.52132*wheat              # 24.74793 in D&S19
festuca~~118.9233 *festuca           # 20.67148 in D&S19
trifolium~~55.64609 *trifolium       # 16.54721 in D&S19
arabidopsis~~304.3423 *arabidopsis    # 26.71112 in D&S19
"
ModelMeasure.fit1<-sem(model=Modelall, data = RGR.all.species)   # note a negative variatiance is estimated, but it only means the scores of Lfestuca (FGfest) are in reverse order.... 
summary(ModelMeasure.fit1,rsquare = T)

# PCA of the RGR data
pca1 <- vegan::rda(RGR.all.species[,2:5])
plot(pca1, scaling = "symmetric", type="n") |>
  text("sites", cex=0.8) |>
  text("species", arrows=TRUE, length=0.02, col="red", cex=0.6)

modificationIndices(ModelMeasure.fit1)

pred = data.frame(RGR.all.species,predict(ModelMeasure.fit1))
colnames(pred)[c(1,6:10)]<- c ("soil.p", "FG" ,"FGwheat" ,"FGfest","FGtrif" ,"FGarab")

#pdf("Factor loadings.pdf",paper="a4")
par(mfrow=c(2,2),mar=c(4,6,1,1))
plot(wheat~FG,pred,ylab=bquote('Wheat RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FG index");abline(lm(wheat~FG,pred));text(-10,70,paste("R2 = ",round(summary(lm(wheat~FG,pred))$r.square,2)))#;plot(wheat~FGwheat,pred)
plot(festuca ~FG,pred,ylab=bquote('Festuca RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FG index");abline(lm(festuca~FG,pred));text(-10,60,paste("R2 = ",round(summary(lm(festuca~FG,pred))$r.square,2)))#;plot(festuca ~FGfest,pred)
plot(trifolium ~FG,pred,ylab=bquote('Trifolium RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FG index");abline(lm(trifolium~FG,pred));text(-10,90,paste("R2 = ",round(summary(lm(trifolium~FG,pred))$r.square,2)))#;plot(trifolium ~FGtrif,pred)
plot(arabidopsis~FG,pred,ylab=bquote('Arabidopsis RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FG index");abline(lm(arabidopsis~FG,pred));text(-10,150,paste("R2 = ",round(summary(lm(arabidopsis~FG,pred))$r.square,2)))#;plot(arabidopsis~FGarab,pred)
par(mfrow=c(1,1))
#dev.off()

par(mfrow=c(2,2),mar=c(4,6,1,1))
plot(wheat~FGwheat,pred,ylab=bquote('Wheat RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FGwheat index");abline(lm(wheat~FGwheat,pred));text(-10,70,paste("R2 = ",round(summary(lm(wheat~FG,pred))$r.square,2)))#;plot(wheat~FGwheat,pred)
plot(festuca ~FGfest,pred,ylab=bquote('Festuca RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FGfest index");abline(lm(festuca~FGfest,pred));text(-10,60,paste("R2 = ",round(summary(lm(festuca~FG,pred))$r.square,2)))#;plot(festuca ~FGfest,pred)
plot(trifolium ~FGtrif,pred,ylab=bquote('Trifolium RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FGtrif index");abline(lm(trifolium~FGtrif,pred));text(-10,90,paste("R2 = ",round(summary(lm(trifolium~FG,pred))$r.square,2)))#;plot(trifolium ~FGtrif,pred)
plot(arabidopsis~FGarab,pred,ylab=bquote('Arabidopsis RGR ('* 'mg' ~ 'g'^-1~'d'^-1*')'),xlab="FGarab index");abline(lm(arabidopsis~FGarab,pred));text(-10,150,paste("R2 = ",round(summary(lm(arabidopsis~FG,pred))$r.square,2)))#;plot(arabidopsis~FGarab,pred)
par(mfrow=c(1,1))


#### Repeatablity

pred$rep <- c(rep(c("A","B"),41)) 

plot(pred$FG[pred$rep=="B"]~pred$FG[pred$rep=="A"],ylab="FG in rep B",xlab="FG in rep A");abline(0,1)



##### Figure responsiveness fertilization

pred.dat <- data.frame(cbind(predict(ModelMesure.fit1),factor("soil_rep" <- c(1:82))))

colnames(pred.dat)[6] <- "soil_rep"
boxplot(FG~soil_rep,pred.dat)

pred.dat$Group <- c(1:82) #factor(RGR.all.species$soil)
pred.dat$Group[c(3:62)] <- "Samples"
pred.dat$Group[c(63:82)] <- "Exp2"
pred.dat$Group[c(1:2,75:76)] <- "Potting soil"
pred.dat$Group[c(77:78)] <- "Sand"

pred.dat$Exp[c(1:82)] <- "Exp1"
pred.dat$Exp[c(63:82)] <- "Exp2"

pred.dat$Fertilized <- c(1:82)
pred.dat$Fertilized[c(63:64,69:70,75:78)] <- "Control"
pred.dat$Fertilized[c(67:68,73,74,81:82)] <- "Low"
pred.dat$Fertilized[c(65:66,71:72,79:80)] <- "High"

pred.dat$soil2 <- c(1:82)
pred.dat$soil2[c(63:68)] <- "Clay"
pred.dat$soil2[c(69:74)] <- "Mix"
pred.dat$soil2[c(77:82)] <- "Sand"
pred.dat$soil2[c(75:76)] <- "Potting soil"

pred.dat$Group <- factor(pred.dat$Group, levels=c('Sand', 'Samples', 'Potting soil'))
pred.dat$Fertilized <- factor(pred.dat$Fertilized, levels=c('Control', 'Low', 'High'))

boxplot(FG~Group,pred.dat[pred.dat$Group!="Exp2",],ylab="Generalized soil fertility (Fg)",xlab="")

dotplot(FG~Fertilized| soil2,data=pred.dat[pred.dat$Exp=="Exp2" & pred.dat$soil2!="Potting soil",],ylab="FG index",xlab="Fertilization")

anova(lm(FG~Fertilized*soil2,data=pred.dat[pred.dat$Exp=="Exp2" & pred.dat$soil2!="Potting soil",]))


