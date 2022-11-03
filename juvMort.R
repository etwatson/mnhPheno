library(reshape2)
library(dplyr)
library(lme4)
library(MuMIn)
library(car)
setwd("~/Dropbox/Cybrids (1)/manuscript/original_code/")
m <- read.table("juvenile_mortality_data.txt", head = TRUE,sep = "\t")
m$Temperature <- factor(m$Temperature)
m$Sex <- factor(m$Sex)
m$mtDNA <- factor(m$mtDNA)
m$nDNA <- factor(m$nDNA)
m$Death.stage <- factor(m$Death.stage)
m$Family <- factor(m$Family)

#count number of individuals for each death stage (nauplius, copepodid, adult)
m.s <- dcast(data = m,formula = mtDNA + nDNA + Temperature ~ Death.stage, value.var = 'longevity')[-7] ;names(m.s) <- c("mtDNA","nDNA","Temperature","n.Copepodid","n.Nauplius","n.Survived")
# Calculate mortality rates for each death stage
m.s$surv.rate <- m.s$n.Survived/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived) 
m.s$cope.rate <- m.s$n.Copepodid/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived)
m.s$naup.rate <- m.s$n.Nauplius/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived)
m.s$larval.rate <- (m.s$n.Nauplius+m.s$n.Copepodid)/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived)
m.s$naup.share <- m.s$naup.rate/m.s$larval.rate
m.s$cope.share <- m.s$cope.rate/m.s$larval.rate
m.s$mtDNA <- factor(m.s$mtDNA,levels = c("SD","LJ","BR","FHL"))
m.s$nDNA <- factor(m.s$nDNA, levels = c("BR","SD"),labels = c("BR nDNA","SD nDNA"))
m.s <- m.s %>% filter(Temperature %in% c("15ºC","25ºC"),nDNA %in% c("BR nDNA","SD nDNA"), mtDNA %in% c("SD","BR","LJ","FHL")) # ensure no NAs
m.l <- melt(m.s, id.vars = c("Temperature","n.Copepodid","n.Nauplius","n.Survived",
                             "mtDNA","nDNA","surv.rate","larval.rate","naup.share","cope.share") )
m.l$variable <- factor(m.l$variable, levels =c("naup.rate","cope.rate"), labels = c("Nauplius","Copepodid"))

# modeling developmental stages
a <- lmer(copeDuration ~ Temperature + mtDNA +nDNA + (1|Plate.ID),data = m %>% filter(copeDuration > 0))
b <- lmer(copeDuration ~ Temperature + mtDNA *nDNA + (1|Plate.ID),data = m %>% filter(copeDuration > 0))
c <- lmer(copeDuration ~ Temperature * nDNA +mtDNA + (1|Plate.ID),data = m %>% filter(copeDuration > 0))
d <- lmer(copeDuration ~ Temperature * nDNA *mtDNA + (1|Plate.ID),data = m %>% filter(copeDuration > 0))
model.sel(a,b,c,d)
mod.cd <- d
Anova(mod.cd, type=3) 

a <- lmer(AgeAtCopepodite ~ Temperature + mtDNA +nDNA + (1|Plate.ID),data = m)
b <- lmer(AgeAtCopepodite ~ Temperature + mtDNA *nDNA + (1|Plate.ID),data = m)
c <- lmer(AgeAtCopepodite ~ Temperature * nDNA +mtDNA + (1|Plate.ID),data = m)
d <- lmer(AgeAtCopepodite ~ Temperature * nDNA *mtDNA + (1|Plate.ID),data = m)
model.sel(a,b,c,d)
mod.ct <- d
Anova(mod.ct, type=3) 

a <- lmer(AgeAtMaturity ~ Temperature + mtDNA +nDNA + (1|Plate.ID),data = m)
b <- lmer(AgeAtMaturity ~ Temperature + mtDNA *nDNA + (1|Plate.ID),data = m)
c <- lmer(AgeAtMaturity ~ Temperature * nDNA +mtDNA + (1|Plate.ID),data = m)
d <- lmer(AgeAtMaturity ~ Temperature * nDNA *mtDNA + (1|Plate.ID),data = m)
model.sel(a,b,c,d)
mod.m <- d
Anova(mod.m, type=3) 
