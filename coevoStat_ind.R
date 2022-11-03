library(lme4)
library(splitstackshape)
library(car)
library(emmeans)
library(reshape2)

fi <- read.table("fertility_data.txt",header = TRUE)
ti <- read.table("longevity_data.txt",header = TRUE)
ti <- ti %>% filter(deathStage == "adult", Sex %in% c("Male","Female"))
m <- read.table("juvenile_mortality_data.txt", head = TRUE,sep = "\t")

m.s <- dcast(data = m,formula = mtDNA + nDNA + Temperature ~ Death.stage, value.var = 'longevity')[-7] ;names(m.s) <- c("mtDNA","nDNA","Temperature","n.Copepodid","n.Nauplius","n.Survived")m.s$cope.rate <- m.s$n.Copepodid/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived);m.s$naup.rate <- m.s$n.Nauplius/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived);m.s$larval.rate <- (m.s$n.Nauplius+m.s$n.Copepodid)/(m.s$n.Copepodid+m.s$n.Nauplius+m.s$n.Survived);m.s$naup.share <- m.s$naup.rate/m.s$larval.rate;m.s$cope.share <- m.s$cope.rate/m.s$larval.rate;m.s$mtDNA <- factor(m.s$mtDNA,levels = c("SD","LJ","BR","FHL"));m.s$nDNA <- factor(m.s$nDNA, levels = c("BR","SD"),labels = c("BR nDNA","SD nDNA"));m.s <- m.s %>% filter(Temperature %in% c("15ºC","25ºC"),nDNA %in% c("BR nDNA","SD nDNA"), mtDNA %in% c("SD","BR","LJ","FHL"));m.l <- melt(m.s, id.vars = c("Temperature","n.Copepodid","n.Nauplius","n.Survived","mtDNA","nDNA","surv.rate","larval.rate","naup.share","cope.share") );
m.l$variable <- factor(m.l$variable, levels =c("naup.rate","cope.rate"), labels = c("Nauplius","Copepodid"))
m.l$coevo <- 
  ifelse(as.character(m.l$mtDNA) == as.character(gsub("\ nDNA","",m.l$nDNA)), "matched","mismatched")
m.s$coevo <- 
  ifelse(as.character(m.s$mtDNA) == as.character(gsub("\ nDNA","",m.s$nDNA)), "matched","mismatched")
m$coevo <- ifelse(as.character(m$nDNA)==as.character(m$mtDNA),"matched","mismatched")

test <- c("Female fertility","Male fertility","Female longevity","Male longevity","Naupliar development rate","Copepodid development rate","Male Naupliar development rate","Female Naupliar development rate","Male Copepodid development rate","Female Copepodid development rate")

coevoStat <- NULL
rm(FHL_t1,SD_t1,LJ_t1,BR_t1)
for(i in test){
  print(paste("Modeling ",i,sep = ""))
  if(i == "Female fertility"){
    FHL_t1 <- fi %>% filter(Cyto == "FHL" | coevo == "matched", Sex == "Female")
    LJ_t1 <- fi %>% filter(Cyto == "LJ" | coevo == "matched", Sex == "Female")
    BR_t1 <- fi %>% filter(Cyto %in% c("BR","SD"), Nuclear != "BR", Sex == "Female")
    SD_t1 <- fi %>% filter(Cyto %in% c("BR","SD"), Nuclear != "SD" , Sex == "Female")
    mod <- lmer(data = LJ_t1, fert ~ coevo + Temperature + (1|Family:Clutch)+ (1|Nuclear))
    mod2 <- lmer(data = BR_t1, fert ~ coevo + Temperature  + (1|Family:Clutch))
    mod3 <- lmer(data = SD_t1, fert ~ coevo + Temperature  + (1|Family:Clutch))
    mod4 <- lmer(data = FHL_t1, fert ~ coevo + Temperature  + (1|Family:Clutch)+ (1|Nuclear))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          temp.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}
  if(i == "Male fertility"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- fi %>% filter(Cyto == "FHL" | coevo == "matched", Sex == "Male")
    LJ_t1 <- fi %>% filter(Cyto == "LJ" | coevo == "matched", Sex == "Male")
    BR_t1 <- fi %>% filter(Cyto %in% c("BR","SD"), Nuclear != "BR", Sex == "Male")
    SD_t1 <- fi %>% filter(Cyto %in% c("BR","SD"), Nuclear != "SD" , Sex == "Male")
    mod <- lmer(data = LJ_t1, fert ~ coevo + Temperature + (1|Family:Clutch)+ (1|Nuclear))
    mod2 <- lmer(data = BR_t1, fert ~ coevo + Temperature  + (1|Family:Clutch))
    mod3 <- lmer(data = SD_t1, fert ~ coevo + Temperature  + (1|Family:Clutch))
    mod4 <- lmer(data = FHL_t1, fert ~ coevo + Temperature  + (1|Family:Clutch)+ (1|Nuclear))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          temp.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}
  if(i == "Female longevity"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- ti %>% filter(Cyto == "FHL" | coevo == "match", Sex == "Female")
    LJ_t1 <- ti %>% filter(Cyto == "LJ" | coevo == "match", Sex == "Female")
    BR_t1 <- ti %>% filter(Cyto %in% c("BR","SD"), Nuclear != "BR", Sex == "Female")
    SD_t1 <- ti %>% filter(Cyto %in% c("BR","SD"), Nuclear != "SD" , Sex == "Female")
    mod <- lmer(data = LJ_t1, longevity ~ coevo + Temperature + (1|Family:Clutch)+ (1|Nuclear))
    mod2 <- lmer(data = BR_t1, longevity ~ coevo + Temperature  + (1|Family:Clutch))
    mod3 <- lmer(data = SD_t1, longevity ~ coevo + Temperature  + (1|Family:Clutch))
    mod4 <- lmer(data = FHL_t1, longevity ~ coevo + Temperature  + (1|Family:Clutch)+ (1|Nuclear))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          temp.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}  
  if(i == "Male longevity"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- ti %>% filter(Cyto == "FHL" | coevo == "match", Sex == "Male")
    LJ_t1 <- ti %>% filter(Cyto == "LJ" | coevo == "match", Sex == "Male")
    BR_t1 <- ti %>% filter(Cyto %in% c("BR","SD"), Nuclear != "BR", Sex == "Male")
    SD_t1 <- ti %>% filter(Cyto %in% c("BR","SD"), Nuclear != "SD" , Sex == "Male")
    mod <- lmer(data = LJ_t1, longevity ~ coevo + Temperature + (1|Family:Clutch)+ (1|Nuclear))
    mod2 <- lmer(data = BR_t1, longevity ~ coevo + Temperature  + (1|Family:Clutch))
    mod3 <- lmer(data = SD_t1, longevity ~ coevo + Temperature  + (1|Family:Clutch))
    mod4 <- lmer(data = FHL_t1, longevity ~ coevo + Temperature  + (1|Family:Clutch)+ (1|Nuclear))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          temp.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}  
  if(i == "Naupliar development rate"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- m %>% filter(mtDNA == "FHL" | coevo == "matched")
    LJ_t1 <- m %>% filter(mtDNA == "LJ" | coevo == "matched")
    BR_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "BR")
    SD_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "SD")
    mod <- lmer(data = LJ_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    mod2 <- lmer(data = BR_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID))
    mod3 <- lmer(data = SD_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID))
    mod4 <- lmer(data = FHL_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          temp.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}  
  if(i == "Copepodid development rate"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- m %>% filter(mtDNA == "FHL" | coevo == "matched")
    LJ_t1 <- m %>% filter(mtDNA == "LJ" | coevo == "matched")
    BR_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "BR")
    SD_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "SD")
    mod <- lmer(data = LJ_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    mod2 <- lmer(data = BR_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID))
    mod3 <- lmer(data = SD_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID))
    mod4 <- lmer(data = FHL_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          temp.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}   
  if(i == "Male Naupliar development rate"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- m %>% filter(mtDNA == "FHL" | coevo == "matched", Sex == "m")
    LJ_t1 <- m %>% filter(mtDNA == "LJ" | coevo == "matched", Sex == "m")
    BR_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "BR", Sex == "m")
    SD_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "SD", Sex == "m")
    mod <- lmer(data = LJ_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    mod2 <- lmer(data = BR_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID))
    mod3 <- lmer(data = SD_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID))
    mod4 <- lmer(data = FHL_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          temp.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)} 
  if(i == "Female Naupliar development rate"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- m %>% filter(mtDNA == "FHL" | coevo == "matched", Sex == "f")
    LJ_t1 <- m %>% filter(mtDNA == "LJ" | coevo == "matched", Sex == "f")
    BR_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "BR", Sex == "f")
    SD_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "SD", Sex == "f")
    mod <- lmer(data = LJ_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    mod2 <- lmer(data = BR_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID))
    mod3 <- lmer(data = SD_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID))
    mod4 <- lmer(data = FHL_t1, AgeAtCopepodite ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          temp.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}  
  if(i == "Male Copepodid development rate"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- m %>% filter(mtDNA == "FHL" | coevo == "matched", Sex == "m")
    LJ_t1 <- m %>% filter(mtDNA == "LJ" | coevo == "matched", Sex == "m")
    BR_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "BR", Sex == "m")
    SD_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "SD", Sex == "m")
    mod <- lmer(data = LJ_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    mod2 <- lmer(data = BR_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID))
    mod3 <- lmer(data = SD_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID))
    mod4 <- lmer(data = FHL_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          temp.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}  
  if(i == "Female Copepodid development rate"){rm(FHL_t1,SD_t1,LJ_t1,BR_t1);
    FHL_t1 <- m %>% filter(mtDNA == "FHL" | coevo == "matched", Sex == "f")
    LJ_t1 <- m %>% filter(mtDNA == "LJ" | coevo == "matched", Sex == "f")
    BR_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "BR", Sex == "f")
    SD_t1 <- m %>% filter(mtDNA %in% c("BR","SD"), nDNA != "SD", Sex == "f")
    mod <- lmer(data = LJ_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA))
    mod2 <- lmer(data = BR_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID))
    mod3 <- lmer(data = SD_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID))
    mod4 <- lmer(data = FHL_t1, copeDuration ~  Temperature + coevo + (1|Plate.ID) + (1|nDNA),control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    coevoStat  <- rbind(coevoStat, 
                        data.frame(
                          intercept.chi = c(unlist(Anova(mod, type = 3))[1],unlist(Anova(mod2, type = 3))[1],unlist(Anova(mod3, type = 3))[1],unlist(Anova(mod4, type = 3))[1]),
                          temp.chi = c(unlist(Anova(mod, type = 3))[2],unlist(Anova(mod2, type = 3))[2],unlist(Anova(mod3, type = 3))[2],unlist(Anova(mod4, type = 3))[2]),
                          coevo.chi = c(unlist(Anova(mod, type = 3))[3],unlist(Anova(mod2, type = 3))[3],unlist(Anova(mod3, type = 3))[3],unlist(Anova(mod4, type = 3))[3]),
                          intercept.p =c(unlist(Anova(mod, type = 3))[7],unlist(Anova(mod2, type = 3))[7],unlist(Anova(mod3, type = 3))[7],unlist(Anova(mod4, type = 3))[7]),
                          temp.p = c(unlist(Anova(mod, type = 3))[8],unlist(Anova(mod2, type = 3))[8],unlist(Anova(mod3, type = 3))[8],unlist(Anova(mod4, type = 3))[8]),
                          coevo.p = c(unlist(Anova(mod, type = 3))[9],unlist(Anova(mod2, type = 3))[9],unlist(Anova(mod3, type = 3))[9],unlist(Anova(mod4, type = 3))[9]),
                          mtDNA = c("LJ","BR","SD","FHL"),trait = i));rm(mod,mod2,mod3,mod4)}   
}




