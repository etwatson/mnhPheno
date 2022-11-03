library(lme4)
library(emmeans)
library(dplyr)
library(MuMIn)

setwd("~/Dropbox/Cybrids (1)/manuscript/original_code/")
ti <- read.table("longevity_data.txt",header = TRUE)
fi <- read.table("fertility_data.txt",header = TRUE)
ti <- ti %>% filter(deathStage == "adult", Sex %in% c("Male","Female"))
ti.l <- split(ti, ti$Temperature)

lon_sexEff <- NULL
for (i in 1:2){
  T <- ifelse(i == 1,"15ºC","25ºC")
  #Model 1: Sex x mtDNA x nDNA
  lon_sexEff[[length(lon_sexEff) + 1]] <- lmer(data=ti.l[[i]], longevity ~ Sex * Cyto * Nuclear + (1|Family:Clutch)); names(lon_sexEff)[length(lon_sexEff)] <- paste('SxMxN',T,"lon", sep = ".");
  #Model 2: Sex x mtDNA + nDNA
  lon_sexEff[[length(lon_sexEff) + 1]] <- lmer(data=ti.l[[i]], longevity ~ Sex * Cyto + Nuclear + (1|Family:Clutch));names(lon_sexEff)[length(lon_sexEff)] <- paste('SxM+N',T, "lon",sep = ".");
  #Model 3: Sex x nDNA + mtDNA
  lon_sexEff[[length(lon_sexEff) + 1]] <- lmer(data=ti.l[[i]], longevity ~ Sex * Nuclear + Cyto + (1|Family:Clutch));names(lon_sexEff)[length(lon_sexEff)] <- paste('SxN+M',T, "lon",sep = ".");
  #Model 4: Sex + mtDNA + nDNA
  lon_sexEff[[length(lon_sexEff) + 1]] <- lmer(data=ti.l[[i]], longevity ~ Sex + Cyto + Nuclear + (1|Family:Clutch));names(lon_sexEff)[length(lon_sexEff)] <- paste('S+M+N',T, "lon",sep = ".");
}

fi.l <- split(fi, fi$Temperature)
fert_sexEff <- NULL
for (i in 1:2){
  T <- ifelse(i == 1,"15ºC","25ºC")
  #Model 1: Sex x mtDNA x nDNA
  fert_sexEff[[length(fert_sexEff) + 1]] <- lmer(data=fi.l[[i]], adj.fert ~ Sex * Cyto * Nuclear + (1|Family:Clutch));names(fert_sexEff)[length(fert_sexEff)] <- paste('SxMxN',T, "fert",sep = ".");
  #Model 2: Sex x mtDNA + nDNA
  fert_sexEff[[length(fert_sexEff) + 1]] <- lmer(data=fi.l[[i]], adj.fert ~ Sex * Cyto + Nuclear + (1|Family:Clutch));names(fert_sexEff)[length(fert_sexEff)] <- paste('SxM+N',T, "fert",sep = ".");
  #Model 3: Sex x nDNA + mtDNA
  fert_sexEff[[length(fert_sexEff) + 1]] <- lmer(data=fi.l[[i]], adj.fert ~ Sex * Nuclear + Cyto + (1|Family:Clutch));names(fert_sexEff)[length(fert_sexEff)] <- paste('SxN+M',T, "fert",sep = ".");
  #Model 4: Sex + mtDNA + nDNA
  fert_sexEff[[length(fert_sexEff) + 1]] <- lmer(data=fi.l[[i]], adj.fert ~ Sex + Cyto + Nuclear + (1|Family:Clutch));names(fert_sexEff)[length(fert_sexEff)] <- paste('S+M+N',T, "fert",sep = ".");
}

sexEff <- c(lon_sexEff,fert_sexEff)
model.sel(sexEff[[1]],sexEff[[2]],sexEff[[3]],sexEff[[4]]) # 15C longevity
model.sel(sexEff[[5]],sexEff[[6]],sexEff[[7]],sexEff[[8]]) # 25C longevity
model.sel(sexEff[[9]],sexEff[[10]],sexEff[[11]],sexEff[[12]])# 15C Fertility
model.sel(sexEff[[13]],sexEff[[14]],sexEff[[15]],sexEff[[16]]) # 25C Fertility

sexEffC <- NULL
for (j in 1:length(sexEff)){
  m <- strsplit(names(sexEff[j]),"\\.")[[1]][1];t <- strsplit(names(sexEff[j]),"\\.")[[1]][2];ft <- strsplit(names(sexEff[j]),"\\.")[[1]][3];
  # get contrasts between the sexes using emmeans 
  if (m == "S+M+N"){df <- data.frame(plot(pairs(emmeans(sexEff[j][[1]],1 ~ Sex)),plotit=FALSE),Temperature = t,term = "Sex",Cyto = NA, Nuclear = NA,Trait = ft)}
    else if (m == "SxM+N"){df <- data.frame(plot(pairs(emmeans(sexEff[j][[1]],1 ~ Sex|Cyto)),plotit=FALSE),Temperature = t,term = "Sex x mtDNA", Nuclear = NA,Trait = ft)}
    else if (m == "SxN+M"){df <- data.frame(plot(pairs(emmeans(sexEff[j][[1]],1 ~ Sex|Nuclear)),plotit=FALSE),Temperature = t,term = "Sex x nDNA", Cyto = NA,Trait = ft)}
    else {df <- data.frame(plot(pairs(emmeans(sexEff[j][[1]],1 ~ Sex|Cyto + Nuclear)),plotit=FALSE),Temperature = t,term = "Sex x mtDNA x nDNA",Trait = ft)}
  sexEffC <- rbind(sexEffC,df)
}

# significant contrasts have 95% CI that do not pass the origin, negative = male biased, positive = female biased 
sexEffC$Nuclear <- factor(sexEffC$Nuclear, levels = c("BR","SD"));sexEffC$Cyto <- factor(sexEffC$Cyto, levels = c("SD","BR","LJ","FHL"));sexEffC$bias <- ifelse(sexEffC$lower.CL < 0 & sexEffC$upper.CL < 0, "male","NS");sexEffC$bias <- ifelse(sexEffC$lower.CL > 0 & sexEffC$upper.CL > 0, "female",sexEffC$bias);

