library(lme4)
library(splitstackshape)
library(car)

fi <- read.table("fertility_data.txt",header = TRUE)
ti <- read.table("longevity_data.txt",header = TRUE)
ti <- ti %>% filter(deathStage == "adult", Sex %in% c("Male","Female"))

fi.l <- split.data.frame(x = fi , interaction(fi$Temperature,fi$Nuclear,fi$Sex)); names(fi.l) <- paste(names(fi.l),"fertility",sep =".");
ti.l <- split.data.frame(x = ti , interaction(ti$Temperature,ti$Nuclear,ti$Sex)); names(ti.l) <- paste(names(ti.l),"longevity",sep =".");

coevoStat <- NULL
if(length(fi.l)==length(ti.l)){n <- length(ti.l)}
for(i in seq(1:n)){
  m1 <- lmer(data = fi.l[[i]], formula = fert ~ coevo + (1|Family:Clutch));  # Fertility model
  m2 <- lmer(data = ti.l[[i]], formula = longevity ~ coevo + (1|Family:Clutch)); # Longevity model
  coevoStat  <- rbind(coevoStat, #Create table from type III ANOVA
               data.frame(chi.squared=Anova(m1, type =3)$Chisq[2],df=Anova(m1, type =3)$Df[2],p.value=Anova(m1, type =3)$`Pr(>Chisq)`[2],row.names = names(fi.l[i])),
               data.frame(chi.squared=Anova(m2, type =3)$Chisq[2],df=Anova(m2, type =3)$Df[2],p.value=Anova(m2, type =3)$`Pr(>Chisq)`[2],row.names = names(ti.l[i])))}
coevoStat$ID <- row.names(coevoStat);coevoStat <- cSplit(coevoStat,"ID", sep = ".");names(coevoStat)[4:7]<- c("Temperature","Nuclear","Sex","Trait")

coevoStat[order(Trait,Temperature,Nuclear,Sex),]
