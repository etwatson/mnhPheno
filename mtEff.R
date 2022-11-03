library(dplyr)
library(lme4)
library(MuMIn)
library(splitstackshape)
library(car)


fi <- read.table("fertility_data.txt",header = TRUE)
ti <- read.table("longevity_data.txt",header = TRUE)

# Male fertility data
rM <- fi %>% filter(Sex == "Male")
a<- lmer(data = rM,  #set fully parameterized model
         fert ~ Temperature * Cyto * Nuclear  + (1 | Family:Clutch),
         na.action = "na.fail" ,REML = FALSE)
dredge(a,rank = "AICc",evaluate = TRUE,fixed = ~ Temperature + Cyto + Nuclear) #  perform model testing to find best-fit, most reduced model
fert.M.lmer <- a
# Get effect size of each model term
rMeff1 <- plot(eff_size(emmeans(fert.M.lmer, 1 ~ Temperature + Nuclear + Cyto), sigma = sqrt(mean(sigma(fert.M.lmer)^2)), edf = df.residual(fert.M.lmer)),plotit = FALSE); rMeff1 <- rMeff1[-length(rMeff1)];

# split the Contrast columns and name the "contrast types" in rMeff1 
rMeff1<- cSplit(indt =rMeff1,splitCols = "contrast" ,sep = " "); rMeff1$contrast_2 <- as.character(rMeff1$contrast_2); rMeff1$contrast_3 <- as.character(rMeff1$contrast_3); rMeff1$contrast_6 <- as.character(rMeff1$contrast_6); rMeff1$contrast_7 <- as.character(rMeff1$contrast_7);
rMeff1$comp <- ifelse(rMeff1$contrast_1 != rMeff1$contrast_5 &rMeff1$contrast_2 == rMeff1$contrast_3 &rMeff1$contrast_6 == rMeff1$contrast_7 &rMeff1$contrast_2 == rMeff1$contrast_6,"Temperature","foobar");
rMeff1$comp <- ifelse(rMeff1$contrast_1 == rMeff1$contrast_5 &rMeff1$contrast_2 != rMeff1$contrast_6 &rMeff1$contrast_3 == rMeff1$contrast_7,"Nuclear",rMeff1$comp);
rMeff1$comp <- ifelse(rMeff1$contrast_1 != rMeff1$contrast_5 &rMeff1$contrast_2 == rMeff1$contrast_6 &rMeff1$contrast_3 == rMeff1$contrast_7,"Temperature",rMeff1$comp);
rMeff1$comp <- ifelse(rMeff1$contrast_1 == rMeff1$contrast_5 & rMeff1$contrast_2 == rMeff1$contrast_6 & rMeff1$contrast_3 != rMeff1$contrast_7, "Mito",rMeff1$comp); rMeff1$contrast_4 <- NULL;
rMeff1$Sex <- "Male"; rMeff1$Trait <- "Fertility"

#Female fertility data 
hN <- fi %>% filter(Sex == "Female")
a<- lmer(data = hN, 
         fert ~ Temperature * Cyto * Nuclear  + (1 | Family:Clutch),
         na.action = "na.fail" ,REML = FALSE)
dredge(a,rank = "AICc",evaluate = TRUE,fixed = ~ Temperature + Cyto + Nuclear)
# Model testing did not fully differentiate top two models
a1 <- lmer(data = hN, 
                    fert ~ Temperature + Cyto + Nuclear  + (1 | Family:Clutch),
                    na.action = "na.fail" ,REML = FALSE)
a2 <- lmer(data = hN, 
          fert ~ Temperature * Nuclear + Cyto + (1 | Family:Clutch),
          na.action = "na.fail" ,REML = FALSE)
model.sel(a1,a2) # adding extra term does not drastically improve fit, keeping more reduced model
fert.F.lmer <- a1
hNeff1 <- plot(eff_size(emmeans(fert.F.lmer, 1 ~ Temperature + Nuclear + Cyto), sigma = sqrt(mean(sigma(fert.F.lmer)^2)), edf = df.residual(fert.F.lmer)),plotit = FALSE); hNeff1 <- hNeff1[-length(hNeff1)];

# now split the Contrast columns and name the "contrast types" in hNeff1 
hNeff1<- cSplit(indt =hNeff1,splitCols = "contrast" ,sep = " "); hNeff1$contrast_2 <- as.character(hNeff1$contrast_2); hNeff1$contrast_3 <- as.character(hNeff1$contrast_3); hNeff1$contrast_6 <- as.character(hNeff1$contrast_6); hNeff1$contrast_7 <- as.character(hNeff1$contrast_7);
hNeff1$comp <- ifelse(hNeff1$contrast_1 != hNeff1$contrast_5 & hNeff1$contrast_2 == hNeff1$contrast_3 &hNeff1$contrast_6 == hNeff1$contrast_7 &hNeff1$contrast_2 == hNeff1$contrast_6,"Temperature","foobar");
hNeff1$comp <- ifelse(hNeff1$contrast_1 == hNeff1$contrast_5 & hNeff1$contrast_2 != hNeff1$contrast_6 &hNeff1$contrast_3 == hNeff1$contrast_7,"Nuclear",hNeff1$comp);
hNeff1$comp <- ifelse(hNeff1$contrast_1 != hNeff1$contrast_5 & hNeff1$contrast_2 == hNeff1$contrast_6 &hNeff1$contrast_3 == hNeff1$contrast_7,"Temperature",hNeff1$comp);
hNeff1$comp <- ifelse(hNeff1$contrast_1 == hNeff1$contrast_5 & hNeff1$contrast_2 == hNeff1$contrast_6 & hNeff1$contrast_3 != hNeff1$contrast_7, "Mito",hNeff1$comp); hNeff1$contrast_4 <- NULL;
hNeff1$Sex <- "Female"; hNeff1$Trait <- "Fertility"

# Now longevity

lM <- ti %>% filter(Sex == "Male", deathStage == "adult")
a<- lmer(data = lM, 
         longevity ~ Temperature * Cyto * Nuclear  + (1 | Family:Clutch),
         na.action = "na.fail" ,REML = FALSE)
dredge(a,rank = "AICc",evaluate = TRUE,fixed = ~ Temperature + Cyto + Nuclear)
a1 <- lmer(data = lM, 
         longevity ~ Temperature * Nuclear + Cyto  + (1 | Family:Clutch),
         na.action = "na.fail" ,REML = FALSE)
lon.M.lmer <- a1
lMeff1 <- plot(eff_size(emmeans(lon.M.lmer, 1 ~ Temperature + Nuclear + Cyto), sigma = sqrt(mean(sigma(lon.M.lmer)^2)), edf = df.residual(lon.M.lmer)),plotit = FALSE); lMeff1 <- lMeff1[-length(lMeff1)];

# now split the Contrast columns and name the "contrast types" in lMeff1 
lMeff1<- cSplit(indt =lMeff1,splitCols = "contrast" ,sep = " "); lMeff1$contrast_2 <- as.character(lMeff1$contrast_2); lMeff1$contrast_3 <- as.character(lMeff1$contrast_3); lMeff1$contrast_6 <- as.character(lMeff1$contrast_6); lMeff1$contrast_7 <- as.character(lMeff1$contrast_7);
lMeff1$comp <- ifelse(lMeff1$contrast_1 != lMeff1$contrast_5 &lMeff1$contrast_2 == lMeff1$contrast_3 &lMeff1$contrast_6 == lMeff1$contrast_7 &lMeff1$contrast_2 == lMeff1$contrast_6,"Temperature","foobar");
lMeff1$comp <- ifelse(lMeff1$contrast_1 == lMeff1$contrast_5 &lMeff1$contrast_2 != lMeff1$contrast_6 &lMeff1$contrast_3 == lMeff1$contrast_7,"Nuclear",lMeff1$comp);
lMeff1$comp <- ifelse(lMeff1$contrast_1 != lMeff1$contrast_5 &lMeff1$contrast_2 == lMeff1$contrast_6 &lMeff1$contrast_3 == lMeff1$contrast_7,"Temperature",lMeff1$comp);
lMeff1$comp <- ifelse(lMeff1$contrast_1 == lMeff1$contrast_5 & lMeff1$contrast_2 == lMeff1$contrast_6 & lMeff1$contrast_3 != lMeff1$contrast_7, "Mito",lMeff1$comp); lMeff1$contrast_4 <- NULL;
lMeff1$Sex <- "Male"; lMeff1$Trait <- "Longevity"


lF <- ti %>% filter(Sex == "Female", deathStage == "adult")
a<- lmer(data = lF, 
         longevity ~ Temperature * Cyto * Nuclear  + (1 | Family:Clutch),
         na.action = "na.fail" ,REML = FALSE)
dredge(a,rank = "AICc",evaluate = TRUE,fixed = ~ Temperature + Cyto + Nuclear)
a1 <- lmer(data = lF, 
           longevity ~ Temperature * Nuclear + Cyto * Nuclear  + (1 | Family:Clutch),
           na.action = "na.fail" ,REML = FALSE)
lon.F.lmer <- a1
lFeff1 <- plot(eff_size(emmeans(lon.F.lmer, 1 ~ Temperature + Nuclear + Cyto), sigma = sqrt(mean(sigma(lon.F.lmer)^2)), edf = df.residual(lon.F.lmer)),plotit = FALSE); lFeff1 <- lFeff1[-length(lFeff1)];

# now split the Contrast columns and name the "contrast types" in lFeff1 
lFeff1<- cSplit(indt =lFeff1,splitCols = "contrast" ,sep = " "); lFeff1$contrast_2 <- as.character(lFeff1$contrast_2); lFeff1$contrast_3 <- as.character(lFeff1$contrast_3); lFeff1$contrast_6 <- as.character(lFeff1$contrast_6); lFeff1$contrast_7 <- as.character(lFeff1$contrast_7);
lFeff1$comp <- ifelse(lFeff1$contrast_1 != lFeff1$contrast_5 &lFeff1$contrast_2 == lFeff1$contrast_3 &lFeff1$contrast_6 == lFeff1$contrast_7 &lFeff1$contrast_2 == lFeff1$contrast_6,"Temperature","foobar");
lFeff1$comp <- ifelse(lFeff1$contrast_1 == lFeff1$contrast_5 &lFeff1$contrast_2 != lFeff1$contrast_6 &lFeff1$contrast_3 == lFeff1$contrast_7,"Nuclear",lFeff1$comp);
lFeff1$comp <- ifelse(lFeff1$contrast_1 != lFeff1$contrast_5 &lFeff1$contrast_2 == lFeff1$contrast_6 &lFeff1$contrast_3 == lFeff1$contrast_7,"Temperature",lFeff1$comp);
lFeff1$comp <- ifelse(lFeff1$contrast_1 == lFeff1$contrast_5 & lFeff1$contrast_2 == lFeff1$contrast_6 & lFeff1$contrast_3 != lFeff1$contrast_7, "Mito",lFeff1$comp); lFeff1$contrast_4 <- NULL;
lFeff1$Sex <- "Female"; lFeff1$Trait <- "Longevity"
eff1 <- rbind(rMeff1, hNeff1, lMeff1, lFeff1); eff1 <- eff1 %>% filter(comp != "foobar"); eff1$comp <- factor(eff1$comp, levels = c("Mito","Nuclear","Temperature"),labels = c("mtDNA","nDNA","Temperature"))

Anova(fert.M.lmer, type =3)
Anova(fert.F.lmer, type =3)
Anova(lon.M.lmer, type =3)
Anova(lon.F.lmer, type =3)
eff1



