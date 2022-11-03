library(dplyr)
library(boot)
library(stringr)
cor.boot <- function(data, k) cor(data[k,],method = "pearson")[1,2]


fi <- read.table("fertility_data.txt",header = TRUE)
ti <- read.table("longevity_data.txt",header = TRUE)
ti <- ti %>% filter(deathStage == "adult", Sex %in% c("Male","Female"))

# summarize data into tables of mean fertility for each Temperature x Nuclear x Cyto x Sex grouping.
maleFertTab <- fi %>% filter(Sex == "Male") %>% group_by(Temperature, Nuclear, Cyto) %>% summarize(fert = mean(fert, na.rm = T));names(maleFertTab)[4]  <- "Male"; maleFertTab$Temperature <- factor(maleFertTab$Temperature);
femaleFertTab <- fi %>% filter(Sex == "Female") %>% group_by(Temperature, Nuclear, Cyto) %>% summarize(fert = mean(fert, na.rm = T)); names(femaleFertTab)[4] <- "Female"; femaleFertTab$Temperature <- factor(femaleFertTab$Temperature)
fert.tab <- merge(maleFertTab,femaleFertTab); l.fertTab <- split.data.frame(fert.tab,interaction(fert.tab$Temperature,fert.tab$Nuclear))
maleLonTab <- ti %>% filter(deathStage == "adult",Sex == "Male") %>% group_by(Temperature, Nuclear, Cyto) %>% summarize(lon = mean(longevity, na.rm = T));names(maleLonTab)[4] <- "Male";
femaleLonTab <- ti %>% filter(deathStage == "adult",Sex == "Female") %>% group_by(Temperature, Nuclear, Cyto) %>% summarize(lon = mean(longevity, na.rm = T));names(femaleLonTab)[4] <- "Female";
lon.tab <- merge(maleLonTab,femaleLonTab); l.lonTab <- split.data.frame(lon.tab,interaction(lon.tab$Temperature,lon.tab$Nuclear))

mIa <- lon.tab[,-5];names(mIa)[4]<- "longevity";mIb <- fert.tab[,-5];names(mIb)<- c("Temperature","Nuclear","Cyto","fertility");mI <- merge(mIa,mIb);l.mI <- split.data.frame(mI,interaction(mI$Temperature,mI$Nuclear));
fIa <- lon.tab[,-4];names(fIa)[4]<- "longevity";fIb <- fert.tab[,-4];names(fIb)<- c("Temperature","Nuclear","Cyto","fertility");fI <- merge(fIa,fIb);l.fI <- split.data.frame(fI,interaction(fI$Temperature,fI$Nuclear));

# get bootstrapped measure of Pearson's r for each intra- and intersexual correlation.
sexCor <- NULL
for (i in seq(1:length(l.mI))){
  t <- str_split(names(l.mI[i]),"\\.")[[1]][1];n <- paste(str_split(names(l.mI[i]),"\\.")[[1]][2],"nDNA",sep = " ");
  b <- boot(l.mI[i][[1]][4:5],statistic = cor.boot,R=1000); ci <- boot.ci(b, type = "norm");   df <- data.frame(Temperature = t,Nuclear = n,Sex = "Intrasex (Males)", Comp = "FvL", coef = ci$t0,lower.CI = ci$normal[2],upper.CI = ci$normal[3],Name = "M-Fert-v-Lon")
  sexCor <- rbind(sexCor,df);}
for (i in seq(1:length(l.fI))){
  t <- str_split(names(l.fI[i]),"\\.")[[1]][1];n <- paste(str_split(names(l.fI[i]),"\\.")[[1]][2],"nDNA",sep = " ");
  b <- boot(l.fI[i][[1]][4:5],statistic = cor.boot,R=1000); ci <- boot.ci(b, type = "norm");
  df <- data.frame(Temperature = t,Nuclear = n,Sex = "Intrasex (Females)", Comp = "FvL", coef = ci$t0,lower.CI = ci$normal[2],upper.CI = ci$normal[3],Name = "F-Fert-v-Lon")
  sexCor <- rbind(sexCor,df);}
for (i in seq(1:length(l.lonTab))){
  t <- str_split(names(l.lonTab[i]),"\\.")[[1]][1];n <- paste(str_split(names(l.lonTab[i]),"\\.")[[1]][2],"nDNA",sep = " ");
  b <- boot(l.lonTab[i][[1]][4:5],statistic = cor.boot,R=1000); ci <- boot.ci(b, type = "norm");
  df <- data.frame(Temperature = t,Nuclear = n,Sex = "Intersex", Comp = "L", coef = ci$t0,lower.CI = ci$normal[2],upper.CI = ci$normal[3],Name = "Longevity")
  sexCor <- rbind(sexCor,df);}
for (i in seq(1:length(l.fertTab))){
  t <- str_split(names(l.fertTab[i]),"\\.")[[1]][1];n <- paste(str_split(names(l.fertTab[i]),"\\.")[[1]][2],"nDNA",sep = " ");
  b <- boot(l.fertTab[i][[1]][4:5],statistic = cor.boot,R=1000); ci <- boot.ci(b, type = "norm");
  df <- data.frame(Temperature = t,Nuclear = n,Sex = "Intersex", Comp = "F", coef = ci$t0,lower.CI = ci$normal[2],upper.CI = ci$normal[3],Name = "Fertility")
  sexCor <- rbind(sexCor,df);}
sexCor$Sex <- factor(sexCor$Sex,levels = c("Intersex","Intrasex (Females)","Intrasex (Males)"));sexCor$Comp <- factor(sexCor$Comp,levels = c("F","L","FvL"));sexCor$Temperature <- factor(sexCor$Temperature);sexCor$Name <- factor(sexCor$Name, levels = c("F-Fert-v-Lon","M-Fert-v-Lon","Fertility","Longevity"))
