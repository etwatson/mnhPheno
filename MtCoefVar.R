library(dplyr)
library(BSDA)
library(boot)
library(splitstackshape)
samplemean <- function(x, d) {return(mean(x[d],na.rm = TRUE))};
samplemedian <- function(x, d) {return(median(x[d]))};

fi <- read.table("fertility_data.txt",header = TRUE)
ti <- read.table("longevity_data.txt",header = TRUE)
ti <- ti %>% filter(deathStage == "adult", Sex %in% c("Male","Female"))

# Fertility
tmp <- fi
tmp.l <- split.data.frame(x = tmp , interaction(tmp$Temperature,tmp$Nuclear,tmp$Sex))
mtBoot <- NULL
for (i in seq(1:8)){ #get bootstrapped haplotype means
  BR <- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "BR"),]$fert,samplemean,R = 1000)$t) 
  SD <- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "SD"),]$fert,samplemean,R = 1000)$t)
  LJ <- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "LJ"),]$fert,samplemean,R = 1000)$t)
  FHL<- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "FHL"),]$fert,samplemean,R = 1000)$t)
  for (n in seq(1:1000)){
    # Calculate CV from each 1000 samples
    Vec <- c(BR[n],SD[n],LJ[n],FHL[n])
    mtBoot <- rbind(mtBoot, data.frame(ID=names(tmp.l[i]),CV = sd(Vec)/mean(Vec))) 
  }
}

mtBoot <- cSplit(mtBoot, 'ID', sep=".", type.convert=FALSE)
names(mtBoot) <- c("CV","Temperature","Nuclear", "Sex")
mtBoot$CV <- as.numeric(as.character(mtBoot$CV))

mtBoot.l <- split.data.frame(mtBoot, interaction(mtBoot$Temperature,mtBoot$Nuclear))
signDat <- NULL
signT <- NULL
for(i in seq(1:4)){ #perform sign test comparing males and females 
  r1 <-  BSDA::SIGN.test(x=mtBoot.l[[i]][with(mtBoot.l[[i]],Sex == "Female"),]$CV,
                         y=mtBoot.l[[i]][with(mtBoot.l[[i]],Sex == "Male"),]$CV,alternative = "two.sided",conf.int = 0.95)
  signDat[[length(signDat) + 1]] <- r1
  signT <- rbind(signT, data.frame(Phenotype=names(mtBoot.l[i]),
                                   p.value=r1$p.value,estimate=r1$estimate[[1]],
                                   LconfINT=r1$conf.int[1],UconfINT=r1$conf.int[2]))
}
names(signDat) <- names(mtBoot.l)
signT <- cSplit(signT, 'Phenotype', sep=".", type.convert=FALSE)
names(signT)[c(5:6)] <- c("Temperature","Nuclear")
signT$log.p <- -log10(signT$p.value)
Fert.signT <- signT

# Longevity 
tmp <- ti
tmp.l <- split.data.frame(x = tmp , interaction(tmp$Temperature,tmp$Nuclear,tmp$Sex))
mtBoot <- NULL
for (i in seq(1:8)){
  BR <- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "BR"),]$longevity,samplemedian,R = 1000)$t)
  SD <- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "SD"),]$longevity,samplemedian,R = 1000)$t)
  LJ <- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "LJ"),]$longevity,samplemedian,R = 1000)$t)
  FHL<- as.numeric(boot(tmp.l[[i]][with(tmp.l[[i]],Cyto == "FHL"),]$longevity,samplemedian,R = 1000)$t)
  for (n in seq(1:1000)){
    # Calculate CV from each 1000 samples
    Vec <- c(BR[n],SD[n],LJ[n],FHL[n])
    mtBoot <- rbind(mtBoot, data.frame(ID=names(tmp.l[i]),CV = sd(Vec)/mean(Vec)))
  }
}

mtBoot <- cSplit(mtBoot, 'ID', sep=".", type.convert=FALSE)
names(mtBoot) <- c("CV","Temperature","Nuclear", "Sex")
mtBoot$CV <- as.numeric(as.character(mtBoot$CV))

mtBoot.l <- split.data.frame(mtBoot, interaction(mtBoot$Temperature,mtBoot$Nuclear))
signDat <- NULL
signT <- NULL
for(i in seq(1:4)){
  r1 <-  BSDA::SIGN.test(x=mtBoot.l[[i]][with(mtBoot.l[[i]],Sex == "Female"),]$CV,
                         y=mtBoot.l[[i]][with(mtBoot.l[[i]],Sex == "Male"),]$CV,alternative = "two.sided",conf.int = 0.95)
  signDat[[length(signDat) + 1]] <- r1
  signT <- rbind(signT, data.frame(Phenotype=names(mtBoot.l[i]),
                                   p.value=r1$p.value,estimate=r1$estimate[[1]],
                                   LconfINT=r1$conf.int[1],UconfINT=r1$conf.int[2]))
}
names(signDat) <- names(mtBoot.l)
signT <- cSplit(signT, 'Phenotype', sep=".", type.convert=FALSE)
names(signT)[c(5:6)] <- c("Temperature","Nuclear")
signT$log.p <- -log10(signT$p.value)
Lon.signT <- signT

Fert.signT$Trait <- "Fertility"
Lon.signT$Trait <- "Median longevity"

signT <- rbind(Fert.signT, Lon.signT)
