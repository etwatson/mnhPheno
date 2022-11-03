library(dplyr)
library(flexsurv)
library(splitstackshape)
library(survival)

fi <- read.table("fertility_data.txt",header = TRUE)
ti <- read.table("longevity_data.txt",header = TRUE)
dat <- ti %>% filter(deathStage == "adult", Sex %in% c("Male","Female"))

# Perform Gompertz regression on cross x nuclear x temperature groupings
fitDat <- NULL
out <- NULL
for (i in unique(dat$Cross)) {
  foi <- dat %>% filter(Cross == i) 
  for (h in unique(dat$Temperature)) {
    toi <- foi %>% filter(Temperature == h) #Temp in Criss
    for (j in unique(dat$Sex)) {
      soi <- toi %>% filter(Sex == j) %>% arrange(longevity) # Sex in Temp in Cross
      if (nrow(soi) == 0) {
        # print("breaking")
        next 
      }
      fitoi <- flexsurvreg(Surv(longevity) ~ 1, data = soi, dist = "gompertz")
      tmp <- t(data.frame(fitoi$coefficients))
      rownames(tmp) <- NULL
      tmp2 <- data.frame(i, h, j, fitoi$N, tmp)
      colnames(tmp2)[1:4] <- c("Cross", "Temperature", "Sex", "N")
      out <- rbind(out, tmp2)
      f <- fitoi$data
      mod <- survfit(f$m$`Surv(longevity)`~1,data = f$m$`Surv(longevity)`)
      fitDat <- rbind(fitDat,data.frame(time=mod$time, surv=mod$surv,upper=mod$upper,lower=mod$lower,cross=i, Temperature=h,Sex=j))
    }
  }
}
# Data wrangling
fitDat <- cSplit(fitDat,"cross", sep = "x");names(fitDat)[7:8]<- c("Cyto","Nuclear");
sum_out.nts <- full_join(out, dat %>% select(Cross), by = "Cross") %>% distinct() %>% na.omit();colnames(sum_out.nts)[5:6] <- c("G", "ln_A");
sum_out.nts <- cSplit(sum_out.nts,"Cross",sep = "x");names(sum_out.nts)[6:7] <- c("Cyto","Nuclear");sum_out <- sum_out.nts;
sum_out.nts <- melt(sum_out, id.vars = c("Cyto","Nuclear", "Temperature","Sex","N"));sum_out.nts$variable <- factor(sum_out.nts$variable, levels = c("ln_A","G"), labels = c("Frailty","Rate of senescence"));

# Perform Kruskal-Wallis rank sum test on various groupings for frailty
kw.ln_A <-  list(kruskal.test(sum_out$ln_A,sum_out$Sex),kruskal.test(sum_out$ln_A,sum_out$Temperature),kruskal.test(sum_out$ln_A,sum_out$Nuclear),kruskal.test(sum_out$ln_A,sum_out$Cyto),kruskal.test(sum_out$ln_A,interaction(sum_out$Temperature,sum_out$Sex)),kruskal.test(sum_out$ln_A,interaction(sum_out$Temperature,sum_out$Nuclear)),kruskal.test(sum_out$ln_A,interaction(sum_out$Temperature,sum_out$Cyto)),kruskal.test(sum_out$ln_A,interaction(sum_out$Sex,sum_out$Cyto)),kruskal.test(sum_out$ln_A,interaction(sum_out$Nuclear,sum_out$Cyto)),kruskal.test(sum_out$ln_A,interaction(sum_out$Sex,sum_out$Nuclear)));
tmp <- sum_out %>% filter(Temperature == "15ºC");
kw.ln_A <-  c(kw.ln_A,list(kruskal.test(tmp$ln_A,tmp$Sex),kruskal.test(tmp$ln_A,tmp$Nuclear),kruskal.test(tmp$ln_A,tmp$Cyto),kruskal.test(tmp$ln_A,interaction(tmp$Nuclear,tmp$Sex)),kruskal.test(tmp$ln_A,interaction(tmp$Cyto,tmp$Nuclear))));
tmp <- sum_out %>% filter(Temperature == "25ºC");
kw.ln_A <-  c(kw.ln_A,list(kruskal.test(tmp$ln_A,tmp$Sex),kruskal.test(tmp$ln_A,tmp$Nuclear),kruskal.test(tmp$ln_A,tmp$Cyto),kruskal.test(tmp$ln_A,interaction(tmp$Nuclear,tmp$Sex)),kruskal.test(tmp$ln_A,interaction(tmp$Cyto,tmp$Nuclear))));
tmp <- sum_out %>% filter(Temperature == "15ºC",Nuclear == "BR");
kw.ln_A <-  c(kw.ln_A,list(kruskal.test(tmp$ln_A,tmp$Sex),kruskal.test(tmp$ln_A,tmp$Cyto)));
tmp <- sum_out %>% filter(Temperature == "25ºC",Nuclear == "BR");
kw.ln_A <-  c(kw.ln_A,list(kruskal.test(tmp$ln_A,tmp$Sex),kruskal.test(tmp$ln_A,tmp$Cyto)));
tmp <- sum_out %>% filter(Temperature == "15ºC",Nuclear == "SD");
kw.ln_A <-  c(kw.ln_A,list(kruskal.test(tmp$ln_A,tmp$Sex),kruskal.test(tmp$ln_A,tmp$Cyto)));
tmp <- sum_out %>% filter(Temperature == "25ºC",Nuclear == "SD");
kw.ln_A <-  c(kw.ln_A,list(kruskal.test(tmp$ln_A,tmp$Sex),kruskal.test(tmp$ln_A,tmp$Cyto)));
names(kw.ln_A) <- c("ln_A.Sex","ln_A.Temperature","ln_A.Nuclear","ln_A.Cyto","ln_A.Temperature.Sex","ln_A.Temperature.Nuclear","ln_A.Temperature.Cyto","ln_A.Sex.Cyto","ln_A.Nuclear.Cyto","ln_A.Sex.Nuclear","ln_A.15.Sex","ln_A.15.Nuclear","ln_A.15.Cyto","ln_A.15.Sex.Nuclear","ln_A.15.Nuclear.Cyto","ln_A.25.Sex","ln_A.25.Nuclear","ln_A.25.Cyto","ln_A.25.Sex.Nuclear","ln_A.25.Nuclear.Cyto","ln_A.BR.15.Sex","ln_A.BR.15.Cyto","ln_A.BR.25.Sex","ln_A.BR.25.Cyto","ln_A.SD.15.Sex","ln_A.SD.15.Cyto","ln_A.SD.25.Sex","ln_A.SD.25.Cyto");

# Perform Kruskal-Wallis rank sum test on various groupings for rate of senescence
kw.G <-  list(kruskal.test(sum_out$G,sum_out$Sex),kruskal.test(sum_out$G,sum_out$Temperature),kruskal.test(sum_out$G,sum_out$Nuclear),kruskal.test(sum_out$G,sum_out$Cyto),kruskal.test(sum_out$G,interaction(sum_out$Temperature,sum_out$Sex)),kruskal.test(sum_out$G,interaction(sum_out$Temperature,sum_out$Nuclear)),kruskal.test(sum_out$G,interaction(sum_out$Temperature,sum_out$Cyto)),kruskal.test(sum_out$G,interaction(sum_out$Sex,sum_out$Cyto)),kruskal.test(sum_out$G,interaction(sum_out$Nuclear,sum_out$Cyto)),kruskal.test(sum_out$G,interaction(sum_out$Sex,sum_out$Nuclear)));
tmp <- sum_out %>% filter(Temperature == "15ºC");
kw.G <-  c(kw.G,list(kruskal.test(tmp$G,tmp$Sex,),kruskal.test(tmp$G,tmp$Nuclear),kruskal.test(tmp$G,tmp$Cyto),kruskal.test(tmp$G,interaction(tmp$Nuclear,tmp$Sex)),kruskal.test(tmp$G,interaction(tmp$Cyto,tmp$Nuclear))));
tmp <- sum_out %>% filter(Temperature == "25ºC");
kw.G <-  c(kw.G,list(kruskal.test(tmp$G,tmp$Sex),kruskal.test(tmp$G,tmp$Nuclear),kruskal.test(tmp$G,tmp$Cyto),kruskal.test(tmp$G,interaction(tmp$Nuclear,tmp$Sex)),kruskal.test(tmp$G,interaction(tmp$Cyto,tmp$Nuclear))));
tmp <- sum_out %>% filter(Temperature == "15ºC",Nuclear == "BR");
kw.G <-  c(kw.G,list(kruskal.test(tmp$G,tmp$Sex),kruskal.test(tmp$G,tmp$Cyto)));
tmp <- sum_out %>% filter(Temperature == "25ºC",Nuclear == "BR");
kw.G <-  c(kw.G,list(kruskal.test(tmp$G,tmp$Sex),kruskal.test(tmp$G,tmp$Cyto)));
tmp <- sum_out %>% filter(Temperature == "15ºC",Nuclear == "SD");
kw.G <-  c(kw.G,list(kruskal.test(tmp$G,tmp$Sex),kruskal.test(tmp$G,tmp$Cyto)));
tmp <- sum_out %>% filter(Temperature == "25ºC",Nuclear == "SD");
kw.G <-  c(kw.G,list(kruskal.test(tmp$G,tmp$Sex),kruskal.test(tmp$G,tmp$Cyto)));
names(kw.G) <- c("G.Sex","G.Temperature","G.Nuclear","G.Cyto","G.Temperature.Sex","G.Temperature.Nuclear","G.Temperature.Cyto","G.Sex.Cyto","G.Nuclear.Cyto","G.Sex.Nuclear","G.15.Sex","G.15.Nuclear","G.15.Cyto","G.15.Sex.Nuclear","G.15.Nuclear.Cyto","G.25.Sex","G.25.Nuclear","G.25.Cyto","G.25.Sex.Nuclear","G.25.Nuclear.Cyto","G.BR.15.Sex","G.BR.15.Cyto","G.BR.25.Sex","G.BR.25.Cyto","G.SD.15.Sex","G.SD.15.Cyto","G.SD.25.Sex","G.SD.25.Cyto");

# Create flat table. Used for Figure 3B and Table S6

kw <- NULL;
if(length(kw.ln_A)==length(kw.G)){n <- length(kw.ln_A)}
for(i in seq(1:n)){
kw  <- rbind(kw,
    data.frame(chi.squared=kw.ln_A[i][[1]]$statistic,df=kw.ln_A[i][[1]]$parameter,p.value=kw.ln_A[i][[1]]$p.value,row.names = names(kw.ln_A[i])),
    data.frame(chi.squared=kw.G[i][[1]]$statistic,df=kw.G[i][[1]]$parameter,p.value=kw.G[i][[1]]$p.value,row.names = names(kw.G[i])))}
kw$ID <- row.names(kw);kw <- cSplit(kw,"ID", sep = ".");names(kw)[4:7]<- c("parameter","term1","term2","term3")

# Kaplan Meier for Figure 3A

surv_L <- survfit(Surv(longevity) ~ Nuclear+Sex, data = dat %>% filter(Temperature == "15ºC"))
surv_H <- survfit(Surv(longevity) ~ Nuclear+Sex, data = dat %>% filter(Temperature == "25ºC"))

# Kruskal-Wallis tests for median longevity for Figure 3A inset

kruskal.test(dat[dat$Temperature == "15ºC" & dat$Nuclear == "BR",]$longevity,dat[dat$Temperature == "15ºC" & dat$Nuclear == "BR",]$Sex)
kruskal.test(dat[dat$Temperature == "15ºC" & dat$Nuclear == "SD",]$longevity,dat[dat$Temperature == "15ºC" & dat$Nuclear == "SD",]$Sex)
kruskal.test(dat[dat$Temperature == "25ºC" & dat$Nuclear == "BR",]$longevity,dat[dat$Temperature == "25ºC" & dat$Nuclear == "BR",]$Sex)
kruskal.test(dat[dat$Temperature == "25ºC" & dat$Nuclear == "SD",]$longevity,dat[dat$Temperature == "25ºC" & dat$Nuclear == "SD",]$Sex)

