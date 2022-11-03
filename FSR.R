library(car)
library(lme4)

dc <- read.table("family_sex_ratio_data.txt",header = TRUE)

# Genotype and environment model
mod0 <- lmer(csr ~ Temperature * Cyto * Nuclear + (1| Family) ,data = dc, weights = Family.size,na.action = 'na.fail',REML = FALSE)
dredge(mod0, evaluate = TRUE,fixed = ~ Cyto + Nuclear + Temperature);
mod.SR  <- lmer(csr ~ Temperature * Nuclear + Cyto * Nuclear + (1| Family) ,data = dc, weights = Family.size,na.action = 'na.fail')
mod.correctedSR  <- lmer(corrected.csr ~ Temperature * Nuclear + Cyto * Nuclear + (1| Family) ,data = dc, weights = Family.size,na.action = 'na.fail')
Anova(mod.SR, type =3)
Anova(mod.correctedSR, type =3)

# Genotype, environment, and coevolutionary status model
mod0 <- lmer(csr ~ Temperature * coevo + Nuclear*coevo + Temperature*Nuclear + Cyto + (1| Family) ,data = dc, weights = Family.size,na.action = 'na.fail',REML = FALSE)
dredge(mod0, evaluate = TRUE,fixed = ~ coevo + Temperature)
mod.CSR <- lmer(csr ~ Nuclear*coevo + Temperature*Nuclear +Cyto +(1| Family) ,data = dc, weights = Family.size,na.action = 'na.fail',REML = FALSE)
mod.correctedCSR <- lmer(corrected.csr ~ Nuclear*coevo + Temperature*Nuclear +Cyto +(1| Family) ,data = dc, weights = Family.size,na.action = 'na.fail',REML = FALSE)
Anova(mod.CSR, type =3)
Anova(mod.correctedCSR, type =3)

# Calculate marginal means (95% CI) for each haplotype.
srTab <- emmip(emmeans(mod.SR, 1~ Cyto|Nuclear+Temperature),Cyto ~ Nuclear+Temperature,CIs = TRUE, plotit = FALSE); srTab$Nuclear <- factor(srTab$Nuclear, levels = c("BR","SD"),labels = c("BR nDNA","SD nDNA"));srTab$coevo <- ifelse(as.character(srTab$Cyto)== as.character(gsub(" nDNA","",srTab$Nuclear)),"matched","mismatched"); srTab$coevo <- factor(srTab$coevo); srTab$Cyto <- factor(srTab$Cyto, levels = c("SD","BR","LJ","FHL"));
csrTab <- emmip(emmeans(mod.CSR, 1~ coevo|Nuclear+Temperature),Nuclear ~ coevo+Temperature,CIs = TRUE,plotit = FALSE); csrTab$Nuclear <- factor(csrTab$Nuclear, levels = c("BR","SD"),labels = c("BR nDNA","SD nDNA"));

