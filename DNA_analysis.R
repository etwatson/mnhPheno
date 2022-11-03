#### mtDNA contamination check
library(reshape2)
library(dplyr)

dep <- read.table("SDv2.2_AllMt.multicov", head = FALSE) # Read depth mapping to SDv2.2 nuclear genome + all mt haplotypes 

names(dep) <- c("chr","start","end","BRxBR","BRxSD","SDxBR","SDxSD","FHLxBR","FHLxSD","LJxBR","LJxSD")
dep <- melt(dep, id.vars = c("chr","start","end")) 
names(dep)[4:5] <- c("sample","depth")
dep$chr <- factor(dep$chr, levels = c("BR4","LJ5","SD202","FHL9"),
                  labels = c("BR","LJ","SD","FHL"))

dep.t <- dep %>% group_by(sample,chr) %>% summarise(depth = mean(depth))
dcast(dep.t, formula = chr ~ sample)
head(dep)
