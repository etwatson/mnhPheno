library(reshape2)
Cy <- read.table("female_fert_data.txt", head = TRUE) 
Cy.cross <- dcast(Cy, Temperature + mtDNA + nDNA  ~ cross, value.var = 'n.offspring',mean) #get mean hatching number for each experimental grouping
Cy.cross$p.out <- Cy.cross$inbred/Cy.cross$outcross 
Cy.cross$mtDNA <- factor(Cy.cross$mtDNA, levels = c("SD","LJ","BR","FHL"))
Cy.cross$nDNA <- factor(Cy.cross$nDNA, levels = c("BR","SD"),labels = c("BR nDNA", "SD nDNA"))

