# Phenotyping of mitonuclear hybrids


**Original code for the article titled “Mitochondrial effects on fertility and longevity in Tigriopus californicus contradict predictions of the mother’s curse hypothesis”**
---
## Background
*Recent empirical work* suggests that, as males do not transmit their mitochondria to subsequent generations, mitochondrial genomes in eukaryotes should accumulate male-harming mutations via a sex-specific selective sieve acting upon female fitness (mother’s curse). Mitochondrial disease variants have been identified in humans (e.g. Leber’s optical neuropathy, infertility) and specific mitochondrial lines causing faster aging or lower fecundity with strong male bias have been identified in a few other species. However, experimental results have been mixed, with some mitonuclear combinations showing no male bias or even the effects of female-biased load. While the logic of mother’s curse may appear indisputable alongside the supporting evidence, many factors may counteract the build-up of male-specific load or compensate for its adverse effects. Inbreeding, kin selection, mitonuclear epistasis, and paternal leakage are a few of the countervailing forces identified. Although maternal transmission of mitochondria is widespread among eukaryotes, it remains unclear how pervasive and important mother’s curse may be across the tree of life. 

*In this paper*, we utilize mitonuclear hybrid constructs in the marine copepod *Tigriopus californicus* to investigate the impact of mitochondrial variation on sex-specific fertility and longevity. *T. californicus* are well known for their extreme mitochondrial genetic divergence among geographically proximate yet isolated populations, making them ideal for such investigations as populations remain interfertile. Additionally, while sex-biased epistatic effects map to X-linked nuclear genes in similar studies using *Drosophila mitonuclear* hybrids, *T. californicus* do not possess sex chromosomes, thus simplifying the interpretation of mitochondrial load. 

*Our study* combined four divergent mtDNA haplotypes and two divergent nuclear lines, with the resulting mitonuclear hybrids being reared in two temperature regimes. We measured sex-specific fertility in 858 age-matched individuals and longevity in 3,072 individuals. Non-parametric modeling was used to measure the impact of genotype, environment, and genotype x environment as well as interactions with sex. We also performed parametric survivorship modeling. With these data, we sought to test the main predictions of the mother’s curse hypothesis.

*We found* that mitochondrial genetic variance was overwhelmingly female-biased, sexually-antagonistic mitochondrial variance was rare, and there was no influence of coevolutionary status on fertility or longevity. We consider the importance of inbreeding as well as mating system and the significance of performing a taxonomically diverse investigation of the generality of mother’s curse. 

## Code

* coevoStat_ind.R - modeling of coevolutionary status
* DNA_analysis.R - measuring potential for DNA contamination 
* FSR.R - modeling of family sex ratio data
* InbVSout.R - measuring the impact of inbreeding vs outbreeding on female hatching number
* ISC.R - bootstrapped Pearsons coefficient for inter- and intrasexual correlation for fertility and longevity
* juvMort.R - measuring juvenile mortality and modeling copepodid duration, age at copepodid transition, and age of sexual maturity from juvenile mortality data
* MtCoefVar.R - bootstrapped sign tests for mitochondrial coefficient of variation
* mtEff.R - best-fit models of fertility and longevity and calculating effect sizes of genetic and environmental variaton. 
* SexEffect.R - model testing of four sex-interaction models, sex contrasts for each model
* SurvReg.R - Gompertz regression analysis, Kruskal Wallis tests for Gompertz parameters, Kaplan Meier plot in Figure 3A

## Citation
> **Watson, Eric T.**; Flanagan, B.; Pascar, J.A.; Edmands, S. *in press*. Proceedings of the Royal Society B: Biological Sciences**
