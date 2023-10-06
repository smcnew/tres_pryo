# Comparing results of NGS with pyro
#
library(tidyr)
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(car)
library(ggplot2)
library(sjPlot)
library(ggthemr)
library(patchwork) # plot layouts


# Data --------------------------------------------------------------------


# 'Meth_Data' is the main dataset with a row per bird and methylation at each site.
crhr <-read.delim("CRH1/Meth_Data.txt")
# 'Align' is location data for where the CPGs line up in the gene
align_crhr<-read.delim("CRH1/Align_to_CRH1.txt")
# 'assays' has info on where the primer pair assays fall with respect to CpG sites
assays<-read.delim("CRH1/CRH1_Assay_Locs.txt")


head(crhr)

crhr_long <- crhr %>% pivot_longer(cols = !c(Band, b.bright, Age, Sample, Treat1,
                                            Treat2, Treatment, Customer.ID,
                                            base1, stress1, dex1, base2, base3,
                                            stress3, dex3),
                                   names_to = "cpg", values_to = "perc_meth") %>%
  mutate(logit_meth = logit(perc_meth, adjust = 0.01))




# Data analysis -----------------------------------------------------------

head(crhr_long)
lmer(perc_meth ~ Treat1 + Treat2 + b.bright + (1|Band) + (1|cpg), data = crhr_long) %>% summary
