
# Analysis of pyrosequencing data from McNew et al. 2023

library(tidyr) # data manipulation
library(dplyr) # data manipulation
library(lme4) # for LMMs
library(lmerTest) # for p values for LMMs
library(car) # stats
library(rptR) # assessing variance explained by fixed and random effects

library(ggplot2) # plotting
library(ggthemr) # plot themes
library(patchwork) # plot layouts
library(sjPlot) # extracting tables


# Plot params -------------------------------------------------------------

ggthemr(palette = "flat", layout = "minimal", text_size = 15)
#plot(swatch())

#
# Data --------------------------------------------------------------------

# Methylation data from pyrosequencing
# Just subset to year of study (2017) and minor data processing
ds <- read.delim("raw_data/Input_Samples.txt")
ds <- subset(ds, ds$Year == 2017)
ds[ds == 0] <- NA
ds <- mutate(ds, Capture = as.factor(Capture))
(ds$Band) %>% unique %>% length # 70 birds

# Coordinate data from samples: what parts of genes were sequenced?
input_samples <- read.delim("raw_data/input_CpGs.txt") %>% rename(cpg = CpG) %>%
  select(-Note)# new = old
head(input_samples)
table(input_samples$Gene) # loci per gene

# Pivot table longer, add information for gene

ds_long <- pivot_longer(ds, cols = !c(SampleID, BibB1, Band, Year, Capture, Treatment),
                        names_to = "cpg", values_to = "perc_meth") %>%
          mutate(gene = gsub('[[:digit:]]+', '', cpg)) %>% # remove number from smps to isolate gene name
          mutate(logit_meth = logit(perc_meth, adjust = 0.01)) %>%
          left_join(., input_samples)

# ds_long$cpg %in% input_samples$cpg %>% table # double check all snp info present
head(ds_long)

# Covariates for fledging success and cort
covars <- read.delim("raw_data/Input_Cort_Fled.txt") %>% rename(Band = fband) %>%
  mutate(Capture = as.factor(Capture)) %>% filter(year == "2017")

ds_long <- left_join(ds_long, covars) %>% as.data.frame() # add covars to meth data


# Add treatment and brightness info to covar set
covars <- left_join(covars, select(ds_long, Band, Treatment, BibB1, Capture)) %>%
  distinct %>% filter(!is.na(Treatment))
aggregate(BibB1 ~ Treatment, covars, mean) # check for variation in starting brightness

# calculate change in methylation between step 1 and 2

ds_wide <- ds_long %>% select (BibB1, Band, Year, Capture, Treatment, cpg,
                               perc_meth, Gene) %>%
  pivot_wider(names_from = Capture, values_from = perc_meth,
              names_prefix = "capture_") %>%
  mutate(deltameth = capture_3 - capture_1) %>%
  mutate(logitdelta = logit(capture_3, adjust = 0.01) - logit(capture_1, adjust = 0.01))


# inspection of logit vs. untransformed data  -----------------------------
# untransformed data
fit2 <- ds_long %>% filter(gene == "CRHR") %>%
  lmer(perc_meth ~ Treatment * Capture * BibB1 + (1|Band) + (1|cpg), data = .)
plot(fit2)
qqnorm(residuals(fit2)) # fit is not great
qqline(residuals(fit2))

par(mfrow=c(2,4))

hist(ds_long$perc_meth[ds_long$gene == "CRHn"])
hist(ds_long$perc_meth[ds_long$gene == "CRHR"])
hist(ds_long$perc_meth[ds_long$gene == "FKBP"])
hist(ds_long$perc_meth[ds_long$gene == "GRn"])


# Logit transformation is indicated for percentage data
hist(logit(ds_long$perc_meth[ds_long$gene == "CRHn"]))
hist(logit(ds_long$perc_meth[ds_long$gene == "CRHR"]))
hist(logit(ds_long$perc_meth[ds_long$gene == "FKBP"]))
hist(logit(ds_long$perc_meth[ds_long$gene == "GRn"]))


# Summary stats -----------------------------------------------------------

# Number of loci sequenced per bird
ds_long %>% filter(Capture ==1 ) %>% group_by(Treatment, Band) %>% tally
# Number of individual samples
paste (ds_long$Band, ds_long$Capture, sep = "_") %>%
  unique  %>% length #121 samples
# Number of birds
ds_long %>% filter(Capture == 1 ) %>%
  select(Treatment, Band) %>%
  distinct %>%
  group_by(Treatment) %>% tally

# gene and CpG information
ds_long %>% select(cpg) %>% unique %>% dim #56 sites
ds_long %>% select(Gene) %>% unique %>% dim #4 genes
ds_long %>% distinct(cpg, .keep_all=T) %>% group_by(Gene) %>% tally # how many sites per gene

# Completeness
table(ds$Band, ds$Capture) %>% apply(., 2, sum) # 70 capture 1 samples, 51 capture 2, total = 121
methylationsummary <- apply(ds[,7:62], 2, summary) %>% t %>% as.data.frame
colnames(methylationsummary) <- c("min", "first", "median", "mean", "third", "max", "nas")
methylationsummary <- methylationsummary %>% mutate(n = 121 - nas )
methylationsummary %>% pull(n) %>% fivenum
methylationsummary %>% pull(mean) %>% fivenum
# write out summary for supplement
methylationsummary %>% select(-nas) %>% write.csv("output_plots/methylationsummary.csv")

# Model percent methylation at each gene ----------------------------------


# Full models

# ds_long %>% filter(gene == "FKBP") %>%
#   lmer(perc_meth ~ Treatment + Capture + BibB1 + (1|Band) + (1|cpg), data = .) %>% summary
# ds_long %>% filter(gene == "GRn") %>%
#   lmer(perc_meth ~ Treatment + Capture + BibB1 + (1|Band) + (1|cpg), data = .) %>% summary
#   ds_long %>% filter(gene == "CRHn") %>%
# lmer(perc_meth ~ Treatment * Capture * BibB1  + (1|Band) + (1|cpg), data = .) %>% summary


# Methylation changes between
# No effect of treatment or capture
CRHn_mod <- ds_long %>% filter(gene == "CRHn") %>%
            lmer(logit_meth ~ Treatment +  Capture + BibB1 + (1|Band) + (1|cpg), data = .)

rep5 <- rpt(logit_meth ~ Treatment +  Capture + BibB1 + (1|Band) + (1|cpg),
             grname = c("Band", "cpg", "Fixed"),
            data = ds_long[ds_long$gene == "CRHn",],
            datatype = "Gaussian",
            nboot = 1000,
            npermut = 0, ratio = TRUE)
rep5

# brightness * capture * Treatment all significant
CRHR_mod <- ds_long %>% filter(gene == "CRHR") %>%
  lmer(logit_meth ~ Treatment * Capture * BibB1 + (1|Band) + (1|cpg), data = .)

rep6 <- rpt(logit_meth ~ Treatment * Capture * BibB1 + (1|Band) + (1|cpg),
            grname = c("Band", "cpg", "Fixed"),
            data = ds_long[ds_long$gene == "CRHR",],
            datatype = "Gaussian",
            nboot = 1000,
            npermut = 0, ratio = TRUE)

# Methylation significantly decreases in capture 3, no effect of treatment
FKBP_mod <- ds_long %>% filter(gene == "FKBP") %>%
  lmer(logit_meth ~ Treatment + Capture  + BibB1 + (1|Band) + (1|cpg), data = .)
rep7 <- rpt(logit_meth ~ Treatment + Capture + BibB1 + (1|Band) + (1|cpg),
            grname = c("Band", "cpg", "Fixed"),
            data = ds_long[ds_long$gene == "FKBP",],
            datatype = "Gaussian",
            nboot = 1000,
            npermut = 0, ratio = TRUE)
rep7

# Methylation significantly increases in capture 3, no effect of treatment
GRn_mod <- ds_long %>% filter(gene == "GRn") %>%
  lmer(logit_meth ~ Treatment + Capture  + BibB1 + (1|Band) + (1|cpg), data = .)
rep8 <- rpt(logit_meth ~ Treatment + Capture + BibB1 + (1|Band) + (1|cpg),
            grname = c("Band", "cpg", "Fixed"),
            data = ds_long[ds_long$gene == "GRn",],
            datatype = "Gaussian",
            nboot = 1000,
            npermut = 0, ratio = TRUE)

rep8
tab_model(CRHn_mod, CRHR_mod,FKBP_mod, GRn_mod, auto.label = FALSE,
          dv.labels = c("CRHn", "CRHR", "FKBP", "GRn"),
          pred.labels = c("Intercept", "Treatment[Dulled]",
                        "Capture number",
                        "Initial brightness",
                        "Treatment:Capture",
                        "Treatment:Brightness",
                        "Capture:Initial brightness",
                        "Treatment:Capture:Initial brightness"))
# Relationship with cort ---------------------------------------------------------------------
# Is methylation related to baseline corticosterone levels?
modb1 <- ds_long %>% filter(Gene == "CRH") %>% lmer(logit_meth ~ base + (1|cpg) + (1|Band), data = . )
modb2 <- ds_long %>% filter(Gene == "CRHR1") %>% lmer(logit_meth ~ base + (1|cpg) + (1|Band), data = . )
modb3 <- ds_long %>% filter(Gene == "FKBP5") %>% lmer(logit_meth ~ base + (1|cpg) + (1|Band), data = . )
modb4 <- ds_long %>% filter(Gene == "GR") %>% lmer(logit_meth ~ base + (1|cpg) + (1|Band), data = . )

tab_model(modb1, modb2, modb3, modb4,
          dv.labels = c("CRH", "CRHR1", "FKBP5", "GR"),
          pred.labels = c("Intercept", "Baseline corticosterone"))

# Relationship with stress induced corticosterone
mods1 <- ds_long %>% filter(Gene == "CRH") %>% lmer(logit_meth ~ stress + (1|cpg) + (1|Band), data = . )
mods2 <- ds_long %>% filter(Gene == "CRHR1") %>% lmer(logit_meth ~ stress + (1|cpg) + (1|Band), data = . )
mods3 <- ds_long %>% filter(Gene == "FKBP5") %>% lmer(logit_meth ~ stress + (1|cpg) + (1|Band), data = . )
mods4 <- ds_long %>% filter(Gene == "GR") %>% lmer(logit_meth ~ stress + (1|cpg) + (1|Band), data = . )

tab_model(mods1, mods2, mods3, mods4,
          dv.labels = c("CRH", "CRHR1", "FKBP5", "GR"),
          pred.labels = c("Intercept", "Stress-induced corticosterone"))
# Relationship with dex
modd1 <- ds_long %>% filter(Gene == "CRH") %>% lmer(logit_meth ~ dex + (1|cpg) + (1|Band), data = . )
modd2 <- ds_long %>% filter(Gene == "CRHR1") %>% lmer(logit_meth ~ dex + (1|cpg) + (1|Band), data = . )
modd3 <- ds_long %>% filter(Gene == "FKBP5") %>% lmer(logit_meth ~ dex + (1|cpg) + (1|Band), data = . )
modd4 <- ds_long %>% filter(Gene == "GR") %>% lmer(logit_meth ~ dex + (1|cpg) + (1|Band), data = . )
tab_model(modd1, modd2, modd3, modd4,
          dv.labels = c("CRH", "CRHR1", "FKBP5", "GR"),
          pred.labels = c("Intercept", "Post-dex corticosterone"))

# interaction plot CRHR -------------------------------------------------------------------




# Plot of interaction between brightness, treatment, and change in methylation
# for CRHR gene
crhr_interaction_plot <- ds_wide %>% filter(Gene == "CRHR1") %>%
  ggplot(aes(x = BibB1, y = deltameth, color = Treatment)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", aes(fill = Treatment)) +
  scale_color_manual(values = c("#e74c3c", "#f1c40f")) +
  scale_fill_manual(values = c("#e74c3c", "#f1c40f"  ))+
  labs(x = "Original plumage brightness (% reflectance)", y = "Change in methylation (%)")

pdf("output_plots/CRHR_plumagexmethxprepost.pdf", height = 7, width = 9)
crhr_interaction_plot
dev.off()


plot_allgenes <-
  ds_long %>%
  ggplot(aes(y= perc_meth, x = Gene, fill = Capture)) +
  labs(x = "Gene", y = "Percent methylation", fill = "Capture") +
  geom_violin(trim = T) +
  geom_boxplot(width = 0.08, position = position_dodge(width = 0.9), show.legend = F) +
  scale_fill_discrete(labels=c('Pre-treatment', 'Post-treatment')) +
  annotate(geom="text", x=2, y=50, label="*", size = 10) +
  annotate(geom="text", x=3, y=104, label="*", size = 10) +
  annotate(geom="text", x=4, y=85, label="*", size = 10) +
  geom_segment(aes(x = 1.8, xend = 2.2, y = 49, yend = 49)) +
  geom_segment(aes(x = 2.8, xend = 3.2, y = 103, yend = 103)) +
  geom_segment(aes(x = 3.8, xend = 4.2, y = 84, yend = 84))
pdf("output_plots/temporal_meth_change.pdf", height = 7, width = 9)
plot_allgenes
dev.off()




# Individual Genes smooth plots --------------------------------------------------------------------
head(ds_long)
p1 <- ds_long %>% filter(gene == "CRHR", Capture == "1") %>%
  ggplot(aes(y = perc_meth, x = FromTSS, color = Treatment )) +
  geom_smooth(aes(fill = Treatment)) +
  geom_point(alpha = 0.4) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Treatment",
       title ="Pre-manipulation") +
  scale_color_manual(values = c("#e74c3c", "#f1c40f"))  +
  scale_fill_manual(values = c("#e74c3c", "#f1c40f")) +
  theme(plot.title = element_text(size=15))


p2 <- ds_long %>% filter(gene == "CRHR", Capture == "3") %>%
  ggplot(aes(y = perc_meth, x = FromTSS, color = Treatment )) +
  geom_smooth(aes(fill = Treatment)) +
  geom_point(alpha = 0.4) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Treatment",
       title ="Post-manipulation") +
  scale_color_manual(values = c("#e74c3c", "#f1c40f"))  +
  scale_fill_manual(values = c("#e74c3c", "#f1c40f"))+
  theme(plot.title = element_text(size=15))

pdf("output_plots/CRHR_pre_post.pdf", width = 15)
p1 + p2
dev.off()


p3 <- ds_long %>% filter(gene == "CRHn") %>%
  ggplot(aes(y = perc_meth, x = FromTSS, color = Capture )) +
  geom_smooth(aes(fill = Capture)) +
  geom_point(alpha = 0.4) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title = "CRH") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment"))+
  theme(legend.position = "none", plot.title = element_text(size=15))

p4 <- ds_long %>% filter(gene == "FKBP") %>%
  ggplot(aes(y = perc_meth, x = FromTSS, color = Capture )) +
  geom_smooth(aes(fill = Capture)) +
  geom_point(alpha = 0.4) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title = "FKBP5") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment"))+
  theme(legend.position = "none", plot.title = element_text(size=15))

p5 <- ds_long %>% filter(gene == "GRn") %>%
  ggplot(aes(y = perc_meth, x = FromTSS, color = Capture )) +
  geom_smooth(aes(fill = Capture)) +
  geom_point(alpha = 0.4) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title ="GR") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment"))+
  theme(plot.title = element_text(size=15))


pdf("output_plots/pre_post_other_genes.pdf", width = 20)
p3 + p4 + p5
dev.off()


# Boxplots methylation per cpg --------------------------------------------
ds_long$gene %>% table

p6 <- ds_long %>% filter(gene == "CRHn") %>%
  ggplot(aes(y = perc_meth, x = as.factor(FromTSS), fill = Capture)) +
  geom_point(aes(color = Capture), position = position_jitterdodge(), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title = "CRH") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment")) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.position = "none", plot.title = element_text(size=15))
p7 <- ds_long %>% filter(gene == "CRHR") %>%
  ggplot(aes(y = perc_meth, x = as.factor(FromTSS), fill = Capture)) +
  geom_point(aes(color = Capture), position = position_jitterdodge(), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title = "CRHR1") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none", plot.title = element_text(size=15))



p8 <- ds_long %>% filter(gene == "FKBP") %>%
  ggplot(aes(y = perc_meth, x = as.factor(FromTSS), fill = Capture)) +
  geom_point(aes(color = Capture), position = position_jitterdodge(), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA) +
   labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title = "FKBP5") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none", plot.title = element_text(size=15))

p9 <- ds_long %>% filter(gene == "GRn") %>%
  ggplot(aes(y = perc_meth, x = as.factor(FromTSS), fill = Capture)) +
  geom_point(aes(color = Capture), position = position_jitterdodge(), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Capture",
       title = "GR") +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(size=15))

pdf("output_plots/boxplotspergene.pdf", width = 20, height = 12)
(p6 + p7) /
(p8 + p9)  + plot_annotation(tag_levels = 'A')
dev.off()


#boxplots pre post for CRHR1
p10 <- ds_long %>% filter(gene == "CRHR", Capture == "1") %>%
  ggplot(aes(y = perc_meth, x = as.factor(FromTSS), fill = Treatment )) +
  geom_point(aes(color = Treatment), position = position_jitterdodge(), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Treatment",
       title ="Pre-manipulation") +
  scale_color_manual(values = c("#e74c3c", "#f1c40f"))  +
  scale_fill_manual(values = c("#e74c3c", "#f1c40f")) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


p11 <- ds_long %>% filter(gene == "CRHR", Capture == "3") %>%
  ggplot(aes(y = perc_meth, x = as.factor(FromTSS), fill = Treatment )) +
  geom_point(aes(color = Treatment), position = position_jitterdodge(), alpha = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Location (relative to TSS)", y = "Percent methylation", color = "Treatment",
       title ="Post-manipulation") +
  scale_color_manual(values = c("#e74c3c", "#f1c40f"))  +
  scale_fill_manual(values = c("#e74c3c", "#f1c40f")) +
  theme(plot.title = element_text(size=15)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# cort plots  -------------------------------------------------------------


# Figure out how to plot this

ds_long %>% filter(Gene == "CRHR1") %>%
  ggplot(aes(x = base, y = perc_meth, color = Treatment )) +
  geom_point() +
  geom_smooth(method = "lm", aes(fill = Treatment)) +
  scale_color_manual(values = c("#e74c3c", "#f1c40f")) +
  scale_fill_manual(values = c("#e74c3c", "#f1c40f"  ))+
  labs(x = "Base cort", y = "Logit Methylation in CRHR1")

cortmethplot <- ds_long %>% filter(Gene == "CRHR1") %>%
  ggplot(aes(x = log(base + 1), y = perc_meth)) +
  geom_point(color = "#d35400" , alpha = 0.4) +
  geom_smooth(method = "lm", color ="#d35400" , fill = "#d35400") +
  labs(x = "Base cort (log scale)", y = "Percent methylation")





pdf("output_plots/cortxmethylation.pdf")
ds_long %>% filter(Gene == "CRHR1") %>%
  ggplot(aes(x = base, y = logit_meth, color = Capture)) +
  geom_point() +
  geom_smooth(method = "lm", aes(fill = Capture)) +
  labs(x = "Log(cort)", y = "Logit (Methylation) in CRHR1" ) +
  scale_color_manual(values = c("#3498db","#2ecc71"),
                     labels = c("Pre-treatment", "Post-treatment"))  +
  scale_fill_manual(values = c("#3498db","#2ecc71"), name = "Capture",
                    labels = c("Pre-treatment", "Post-treatment"))
dev.off()


# Plot1 and Plot2  --------------------------------------------------------



pdf("output_plots/fig1.pdf", height = 10, width = 14)
plot_allgenes / (p3 + p4 + p5) + plot_annotation(tag_levels = "A")
dev.off()

pdf("output_plots/fig1_alternative.pdf", height = 10, width = 14)
plot_allgenes / (p6 + p8 + p9) + plot_annotation(tag_levels = "A")
dev.off()

pdf("output_plots/fig2.pdf", height = 10, width = 12)
(p1 + p2) / (crhr_interaction_plot + cortmethplot ) + plot_annotation(tag_levels = 'A')
  dev.off()
pdf("output_plots/fig2_alternative.pdf", height = 10, width = 12 )
(p10 + p11)/(crhr_interaction_plot + cortmethplot ) + plot_annotation(tag_levels = 'A')
dev.off()


