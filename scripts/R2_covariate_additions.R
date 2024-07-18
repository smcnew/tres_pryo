# McNew et al. 2024
# Manipulation of a social signal affects DNA methylation of a
# stress-related gene in a free-living bird
#
# Journal of Experimental Biology R2
#
# Script to compile additional covariates on female age and reproductive effort
# in response to reviewer's question.

library(dplyr)


# Load data ---------------------------------------------------------------

# Existing covariates
covars <- read.delim("raw_data/Input_Cort_Fled.txt") %>% rename(Band = fband) %>%
  mutate(Capture = as.factor(Capture)) %>% filter(year == "2017")


# Nest data
nestdat <- read.csv("unshared_data/Nest_Records 11.18.2020.csv") %>%
          filter(Exp_Year == "2017") %>%
          filter(Location == "Ithaca") %>%
          rename(Band = Female_ID) %>%
          select(Band, Clutch_Size, Attempt_Number, Site, Nest, Nest_Experiment)
covars <- left_join(covars, nestdat)

# Female capture data
femdat <- read.csv("unshared_data/Captures_Hormone_Bleeding_Blood_DNA 11.18.2020.csv") %>%
  filter(Exp_Year == "2017") %>%
  filter(Location == "Ithaca") %>%
  filter(Adult_or_Nestling == "Adult") %>%
  filter(Sex == "F") %>%
  rename(Band = Individual_Band) %>%
  filter(Capture_Number == 1 ) %>%
  select(Band, Age, Experiment)

covars <- left_join(covars, femdat)

write.csv(covars, "raw_data/covars_updated.csv", row.names = F)



# Sample sizes for captures -----------------------------------------------

femdat2 <- read.csv("unshared_data/Captures_Hormone_Bleeding_Blood_DNA 11.18.2020.csv") %>%
  filter(Exp_Year == "2017") %>%
  filter(Location == "Ithaca") %>%
  filter(Adult_or_Nestling == "Adult") %>%
  filter(Sex == "F") %>%
  rename(Band = Individual_Band) %>%
  filter(Experiment == "Color")
table(femdat2$Capture_Number)
filter(femdat2, Capture_Number == 3) %>% filter(!(Band %in% covars$Band)) %>% dim
121-70
70-55

