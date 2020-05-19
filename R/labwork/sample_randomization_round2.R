

#### Sample randomization round 2 ######

# Read data
tissue <- read_excel("./data/wetlab/tr010_tissue.xlsx", na = "NA") %>%

  print()


# Randomize samples for participant 19, 21, 22, 23
set.seed(1) # seed set to 1


extr_nr <- tissue %>%
  filter(participant %in% c("P19", "P21", "P22", "P23"), 
         sample %in% c("mRNA1", "mRNA2")) %>%
  dplyr::select(participant, leg, time, sample, samplenr, tissue_weight) %>%
  group_by(participant) %>%
  mutate(extraction_nr = sample(c(paste0("I:", seq(1:16)), paste0("II:", seq(1:16))), 
                                32, 
                                replace = FALSE)) %>%
  mutate(extraction_nr = factor(extraction_nr, levels = c(paste0("I:", seq(1:16)), paste0("II:", seq(1:16))))) %>%
  arrange(participant, extraction_nr) %>%
  print()
  



write_csv2(extr_nr, path = "./data/written/extraction_numbers_round2.csv")


