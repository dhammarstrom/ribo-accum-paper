######## Western blot analysis analysis ####################

# Filename: ./R/wester_analysis.R
# Author: Daniel H
# Notes:

# Purpose: Analyze western blot data

# Output: Model estimates

######## -----------------------####################

## Packages ####
source("./R/libs.R")


### Load data  round 1 ####

western_data_round1 <- read_excel("./data/wetlab/western/tr010_western_round1.xlsx", na = "NA", sheet = "total_protein") %>%
  filter(participant != "LADDER") %>%
  dplyr::select(participant, ExtractionNR, round, gel, well, mean.gray1:total.protein2) %>%
  inner_join(read_excel("./data/wetlab/tr010_mRNASamples_round1.xlsx", na = "NA") %>%
               dplyr::select(participant, leg, time, ExtractionNR)) %>%
  inner_join(read_excel("./data/wetlab/western/tr010_western_round1.xlsx", na = "NA", sheet = "ecl") %>%
               dplyr::select(round, gel, well, target, signal)) %>%
  mutate(gel = paste0(round,"_", gel), 
         sample = paste0(participant, "_", ExtractionNR)) %>%
  rowwise() %>%
  ### memcode stain - background = total protein ### 
  mutate(total.protein = mean(c(total.protein1, total.protein2)) - mean(c(mean.gray1, mean.gray2, mean.gray3))) %>%
  ungroup() %>%
  dplyr::select(participant,sample, leg, time, gel, well, target, total.protein, signal, round) %>%
  group_by(gel) %>%
  ### Normalizes the total protein stain per 
  mutate(tp.factor = total.protein) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  dplyr::select(participant,sample, leg, cond, sex, time, gel, well, target, tp.factor, signal, round) %>%
  ungroup() %>%
  mutate(expression = signal / tp.factor) %>%
  filter(target %in% c("t-s6", "t-UBF")) %>%
  print()


western_data_round2 <- read_excel("./data/wetlab/western/tr010_western_round2.xlsx", na = "NA", sheet = "total_protein") %>%

  dplyr::select(participant, ser, num, round, gel, well, total.protein_1:meangray_4) %>%
  mutate(num = as.character(num)) %>%
  inner_join(read_csv2("./data/written/extraction_numbers_round2.csv", na = "NA") %>%
               dplyr::select(participant, leg, time, extraction_nr) %>%
               separate(extraction_nr, into = c("ser", "num"), sep = ":") %>%
               mutate(ser = as.numeric(if_else(ser == "I", 1, 2)))) %>%
  inner_join(read_excel("./data/wetlab/western/tr010_western_round2.xlsx", na = "NA", sheet = "ecl") %>%
               dplyr::select(round, gel, well, target = Image, signal = Signal)) %>%
  mutate(target = gsub("_.*", "", target), 
         target = if_else(target == "rpS6", "t-s6", "t-UBF"),
         gel = paste0(round,"_", gel), 
         sample = paste0(participant, "_", ser, "_", num)) %>%
  rowwise() %>%
  ### memcode stain - background = total protein ### 
  mutate(total.protein = mean(c(total.protein_1, total.protein_2)) - mean(c(meangray_1, meangray_2, meangray_3, meangray_4))) %>%
  ungroup() %>%
  dplyr::select(participant, sample, leg, time, gel, well, target, total.protein, signal, round) %>%
  group_by(gel) %>%
  ### Normalizes the total protein stain per 
  mutate(tp.factor = total.protein) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  dplyr::select(participant, sample, leg, cond, sex, time, gel, well, target, tp.factor, signal, round) %>%
  ungroup() %>%
  mutate(expression = signal / tp.factor) %>%
  filter(target %in% c("t-s6", "t-UBF")) %>%
  print()




# Per gel normalization:
# Due to differences in memcode/ecl between gels expression is normalized to 
# gel averages

western_data <- rbind(western_data_round1, western_data_round2) %>%
  group_by(gel, target) %>%
  mutate(tp.factor = tp.factor / max(tp.factor, na.rm = TRUE), 
         signal = signal / max(signal, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(expression = signal / tp.factor) %>%
  dplyr::select(participant,sample, leg, cond, time, gel, target, expression) %>%
  print()




# Save data 
saveRDS(western_data, "./data/derivedData/western-compile/western_data.RDS")



