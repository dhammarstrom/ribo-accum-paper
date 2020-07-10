#### qPCR analysis ########################






qpcrdat1 <- readRDS(file = "./data/derivedData/qpcr/qpcr_run1.RDS")

### Samples are loaded from extractions setup
samples <-  read_excel("./data/tr010_mRNASamples.xlsx", na = "NA") %>%
  mutate(sample = paste("S", ExtractionNR, sep = "")) %>%
  dplyr::select(participant, leg, time, sample, tissue_weight) %>%
  print()



### Combine data sets, qpcr, samples and leg randomizations
qpcr.dat <- qpcrdat1 %>%
  inner_join(samples) %>%
  inner_join(read_excel("./data/leg_randomization.xlsx")) %>%
  print()



qpcr.dat %>%
  ggplot(aes(target, cq, color = participant)) + geom_point() + coord_flip()



qpcr.dat %>%
  filter(target %in% c("rRNA28S F2R2", "rRNA18S F2R2")) %>%
  filter(cq > 20) %>%
  print()
