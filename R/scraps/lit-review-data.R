

temp <- read_excel("./review/rna-literature-review.xlsx") %>%
  filter(subject != "all") %>%
  print()





## Digitize 
# library(digitize)
# bickel <- digitize("./review/digitize/Bickel2005.tif")




read_excel("./review/rna-literature-review.xlsx") %>%
  filter(statistics != "individual") %>%
  group_by(subject, study_name, n_sessions) %>%
  summarise(rna_fc = mean(rna_fc, na.rm = TRUE)) %>%
  arrange(n_sessions, study_name) %>%
  ggplot(aes(n_sessions, rna_fc, group = paste(study_name, subject))) + 
  geom_point() + geom_line() 

