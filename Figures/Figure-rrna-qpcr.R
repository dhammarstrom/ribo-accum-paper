

#### Figure X RNA/qPCR ###################



## Panels
# A: rRNA primers 





### rRNA and primers #########


# Sequences can be found : https://www.ncbi.nlm.nih.gov/nuccore/NR_146144.1
# Add primer sequences to plot and annotate mature segments

# Primer locations
primers <- read_excel("./data/tr010_primers.xlsx") %>%
  filter(primer_id %in% c("47S F1R1", 
                          "45S F5R5", 
                          "45SITS F12R12", 
                          "5.8S F2R2", 
                          "28S F2R2", 
                          "18S F2R2")) %>%
  print()



data.frame(x = c(1, 13500), 
           y = c(1, 2)) %>%
  ggplot(aes(x, y)) +
  
  # Full 45S segment 
  annotate("segment", x = 1, xend = 13351, y = 1.5, yend = 1.5, size = 0.5) +
  # 18S segment
  annotate("segment", x = 3655, xend = 5523, y = 1.5, yend = 1.5, size = 2.5) +
  # 5.8S segment
  annotate("segment", x = 6601, xend = 6757, y = 1.5, yend = 1.5, size = 2.5) +
  # 28S segment
  annotate("segment", x = 7925, xend = 12990, y = 1.5, yend = 1.5, size = 2.5) +
  
  geom_segment(data = primers, aes(y = c(1.51, 1.52, 1.53, 1.54, 1.55, 1.56), 
                                   yend = c(1.51, 1.52, 1.53, 1.54, 1.55, 1.56), 
                                   x= on45S_start, xend = on45S_end), 
               size = 4) +
  scale_y_continuous(limits = c(1.48, 1.58)) +
  plot_theme() + 
  theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank())






