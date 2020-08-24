

#### Figure X RNA/qPCR ###################



## Panels
# A: rRNA primers 





### rRNA and primers #########


# Sequences can be found : https://www.ncbi.nlm.nih.gov/nuccore/NR_146144.1
# Add primer sequences to plot and annotate mature segments



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
  annotate("segment", x = 7925, xend = 12990, y = 1.5, yend = 1.5, size = 2.5) 






