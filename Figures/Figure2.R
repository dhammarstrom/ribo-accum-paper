############# Figure 2 ##################################


## Panels

# Use Bayesian estimates?
bayes <- TRUE


source("./R/figure-source.R")
source("./R/libs.R")


# 45S primers graphics #############################

### rRNA and primers #########


# Sequences can be found : https://www.ncbi.nlm.nih.gov/nuccore/NR_146144.1
# Add primer sequences to plot and annotate mature segments

# Primer locations
primers <- read_excel("./data/tr010_primers.xlsx") %>%
  filter(primer_id %in% c("rRNA47S F1R1", 
                          "rRNA45S F5R5", 
                          "rRNA45SITS F12R12", 
                          "rRNA5.8S F2R2", 
                          "rRNA28S F2R2", 
                          "rRNA18S F2R2")) %>%
  print()



primers_fig <- data.frame(x = c(1, 13500), 
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
  
  geom_segment(data = primers, aes(y = c(1.54, 1.52, 1.53, 1.54, 1.55, 1.56), 
                                   yend = c(1.54, 1.52, 1.53, 1.54, 1.55, 1.56), 
                                   x= on45S_start, xend = on45S_end), 
               size = 4) +
  
  geom_text_repel(data = primers, aes(y = c(1.54, 1.52, 1.53, 1.54, 1.55, 1.56),
                                      x= on45S_start, label = symbol_ext), 
                  size = 2.2, 
                  nudge_x = 100) +
  
  scale_y_continuous(limits = c(1.48, 1.58)) +
  plot_theme() + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank())




# Total RNA analysis and plot #################################################

# load model comparisons between 

if(bayes == TRUE){
  
comp_rna <- readRDS("./data/derivedData/total-rna-analysis/comp_rna_bayes.RDS")
  
}

if(bayes == FALSE){
  
  comp_rna <- readRDS("./data/derivedData/total-rna-analysis/comp_rna_freq.RDS")
  
}






###### Panel X Intervention vs. control qPCR rRNA + total RNA  ##################


qpcr_res_int_con <- readRDS("./data/derivedData/qpcr-analysis-bayes/qpcr_res_int_con.RDS")


complete_rrna_comp <- comp_rna %>%
  mutate(comparison = time_group, 
         target = "totalRNA") %>%
  dplyr::select(target, comparison, estimate:upper.CL) %>%
  rbind(qpcr_res_int_con %>%
          dplyr::select(target, comparison, estimate = Estimate, lower.CL = CI.Lower, upper.CL = CI.Upper)) %>%
  print()



# Create an annotation data frame

anno.df <- data.frame(target = unique(complete_rrna_comp$target), 
           label = c("Total RNA", 
                     "rRNA 18S", 
                     "rRNA 45S ETS", 
                     "rRNA 47S ETS", 
                     "rRNA 5S", 
                     "rRNA 28S", 
                     "rRNA 45S ITS", 
                     "rRNA 5.8S"), 
           time = factor("S1", levels = c("S1", "post")), 
           estimate = 3.5) %>%
  mutate(target = factor(target, levels = c("totalRNA", 
                                            "rRNA18SF2R2",
                                            "rRNA5.8SF2R2",
                                            "rRNA28SF2R2",
                                            "rRNA5SF3R3", 
                                            "rRNA45SITSF12R12",
                                            "rRNA45SF5R5",     
                                            "rRNA47SF1R1")),
         target = fct_rev(target))



interaction_effects <- complete_rrna_comp %>%
  

  
  
  filter(comparison %in% c("inter:S1", "inter:post", "inter:post1w","")) %>%
  
  mutate(target = factor(target, levels = c("totalRNA", 
                                            "rRNA18SF2R2",
                                            "rRNA5.8SF2R2",
                                            "rRNA28SF2R2",
                                            "rRNA5SF3R3", 
                                            "rRNA45SITSF12R12",
                                            "rRNA45SF5R5",     
                                            "rRNA47SF1R1")),
         
         comparison = gsub("inter:", "", comparison), 
         comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                             labels = c("Session 1", 
                                        "Post-training", 
                                        "Post-trainin \n+ De-training")),
         estimate = exp(estimate), 
         lower.CL = exp(lower.CL), 
         upper.CL = exp(upper.CL), 
         robust = if_else(lower.CL  > 1, "robust", "notrobust")) %>%
  
  ggplot(aes(estimate, 
             target, color = robust)) + 
  
  
  labs(x = "Fold change compared to Control") +
  
  geom_vline(xintercept = 1, color = "gray85", lty = 2) +
  
  geom_point() +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) + 
  facet_wrap(  comparison ~ .) +
  
  scale_color_manual(values = c("blue", "red")) +
  
  plot_theme() +
  theme(strip.background = element_rect(color = "white", fill = "white"), 
        strip.text = element_text(size = 8),
        axis.title.y = element_blank()  , 
        axis.text.y = element_blank(), 
        legend.position = "none")
  
  


fold_changes <- complete_rrna_comp %>%
  
  filter(!(comparison %in% c("inter:S1", "inter:post", "inter:post1w",""))) %>%
  separate(comparison, into = c("group", "time"), sep = "_") %>%
  
  
  
  mutate(target = factor(target, levels = c("totalRNA", 
                                            "rRNA18SF2R2",
                                            "rRNA5.8SF2R2",
                                            "rRNA28SF2R2",
                                            "rRNA5SF3R3", 
                                            "rRNA45SITSF12R12",
                                            "rRNA45SF5R5",     
                                            "rRNA47SF1R1")),
         target = fct_rev(target),

         estimate = exp(estimate), 
         lower.CL = exp(lower.CL), 
         upper.CL = exp(upper.CL), 
         group = if_else(time == "post1w", "int_detrain", group),
         time = if_else(time == "post1w", "post", time),
         time = factor(time, levels = c("S1", "post"))) %>%
  
  ggplot(aes(time, estimate, fill = group)) + 
  
  geom_text(data = anno.df, aes(time, estimate, label = label, fill = NULL), 
            position = position_nudge(x = -0.5), 
            hjust = 0,
            size = 2.2) +
  geom_hline(yintercept = 1, lty = 2, color = "gray80") +

  # geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.15) + 
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                position = position_dodge(width = 0.3), 
                width = 0) + 
  
  geom_point(position = position_dodge(width = 0.3), shape = 21) +
  
  scale_x_discrete(limits = c( "S1", "post"), labels = c("Session 1", "Post-\ntraining")) +
  
  labs(y = "Fold change from Baseline") +
  plot_theme() +
  
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank()) +
  
  facet_grid(target ~ .)


rnra_ctrl_vs_int <- cowplot::plot_grid(
  cowplot::plot_grid(NULL, fold_changes, NULL, ncol = 1, rel_heights = c(0.08, 1, 0.06)),
  interaction_effects, ncol = 2, rel_widths = c(0.6, 1))



# Time course in the intervention group #####

cond_eff_rna_tc <- readRDS("./data/derivedData/total-rna-analysis/cond_eff_rna_tc.RDS")



rna_tc_fig <-  cond_eff_rna_tc %>%
  mutate(cond = factor(cond, levels = c("const", "var"), 
                       labels = c("Constant volume", 
                                  "Variable volume")), 
         cond = fct_rev(cond)) %>%
  ggplot(aes(time.c, estimate, color = cond, fill = cond)) +
  geom_line() + 
   geom_ribbon(aes(ymin = lower, ymax = upper, color = NULL), alpha = 0.2) +
   
   labs(x = "Session", 
        y = "RNA ng \U00D7 mg<sup>-1</sup>") +
   
   scale_x_continuous(limits = c(0, 12), expand = c(0,0), 
                      breaks = c(0, 3, 6, 9, 12)) +
   
   scale_y_continuous(limits = c(250, 650), expand = c(0,0), 
                      breaks = c(250, 300, 350, 400, 450, 500, 550, 600, 650), 
                      labels = c("",  300, "",  400, "",  500, "" , 600, "")) +
   
   plot_theme() +
   theme(legend.position = c(0.5, 0.15), 
         legend.key.size = unit(0.3, "cm"), 
         axis.title.y = element_markdown(), 
         ) 
   

   














############## Combine panels  ###################



figure2 <- plot_grid( 
                     plot_grid(plot_grid(primers_fig, rna_tc_fig, ncol = 1, rel_heights = c(0.5, 1.5)),
                               rnra_ctrl_vs_int, rel_widths = c(0.8, 1),
                               ncol = 2),
                     plot_grid(rna_tc_fig, NULL, ncol = 2), 
                     ncol = 1, rel_heights = c(1, 0.5)) +
  
  
  draw_plot_label(label=c("A",  "B",  "C", "D"),
                  x =   c(0.02, 0.45, 0.5, 0.75), 
                  y =   c(0.98, 0.98, 0.33, 0.33),
                  hjust=.5, vjust=.5, size = label.size)






# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure2.pdf", plot = figure2, width = 8.9*2, height = 23 * 0.75, 
       dpi = 600,
       units = "cm", device=cairo_pdf)

