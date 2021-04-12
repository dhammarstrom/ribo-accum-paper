############# Figure 3 ##################################


## Panels

# Use Bayesian estimates?
bayes <- TRUE


source("./R/figure-source.R")
# source("./R/libs.R")


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
  
  geom_segment(data = primers, aes(y = c(1.54, 1.52, 1.528, 1.54, 1.55, 1.56), 
                                   yend = c(1.54, 1.52, 1.528, 1.54, 1.55, 1.56), 
                                   x= on45S_start, xend = on45S_end), 
               size = 4) +
  
  geom_text(data = primers, aes(y = c(1.54, 1.52, 1.528, 1.544, 1.55, 1.56),
                                      x= on45S_end, label = symbol_ext), 
                  size = 2.2, 
          hjust = 0,
            position = position_nudge(x = 200)) +
  
  # Annotations 
  annotate("text", x = c(3655/2, 5523 - (5523-3655)/2, 
                         6757 - (6757-6601)/2,
                         12990 - (12990-7925)/2), 
                   y = rep(1.48, 4), 
           color = "gray50",
           label = c("47/45S ETS", "18S", "5.8S", "28S"), 
           fontface = "italic",
  
           size = 2.2) +
  
  
  
  scale_y_continuous(limits = c(1.46, 1.58)) +
  plot_theme() + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_blank())



###### Panel X Intervention vs. control qPCR rRNA   ##################

qpcr_res_int_con <- readRDS("./data/derivedData/qpcr-analysis-bayes2/qpcr_res_int_con.RDS")


complete_rrna_comp <- qpcr_res_int_con %>%
          filter(model == "tissue") %>%
          dplyr::select(target, comparison = contrast, estimate, lower.CL , upper.CL) %>%

  filter(!(target %in% c("UBTF F4R4",
                        "UBTF F6R6", 
                        "rpS6 F2R2", 
                        "MyHC1 F1R1", 
                        "MyHC2A F5R5", 
                        "MyHC2X F5R5"))) %>%
  print()

 




# Create an annotation data frame

anno.df <- data.frame(target = unique(complete_rrna_comp$target), 
           label = c("rRNA 18S",
                     "rRNA 28S",  
                     "rRNA 45S ETS",
                     "rRNA 45S ITS", 
                     "rRNA 47S ETS", 
                     "rRNA 5.8S",
                     "rRNA 5S" ), 
           time = factor("S1", levels = c("S1", "post")), 
           estimate = 3.5) %>%
  mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                            "rRNA5.8S F2R2",
                                            "rRNA28S F2R2",
                                            "rRNA5S F3R3", 
                                            "rRNA45SITS F12R12",
                                            "rRNA45S F5R5",     
                                            "rRNA47S F1R1")), 
         comparison = "S1",
         comparison = factor(comparison, levels = c("S1", "post", "post1w"), 
                             labels = c("Session 1", 
                                        "Post-training", 
                                        "Post-trainin \n+ De-training")),
         target = fct_rev(target))



interaction_effects <- complete_rrna_comp %>%
  

  filter(comparison %in% c("inter:S1", "inter:post", "inter:post1w","")) %>%
  
  mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                            "rRNA5.8S F2R2",
                                            "rRNA28S F2R2",
                                            "rRNA5S F3R3", 
                                            "rRNA45SITS F12R12",
                                            "rRNA45S F5R5",     
                                            "rRNA47S F1R1")),
         
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
  
  geom_text(data = anno.df %>%
              mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                                        "rRNA5.8S F2R2",
                                                        "rRNA28S F2R2",
                                                        "rRNA5S F3R3", 
                                                        "rRNA45SITS F12R12",
                                                        "rRNA45S F5R5",     
                                                        "rRNA47S F1R1"))),
            aes(estimate - 2, target,  label = label, color = NULL), 
            position = position_nudge(y = 0.25), 
            hjust = 0,
            size = 2.2) +
  
  
  labs(x = "Fold change compared to Control") +
  
  geom_vline(xintercept = 1, color = "gray85", lty = 2) +
  
  geom_point() +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2) + 
  facet_wrap(  comparison ~ .) +
  
  scale_color_manual(values = c("gray50", "gray10")) +
  
  plot_theme() +
  theme(strip.background = element_rect(color = "white", fill = "white"), 
        strip.text = element_text(size = 7),
        axis.title.y = element_blank()  , 
        axis.text.y = element_blank(), 
        legend.position = "none")
  
  


fold_changes <- complete_rrna_comp %>%
  
  filter(!(comparison %in% c("inter:S1", "inter:post", "inter:post1w",""))) %>%
  separate(comparison, into = c("group", "time"), sep = "_") %>%
  
  
  
  mutate(target = factor(target, levels = c(
                                            "rRNA18S F2R2",
                                            "rRNA5.8S F2R2",
                                            "rRNA28S F2R2",
                                            "rRNA5S F3R3", 
                                            "rRNA45SITS F12R12",
                                            "rRNA45S F5R5",     
                                            "rRNA47S F1R1")),
         target = fct_rev(target),

         estimate = exp(estimate), 
         lower.CL = exp(lower.CL), 
         upper.CL = exp(upper.CL), 
         group = if_else(time == "post1w", "int_detrain", group),
         
         group = factor(group, levels = c("con", "int", "int_detrain"), 
                        labels = c("Control", "Training", "Training\n+De-training")),
         
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
  
  scale_y_continuous(limits = c(0.5, 3.8), breaks = c(1, 2, 3)) +
  scale_x_discrete(limits = c( "S1", "post"), labels = c("Session 1", "Post-\ntraining")) +
  
  scale_fill_manual(values = c(color.scale[5], color.scale[6],color.scale[8])) +
  
  
  labs(y = "Fold change from Baseline") +
  plot_theme() +
  
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size = 7, margin = margin(t = 0.1, b= 0.1,r = 0.1,l = 0.1, unit = "pt")),
        legend.key.size = unit(0, "cm"),
        legend.margin = margin(-0.5, 0, -0.5, 0, "cm"),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.direction = "vertical",
        axis.title.x = element_blank()) +
  
  facet_grid(target ~ .)






### time course rrna ##############

# "trace plots" of each qpcr-estimate


rrna_data <- readRDS("./data/derivedData/qpcr-analysis-bayes2/time_course_qpcr.RDS")

# Get Estimates and create a data set for statistical robustness indicators
# Statistical robust is defined as when 0 is not within the 95% credible interval


rrna_timecourse <- rrna_data$estimated_means %>%
  filter(model == "tissue", 
         target %in% c("rRNA18S F2R2",
                       "rRNA5.8S F2R2",
                       "rRNA28S F2R2",
                       "rRNA5S F3R3", 
                       "rRNA45SITS F12R12",
                       "rRNA45S F5R5",     
                       "rRNA47S F1R1")) %>%
  
  mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                            "rRNA5.8S F2R2",
                                            "rRNA28S F2R2",
                                            "rRNA5S F3R3", 
                                            "rRNA45SITS F12R12",
                                            "rRNA45S F5R5",     
                                            "rRNA47S F1R1")),
                target = fct_rev(target), 
         time.c = gsub("S", "", time), 
         time.c = as.numeric(if_else(time.c == "post1w", "12", time.c))) %>%
  print()


diffs <- rrna_data$estimated_diff %>%
  
  filter(model == "tissue") %>%
  
  mutate(target = factor(target, levels = c("rRNA18S F2R2",
                                            "rRNA5.8S F2R2",
                                            "rRNA28S F2R2",
                                            "rRNA5S F3R3", 
                                            "rRNA45SITS F12R12",
                                            "rRNA45S F5R5",     
                                            "rRNA47S F1R1")),
         target = fct_rev(target), 
         time.c = gsub("S", "", time), 
         time.c = as.numeric(if_else(time.c == "post1w", "12", time.c))) %>%
  filter(robust == "robust") %>%
  mutate(cond = "var") %>%
  left_join(rrna_timecourse %>% 
              group_by(target, time) %>%
              summarise(emmean = max(emmean))) %>% 
  mutate(emmean = emmean * 1.05, 
         sign = "\U0002A") %>%
  print()


targets <- unique(rrna_timecourse$target)


plots <- list()


for(i in 1:length(targets)) {
  
  dat <- rrna_timecourse %>%
    filter(target == targets[i]) 
   
  
  
  diffs2 <- diffs %>%
    filter(target == targets[i])
      
  plot  <-   dat %>%
    filter(time != "post1w" )  %>%
    ggplot(aes(time.c, emmean, color = cond)) + 
    geom_line() +
    geom_text(data = diffs2, 
              aes(time.c, emmean, label = sign, color = NULL), 
              show.legend = FALSE) +
    
    geom_point(data = filter(dat, time == "post1w"), 
               position = position_nudge(x = 1)) +
    
    scale_color_manual(values = c(color.scale[1], color.scale[2])) +
    
    scale_y_continuous(limits = c(floor(min(dat$emmean, na.rm = TRUE)), 
                                  ceiling(max(dat$emmean, na.rm = TRUE))),
                       breaks = c(floor(min(dat$emmean, na.rm = TRUE)),
                                  ceiling(max(dat$emmean, na.rm = TRUE)) + ((floor(min(dat$emmean, na.rm = TRUE)) - ceiling(max(dat$emmean, na.rm = TRUE)))/2), 
                       ceiling(max(dat$emmean, na.rm = TRUE)))) +
    
    labs(x = "Session", 
         y = "Estimated log-abundance per tissue weight (AU)")  +
    plot_theme() +
    
    ggtitle(anno.df[anno.df$target == targets[i], 2]) + 
    
    theme(strip.background = element_blank(), 
          strip.text = element_blank(), 
          plot.title = element_text(size = 7),

          legend.position = "none", 
          axis.title = element_blank(), 
          axis.text.x = element_blank(), 
          axis.line.x = element_blank(), 
          axis.ticks.x = element_blank()) 
  

  plots[[i]]  <- plot
  
}


names(plots) <- unique(rrna_timecourse$target)
 


## Make axis label for subplots
library(grid)
library(gridExtra)

axis_lab <- textGrob("Estimated log-counts per tissue weight (AU)", 
                     gp=gpar( fontsize=7),
          rot=90)

### create a plot with only x axis

x_axis <- data.frame(Session = c(0, 3, 6, 9, 12), 
           y       = c(0, 0, 0, 0, 0)) %>%
  ggplot(aes(Session, y)) + geom_blank() + 
  plot_theme() + 
  plot_theme() +
  
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12), 
                     limits = c(0, 12.5)) +
  
  labs(x = "Session") +
  
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank()) 



time_traces <-  plot_grid(axis_lab, 
           plot_grid(plots[[5]], 
          plots[[3]], 
          plots[[4]], 
          plots[[7]], 
          plots[[2]], 
          plots[[6]], 
          plots[[1]], 
          x_axis,
          
          rel_heights = c(rep(1, 7), 0.5),
          
          align = "v",
          
          ncol = 1), ncol = 2, rel_widths = c(0.1, 1))




rnra_ctrl_vs_int <- cowplot::plot_grid(
  cowplot::plot_grid(NULL, fold_changes, NULL, ncol = 1, rel_heights = c(0.09, 1, 0.06)),
  interaction_effects,

  ncol = 2, rel_widths = c(0.6, 1))


# Total RNA analysis and plot #################################################

# load model comparisons between 

if(bayes == TRUE){
  
  comp_rna <- readRDS("./data/derivedData/total-rna-analysis/comp_rna_bayes.RDS")
  
}

if(bayes == FALSE){
  
  comp_rna <- readRDS("./data/derivedData/total-rna-analysis/comp_rna_freq.RDS")
  
}


tot_rna_interaction <- comp_rna %>%
  mutate(comparison = time_group, 
         target = "totalRNA") %>%
  dplyr::select(target, comparison, estimate:upper.CL) %>%
  filter(comparison %in% c("inter:S1", "inter:post", "inter:post1w","")) %>%
  
  mutate(
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
  
  scale_color_manual(values = c("gray50", "gray10")) +
  scale_x_continuous(limits = c(0.8, 2), 
                     expand = c(0, 0), 
                     breaks = c(0.8, 1, 1.2, 1.4, 1.6, 1.8, 2), 
                     labels = c("",  1, "",   1.4, "", 1.8, "")) +
  
  plot_theme() +
  theme(strip.background = element_rect(color = "white", fill = "white"), 
        strip.text = element_text(size = 7),
        axis.title.y = element_blank()  , 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        legend.position = "none")


  
tot_rna_fold_change <- comp_rna %>%
  mutate(comparison = time_group, 
         target = "totalRNA") %>%
  dplyr::select(target, comparison, estimate:upper.CL) %>%

  filter(!(comparison %in% c("inter:S1", "inter:post", "inter:post1w",""))) %>%
  separate(comparison, into = c("group", "time"), sep = "_") %>%

  mutate(  estimate = exp(estimate), 
    lower.CL = exp(lower.CL), 
    upper.CL = exp(upper.CL), 
    group = if_else(time == "post1w", "int_detrain", group),
    time = if_else(time == "post1w", "post", time),
    time = factor(time, levels = c("S1", "post"))) %>%
  
  ggplot(aes(time, estimate, fill = group)) + 
  
  geom_hline(yintercept = 1, lty = 2, color = "gray80") +
  
  # geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.15) + 
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                position = position_dodge(width = 0.3), 
                width = 0) + 
  
  geom_point(position = position_dodge(width = 0.3), shape = 21) +
  
  scale_x_discrete(limits = c( "S1", "post"), labels = c("Session 1", "Post-\ntraining")) +
  scale_y_continuous(limits = c(0.8, 2), breaks = c(1, 1.2, 1.4, 1.6, 1.8, 2), 
                   labels = c(1, 1.2, 1.4, 1.6, 1.8, 2), 
                   expand = c(0, 0)) +
  
  scale_fill_manual(values = c(color.scale[5], color.scale[6],color.scale[8])) +
  labs(y = "Fold change from Baseline") +
  plot_theme() +
  
  theme(strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank()) 











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
        y = bquote('Total RNA (mg'~'ng'^-1*')')) +
   
   scale_x_continuous(limits = c(0, 12), expand = c(0,0), 
                      breaks = c(0, 3, 6, 9, 12)) +
   
   scale_y_continuous(limits = c(250, 650), expand = c(0,0), 
                      breaks = c(250, 300, 350, 400, 450, 500, 550, 600, 650), 
                      labels = c("",  300, "",  400, "",  500, "" , 600, "")) +
  
  scale_color_manual(values = c(color.scale[2], color.scale[1]), guide = NULL) +
  scale_fill_manual(values = c(color.scale[2], color.scale[1])) +
   
   plot_theme() +
   theme(legend.position = c(0.7, 0.15), 
         legend.key.size = unit(0.3, "cm"), 
         axis.title.y = element_text() ) 
   

   
tot_rna_fig <- plot_grid(tot_rna_fold_change, 
                                   tot_rna_interaction, 
                                   rna_tc_fig, 
                                   ncol = 3,
                     
          rel_widths = c(1, 1, 1))


############## Combine panels  ###################



                  
figure3 <- plot_grid( 
              plot_grid(
                plot_grid(
                    plot_grid(NULL, primers_fig,NULL,ncol = 3, rel_widths = c(0.2, 2, 0.2)),
                    rnra_ctrl_vs_int, ncol = 1, rel_heights = c(0.2, 1)), 
              time_traces, ncol = 2, rel_widths = c(0.7, 0.3)),
              
           tot_rna_fig, 
           rel_heights = c(1.2, 0.5),
           ncol = 1) +
  
  draw_plot_label(label=c("A",  "B",   "C",  "D", "E",   "F",  "G"),
                  x =   c(0.02, 0.02, 0.27, 0.75, 0.02, 0.33, 0.71), 
                  y =   c(0.98, 0.83, 0.83, 0.98, 0.3, 0.27, 0.3),
                  hjust=.5, vjust=.5, size = label.size)





saveRDS(figure3, "./Figures/rds/figure3.RDS")

# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure3.pdf", plot = figure3, width = 8.9*2, height = 23 * 0.75, 
       dpi = 600,
       units = "cm", device=cairo_pdf)

