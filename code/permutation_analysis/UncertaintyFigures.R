################################################
# Paleo Food Webs Project: figures for         #
# uncertainty analysis                         #
# Dec 2023                                     #
# Code by E.M. Beasley                         #
# In collab. with F. Banville, C. Soucy, G.    #
# Ramirez Guerro, & C. Cameron                 #
################################################

# Load packages ---------------------------
library(tidyverse)
library(patchwork)
library(grid)

# Load data -----------------------------
setwd("c:/users/beasl/documents/paleo-foodwebs/code/permutation_analysis")

uncertain.data <- read.csv("uncertain.csv")
random.data <- read.csv("random.csv")
real.data <- read.csv("real.csv")

# Function to manipulate data frames ----------------
cleanup <- function(frame){
  frame$rep <- rep(1:100, 9)
  
  assemblage <- c("Fezouata", "Burgess", "Chengjiang")
  perc <- c("0.1", "0.2", "0.5")
  
  frame$assemblage <- rep(assemblage, each = 300)
  frame$perc <- rep(perc, each = 100, length_out = 900)
  
  frame <- pivot_longer(frame, cols = -c("rep", "assemblage", "perc"),
                        names_to = "Metric")
  
  return(frame)
}

uncertain <- cleanup(frame = uncertain.data)
random <- cleanup(frame = random.data)

# Replicate Dunne figs -----------------------
uncertain.means <- uncertain %>%
  group_by(assemblage, perc, Metric) %>%
  summarise(avg = mean(value)) %>%
  mutate(type = "uncertain") %>%
  mutate(perc = as.numeric(perc))

random.means <- random %>%
  group_by(assemblage, perc, Metric) %>%
  summarise(avg = mean(value)) %>%
  mutate(type = "random") %>%
  mutate(perc = as.numeric(perc))

means <- rbind(uncertain.means, random.means)

plotlist <- list()
for(i in 1:length(unique(means$Metric))){
  means.sub <- filter(means, Metric == unique(means$Metric)[i])
  
  plotlist[[i]] <- ggplot(data = means.sub, aes(x = perc, y = avg, 
                                                fill = type,
                          shape = assemblage, color = type))+
    geom_point(size = 3, color = "black")+
    geom_point(aes(color = type),size = 2)+
    geom_line(color = "black")+
    scale_shape_manual(values = 21:23, name = "Assemblage")+
    scale_fill_manual(values = c("black", "white"), name = "Type")+
    scale_color_manual(values = c("black", "white"), name = "Type")+
    labs(x = "Proportion Linkages Removed", 
         y = unique(means$Metric)[i])+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.title.x = element_blank())
}

dunne.plt <- wrap_plots(plotlist[1:18])+
  plot_layout(ncol = 6, guides = "collect")

# ggsave(filename = "dunne_repro.png", plot = dunne.plt, width = 15, 
#        height = 8, units = "in", dpi = 600)

# Clean full network result ----------------
real.data$site = c("Fezouata", "Burgess", "Chengjiang")

real_data <- real.data %>%
  pivot_longer(cols = -c("site"), names_to = "metric")

# Make figures --------------------------
# Issue in this function:
figs <- function(frames, metric, site){
  framelist <- list()
  
  for(i in 1:length(frames)){
    for(j in 1:length(metric)){
      for(k in 1:length(site)){
        smol.frame <- dplyr::filter(frames[[i]], Metric == metric[j] &
                                 assemblage == site[k])
        smol.frame$name <- paste(names(frames)[i], site[k],
                                     metric[j], sep = "_")
        smol.frame$fill_color <- 
          case_when(site[k] == "Fezouata" ~ "orange2",
                    site[k] == "Burgess" ~ "steelblue1",
                    TRUE ~ "Palevioletred")
        framelist <- append(framelist, list(smol.frame))
      }
    }
  }
  
  plotlist <- list()
  
  for(i in 1:length(framelist)){
    tr <- real_data$value[which(real_data$site == 
                            unique(framelist[[i]]$assemblage) &
                          real_data$metric ==
                            unique(framelist[[i]]$Metric))]
    
    plotlist[[i]] <- ggplot(data = framelist[[i]], 
                            aes(x = value, alpha = perc))+
      geom_density(color = "black", 
                   fill = unique(framelist[[i]]$fill_color))+
      geom_vline(xintercept = tr, color = "red", linetype = "dashed")+
      labs(x = unique(framelist[[i]]$Metric), y = "Density")+
      theme_bw()+
      theme(panel.grid = element_blank())
  }
  
  return(plotlist)
}

unique(uncertain$Metric)

Top.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Top", site = unique(uncertain$assemblage))

Bas.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Bas", site = unique(uncertain$assemblage))

Int.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Int", site = unique(uncertain$assemblage))

Can.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Can", site = unique(uncertain$assemblage))

Herb.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Herb", site = unique(uncertain$assemblage))

Omn.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Omn", site = unique(uncertain$assemblage))

Loop.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Loop", site = unique(uncertain$assemblage))

ChLen.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "ChLen", 
                   site = unique(uncertain$assemblage))

ChSD.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "ChSD", site = unique(uncertain$assemblage))

ChNum.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random),
                  metric = "ChNum", site = unique(uncertain$assemblage))

TL.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                metric = "TL", site = unique(uncertain$assemblage))

MxSim.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "MxSim", 
                   site = unique(uncertain$assemblage))

VulSD.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "VulSD", 
                   site = unique(uncertain$assemblage))

GenSD.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "GenSD", 
                   site = unique(uncertain$assemblage))

LinkSD.figs <- figs(frames = list(uncertain = uncertain, 
                                  random = random), 
                    metric = "LinkSD", 
                    site = unique(uncertain$assemblage))

Path.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Path", site = unique(uncertain$assemblage))

# Make 'em pretty ---------------
fig.grid <- function(figlist){
  row_label_1 <- wrap_elements(panel = textGrob('Fezouata', rot=90))
  row_label_2 <- wrap_elements(panel = textGrob('Burgess', rot=90))
  row_label_3 <- wrap_elements(panel = textGrob('Chengiang', rot=90))
  
  col_label_1 <- wrap_elements(panel = textGrob('Uncertain'))
  col_label_2 <- wrap_elements(panel = textGrob('Random'))
  
  big_ass_plot <- 
    ((plot_spacer() / row_label_1 / row_label_2 / row_label_3) |
      (col_label_1 / figlist[[1]] / figlist[[2]] / figlist[[3]]) |
      (col_label_2 / figlist[[4]] / figlist[[5]] / figlist[[6]]))+
    plot_layout(widths = c(0.5,1,1))
  
  return(big_ass_plot)
}

fig.grid(Bas.figs)
# ggsave(filename = "bas_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(ChLen.figs)
# ggsave(filename = "ChLen_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(ChNum.figs)
# ggsave(filename = "ChNum_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(GenSD.figs)
# ggsave(filename = "GenSD_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Herb.figs)
# ggsave(filename = "Herb_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Int.figs)
# ggsave(filename = "Int_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Omn.figs)
# ggsave(filename = "Omn_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(TL.figs)
# ggsave(filename = "TL_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Top.figs)
# ggsave(filename = "Top_dists.png", width = 8, height = 6,
#        units = "in")

