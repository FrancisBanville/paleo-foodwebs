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

# Read in data as needed
uncertain.data <- read.csv("results/uncertainty_analysis/dunne_method/uncertain.csv")
random.data <- read.csv("results/uncertainty_analysis/dunne_method/random.csv")
real.data <- read.csv("results/uncertainty_analysis/dunne_method/real.csv")
unc.roles.raw <- read.csv("results/uncertainty_analysis/removal_by_trophic_levels/uncertain_fez_roles.csv")
rand.roles.raw <- read.csv("results/uncertainty_analysis/removal_by_trophic_levels/random_fez_roles.csv")

# Function to manipulate data frames ----------------
# Write it:
perc <- c("0.1", "0.2", "0.5")

cleanup <- function(frame, percentages){
  frame$rep <- rep(1:100, 3)
  
  frame$perc <- rep(perc, each = 100, length_out = 300)
  
  frame <- pivot_longer(frame, cols = -c("rep", "perc"),
                        names_to = "Metric")
  
  return(frame)
}

# Run it:
uncertain <- cleanup(frame = uncertain.data, 
                     percentages = perc) 

random <- cleanup(frame = random.data, 
                  percentages = perc) 

# Replicate Dunne figs -----------------------
uncertain.means <- uncertain %>%
  group_by(perc, Metric) %>%
  summarise(avg = mean(value)) %>%
  mutate(type = "uncertain") %>%
  mutate(perc = as.numeric(perc))

random.means <- random %>%
  group_by(perc, Metric) %>%
  summarise(avg = mean(value)) %>%
  mutate(type = "random") %>%
  mutate(perc = as.numeric(perc))

means <- rbind(uncertain.means, random.means) %>%
  filter(!Metric %in% c("Can", "Loop"))

plotlist <- list()
for(i in 1:length(unique(means$Metric))){
  means.sub <- filter(means, Metric == unique(means$Metric)[i])
  
  plotlist[[i]] <- ggplot(data = means.sub, aes(x = perc, y = avg, 
                                                fill = type, 
                          color = type))+
    geom_point(size = 3, color = "black")+
    geom_point(aes(color = type),size = 2)+
    geom_line(color = "black")+
    # scale_shape_manual(values = 21:23, name = "Assemblage")+
    scale_fill_manual(values = c("black", "white"), name = "Type")+
    scale_color_manual(values = c("black", "white"), name = "Type")+
    labs(x = "Proportion Linkages Removed", 
         y = unique(means$Metric)[i])+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.title.x = element_blank())
}

dunne.plt <- wrap_plots(plotlist[1:17])+
  plot_layout(ncol = 5, guides = "collect")

# ggsave(filename = "figures/uncertainty_analysis/dunne_figs/dunne_repro.png", plot = dunne.plt,
#        width = 10, height = 8, units = "in", dpi = 600)

# Clean full network result ----------------
real_data <- real.data %>%
  pivot_longer(cols = everything(), names_to = "metric") %>%
  filter(!metric %in% c("Can", "Loop"))
  
# Make figures --------------------------
figs <- function(frames, metric){
  framelist <- list()
  
  for(i in 1:length(frames)){
    for(j in 1:length(metric)){
        smol.frame <- dplyr::filter(frames[[i]], Metric == metric[j])
        smol.frame$name <- paste(names(frames)[i],
                                     metric[j], sep = "_")
        smol.frame$fill_color <- "orange2"
        framelist <- append(framelist, list(smol.frame))
      }
    }

  plotlist <- list()
  
  for(i in 1:length(framelist)){
    tr <- real_data$value[which(real_data$metric ==
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

Top.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Top")

Bas.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Bas")

Int.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Int")

Herb.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Herb")

Omn.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Omn")

ChLen.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "ChLen")

ChSD.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "ChSD")

ChNum.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random),
                  metric = "ChNum")

TL.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                metric = "TL")

MxSim.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "MxSim")

VulSD.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "VulSD")

GenSD.figs <- figs(frames = list(uncertain = uncertain, 
                                 random = random), 
                   metric = "GenSD")

LinkSD.figs <- figs(frames = list(uncertain = uncertain, 
                                  random = random), 
                    metric = "LinkSD")

Path.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Path")

Clust.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Clust")

# Make 'em pretty ---------------
fig.grid <- function(figlist){
  # row_label_1 <- wrap_elements(panel = textGrob('Fezouata', rot=90))
  
  col_label_1 <- wrap_elements(panel = textGrob('Uncertain'))
  col_label_2 <- wrap_elements(panel = textGrob('Random'))
  
  big_ass_plot <- (figlist[[1]] | figlist[[2]])&
    plot_annotation(tag_levels = "a")
    # ((plot_spacer() / row_label_1) |
    #   (col_label_1 / figlist[[1]] / figlist[[2]] / figlist[[3]]) |
    #   (col_label_2 / figlist[[4]] / figlist[[5]] / figlist[[6]]))+
    # plot_layout(widths = c(0.5,1,1))
    # 
  return(big_ass_plot)
}

fig.grid(Bas.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/bas_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(ChLen.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/ChLen_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(ChNum.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/ChNum_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(GenSD.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/GenSD_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Herb.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/Herb_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Int.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/Int_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Omn.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/Omn_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(TL.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/TL_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Top.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/Top_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Clust.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/Clust_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(LinkSD.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/LinkSD_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(VulSD.figs)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs/VulSD_dists.png", width = 8, height = 6,
#        units = "in")

# Manipulate trophic data frames ----------------
cleanup_roles <- function(frame){
  frame$rep <- rep(1:100, 9)
  
  roles <- c("Basal", "Herb", "Omn")
  perc <- c("0.1", "0.2", "0.5")
  
  frame$perc <- rep(perc, each = 300, length_out = 900)
  frame$roles <- rep(roles, each = 100, length_out = 900)
  
  frame <- pivot_longer(frame, cols = -c("rep", "perc", "roles"),
                        names_to = "Metric")
  
  return(frame)
}

uncertain.roles <- cleanup_roles(frame = unc.roles.raw)

random.roles <- cleanup_roles(frame = rand.roles.raw) 

# Dunne fig trophic roles --------------
replicate_dunne <- function(datas, role){
  uncertain.means <- datas[[1]] %>%
    filter(roles == role) %>% 
    group_by(perc, Metric) %>%
    summarise(avg = mean(value)) %>%
    mutate(type = "uncertain") %>%
    mutate(perc = as.numeric(perc))
  
  random.means <- datas[[2]] %>%
    filter(roles == role) %>%
    group_by(perc, Metric) %>%
    summarise(avg = mean(value)) %>%
    mutate(type = "random") %>%
    mutate(perc = as.numeric(perc))
  
  means <- rbind(uncertain.means, random.means) %>%
    filter(!Metric %in% c("Can", "Loop"))
  
  plotlist <- list()
  for(i in 1:length(unique(means$Metric))){
    means.sub <- filter(means, Metric == unique(means$Metric)[i])
    
    plotlist[[i]] <- ggplot(data = means.sub, 
                            aes(x = perc, y = avg, fill = type,
                                color = type))+
      geom_point(size = 3, color = "black")+
      geom_point(aes(color = type),size = 2)+
      geom_line(color = "black")+
      scale_fill_manual(values = c("black", "white"), name = "Type")+
      scale_color_manual(values = c("black", "white"), name = "Type")+
      labs(x = "Proportion Linkages Removed", 
           y = unique(means$Metric)[i])+
      theme_bw()+
      theme(panel.grid = element_blank(), 
            axis.title.x = element_blank())
  }
  
  dunne.plt <- wrap_plots(plotlist[1:17])+
    plot_layout(ncol = 5, guides = "collect")
  
  return(dunne.plt)
}

basal.dunne <- replicate_dunne(datas = list(uncertain.roles, 
                                            random.roles),
                role = "Basal")
# ggsave(filename = "figures/uncertainty_analysis/dunne_figs/dunne_basal.jpeg", plot = basal.dunne, width = 10,
#        height = 8, units = "in", dpi = 600)

herb.dunne <- replicate_dunne(datas = list(uncertain.roles, 
                                            random.roles),
                               role = "Herb")
# ggsave(filename = "figures/uncertainty_analysis/dunne_figs/dunne_herb.jpeg", plot = herb.dunne, width = 10,
#        height = 8, units = "in", dpi = 600)

omn.dunne <- replicate_dunne(datas = list(uncertain.roles, 
                                           random.roles),
                              role = "Omn")
# ggsave(filename = "figures/uncertainty_analysis/dunne_figs/dunne_omn.jpeg", plot = omn.dunne, width = 10,
#        height = 8, units = "in", dpi = 600)

# Distributions: basal ------------------
# Filter basal removals
uncertain.basal <- uncertain.roles %>%
  filter(roles == "Basal")

random.basal <- random.roles %>%
  filter(roles == "Basal")

# Make figures and save
Bas.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                               random = random.basal), 
                 metric = "Bas")

fig.grid(Bas.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/bas_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

Herb.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                     random = random.basal), 
                       metric = "Herb")

fig.grid(Herb.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/herb_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

ChNum.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                     random = random.basal), 
                       metric = "ChNum")

fig.grid(ChNum.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/ChNum_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

Omn.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                       random = random.basal), 
                         metric = "Omn")

fig.grid(Omn.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/Omn_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

ChLen.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                     random = random.basal), 
                       metric = "ChLen")

fig.grid(ChLen.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/ChLen_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

VulSD.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                     random = random.basal), 
                       metric = "VulSD")

fig.grid(VulSD.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/VulSD_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

# Distributions: Herbivores ----------------------
# Filter herbivore removals
uncertain.herb <- uncertain.roles %>%
  filter(roles == "Herb")

random.herb <- random.roles %>%
  filter(roles == "Herb")

# Create & save figs
bas.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "Bas")

fig.grid(bas.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/bas_dists_herbroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

herb.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "Herb")

fig.grid(herb.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/herb_dists_herbroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

ChNum.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "ChNum")

fig.grid(ChNum.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/ChNum_dists_herbroles.png", width = 8,
#        height = 6, units = "in", dpi = 600)

Int.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "Int")

fig.grid(Int.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/int_dists_herbroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

VulSD.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "VulSD")

fig.grid(VulSD.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/VulSD_dists_herbroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

# Distributions: Omnivores -----------------------
# Filter omnivore removals
uncertain.omn <- uncertain.roles %>%
  filter(roles == "Omn")

random.omn <- random.roles %>%
  filter(roles == "Omn")

# Create & save plots
herb.figs.roles <- figs(frames = list(uncertain = uncertain.omn, 
                                      random = random.omn), 
                        metric = "Herb")

fig.grid(herb.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/herb_dists_omnroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

omn.figs.roles <- figs(frames = list(uncertain = uncertain.omn, 
                                      random = random.omn), 
                        metric = "Omn")

fig.grid(omn.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/omn_dists_omnroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

VulSD.figs.roles <- figs(frames = list(uncertain = uncertain.omn, 
                                     random = random.omn), 
                       metric = "VulSD")

fig.grid(VulSD.figs.roles)
# ggsave(filename = "figures/uncertainty_analysis/dist_figs_roles/VulSD_dists_omnroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)
