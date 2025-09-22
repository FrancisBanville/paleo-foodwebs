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
# Change working directory as needed
# setwd("c:/users/beasl/documents/paleo-foodwebs/code/permutation_analysis")

# Read in data as needed
uncertain.data <- read.csv("uncertain.csv")
random.data <- read.csv("random.csv")
real.data <- read.csv("real.csv")
unc.roles.raw <- read.csv("uncertain_roles.csv")
rand.roles.raw <- read.csv("random_roles.csv")

# Function to manipulate data frames ----------------
# Write it:
assemblage <- c("Fezouata", "Burgess", "Chengjiang")
perc <- c("0.1", "0.2", "0.5")

cleanup <- function(frame, assemblages, percentages){
  frame$rep <- rep(1:100, 9)
  
  frame$assemblage <- rep(assemblage, each = 300)
  frame$perc <- rep(perc, each = 100, length_out = 900)
  
  frame <- pivot_longer(frame, cols = -c("rep", "assemblage", "perc"),
                        names_to = "Metric")
  
  return(frame)
}

# Run it:
uncertain <- cleanup(frame = uncertain.data, assemblages=assemblage,
                     percentages = perc) %>%
  # Include this line for Fezouata only:
  filter(assemblage == "Fezouata")

random <- cleanup(frame = random.data, assemblages = assemblage,
                  percentages = perc) %>%
  filter(assemblage == "Fezouata")

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

means <- rbind(uncertain.means, random.means) %>%
  filter(!Metric %in% c("Can", "Loop"))

plotlist <- list()
for(i in 1:length(unique(means$Metric))){
  means.sub <- filter(means, Metric == unique(means$Metric)[i])
  
  plotlist[[i]] <- ggplot(data = means.sub, aes(x = perc, y = avg, 
                                                fill = type,
                          #shape = assemblage, 
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

# ggsave(filename = "dunne_repro.png", plot = dunne.plt,
#        width = 10, height = 8, units = "in", dpi = 600)

# Clean full network result ----------------
real.data$site = c("Fezouata", "Burgess", "Chengjiang")

real_data <- real.data %>%
  pivot_longer(cols = -c("site"), names_to = "metric") %>%
  filter(!metric %in% c("Can", "Loop")) %>%
  filter(site=="Fezouata")

# Make figures --------------------------
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
          case_when(site[k] == "Fezouata" ~ "orange2"#,
                    # site[k] == "Burgess" ~ "steelblue1",
                    # TRUE ~ "Palevioletred"
                    )
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

Top.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Top", site = unique(uncertain$assemblage))

Bas.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Bas", site = unique(uncertain$assemblage))

Int.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Int", site = unique(uncertain$assemblage))

Herb.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Herb", site = unique(uncertain$assemblage))

Omn.figs <- figs(frames = list(uncertain = uncertain, random = random), 
                 metric = "Omn", site = unique(uncertain$assemblage))

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

Clust.figs <- figs(frames = list(uncertain = uncertain, random = random),
                  metric = "Clust", site = unique(uncertain$assemblage))

# Make 'em pretty ---------------
fig.grid <- function(figlist){
  # row_label_1 <- wrap_elements(panel = textGrob('Fezouata', rot=90))
  # row_label_2 <- wrap_elements(panel = textGrob('Burgess', rot=90))
  # row_label_3 <- wrap_elements(panel = textGrob('Chengiang', rot=90))
  
  col_label_1 <- wrap_elements(panel = textGrob('Uncertain'))
  col_label_2 <- wrap_elements(panel = textGrob('Random'))
  
  big_ass_plot <- (figlist[[1]] | figlist[[2]])&
    plot_annotation(tag_levels = "a")
    # ((plot_spacer() / row_label_1 / row_label_2 / row_label_3) |
    #   (col_label_1 / figlist[[1]] / figlist[[2]] / figlist[[3]]) |
    #   (col_label_2 / figlist[[4]] / figlist[[5]] / figlist[[6]]))+
    # plot_layout(widths = c(0.5,1,1))
    # 
  return(big_ass_plot)
}

fig.grid(Bas.figs)
# ggsave(filename = "./dist_figs/bas_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(ChLen.figs)
# ggsave(filename = "./dist_figs/ChLen_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(ChNum.figs)
# ggsave(filename = "./dist_figs/ChNum_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(GenSD.figs)
# ggsave(filename = "./dist_figs/GenSD_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Herb.figs)
# ggsave(filename = "./dist_figs/Herb_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Int.figs)
# ggsave(filename = "./dist_figs/Int_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Omn.figs)
# ggsave(filename = "./dist_figs/Omn_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(TL.figs)
# ggsave(filename = "./dist_figs/TL_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Top.figs)
# ggsave(filename = "./dist_figs/Top_dists.png", width = 8, height = 6,
#        units = "in")

fig.grid(Clust.figs)
# ggsave(filename = "./dist_figs/Clust_dists.png", width = 8, height = 6,
#        units = "in")

# Manipulate trophic data frames ----------------
cleanup_roles <- function(frame){
  frame$rep <- rep(1:100, 27)
  
  assemblage <- c("Fezouata", "Burgess", "Chengjiang")
  roles <- c("Basal", "Herb", "Omn")
  perc <- c("0.1", "0.2", "0.5")
  
  frame$assemblage <- rep(assemblage, each = 900)
  frame$perc <- rep(perc, each = 300, length_out = 2700)
  frame$roles <- rep(roles, each = 100, length_out = 2700)
  
  frame <- pivot_longer(frame, cols = -c("rep", "assemblage", "perc",
                                         "roles"),
                        names_to = "Metric")
  
  return(frame)
}

uncertain.roles <- cleanup_roles(frame = unc.roles.raw) %>%
  filter(assemblage=="Fezouata")

random.roles <- cleanup_roles(frame = rand.roles.raw) %>%
  filter(assemblage=="Fezouata")

# Dunne fig trophic roles --------------
replicate_dunne <- function(datas, role){
  uncertain.means <- datas[[1]] %>%
    filter(roles == role) %>% 
    group_by(assemblage, perc, Metric) %>%
    summarise(avg = mean(value)) %>%
    mutate(type = "uncertain") %>%
    mutate(perc = as.numeric(perc))
  
  random.means <- datas[[2]] %>%
    filter(roles == role) %>%
    group_by(assemblage, perc, Metric) %>%
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
                                # shape = assemblage, 
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
# ggsave(filename = "dunne_basal.jpeg", plot = basal.dunne, width = 10,
#        height = 8, units = "in", dpi = 600)

herb.dunne <- replicate_dunne(datas = list(uncertain.roles, 
                                            random.roles),
                               role = "Herb")
# ggsave(filename = "dunne_herb.jpeg", plot = herb.dunne, width = 10,
#        height = 8, units = "in", dpi = 600)

omn.dunne <- replicate_dunne(datas = list(uncertain.roles, 
                                           random.roles),
                              role = "Omn")
# ggsave(filename = "dunne_omn.jpeg", plot = omn.dunne, width = 10,
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
                 metric = "Bas", site = unique(uncertain.basal$assemblage))

fig.grid(Bas.figs.roles)
# ggsave(filename = "./dist_figs_roles/bas_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

Herb.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                     random = random.basal), 
                       metric = "Herb", 
                       site = unique(uncertain.basal$assemblage))

fig.grid(Herb.figs.roles)
# ggsave(filename = "./dist_figs_roles/herb_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

ChNum.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                     random = random.basal), 
                       metric = "ChNum", 
                       site = unique(uncertain.basal$assemblage))

fig.grid(ChNum.figs.roles)
# ggsave(filename = "./dist_figs_roles/ChNum_dists_basroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

Omn.figs.roles <- figs(frames = list(uncertain = uncertain.basal, 
                                       random = random.basal), 
                         metric = "Omn", 
                         site = unique(uncertain.basal$assemblage))

fig.grid(Omn.figs.roles)
# ggsave(filename = "./dist_figs_roles/Omn_dists_basroles.png", width = 8, height = 6,
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
                       metric = "Bas", 
                       site = unique(uncertain.herb$assemblage))

fig.grid(bas.figs.roles)
# ggsave(filename = "./dist_figs_roles/bas_dists_herbroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

herb.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "Herb", 
                       site = unique(uncertain.herb$assemblage))

fig.grid(herb.figs.roles)
# ggsave(filename = "./dist_figs_roles/herb_dists_herbroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

ChNum.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "ChNum", 
                       site = unique(uncertain.herb$assemblage))

fig.grid(ChNum.figs.roles)
# ggsave(filename = "./dist_figs_roles/ChNum_dists_herbroles.png", width = 8,
#        height = 6, units = "in", dpi = 600)

Int.figs.roles <- figs(frames = list(uncertain = uncertain.herb, 
                                     random = random.herb), 
                       metric = "Int", 
                       site = unique(uncertain.herb$assemblage))

fig.grid(Int.figs.roles)
# ggsave(filename = "./dist_figs_roles/int_dists_herbroles.png", width = 8, height = 6,
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
                        metric = "Herb", 
                        site = unique(uncertain.omn$assemblage))

fig.grid(herb.figs.roles)
# ggsave(filename = "./dist_figs_roles/herb_dists_omnroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)

omn.figs.roles <- figs(frames = list(uncertain = uncertain.omn, 
                                      random = random.omn), 
                        metric = "Omn", 
                        site = unique(uncertain.omn$assemblage))

fig.grid(omn.figs.roles)
# ggsave(filename = "./dist_figs_roles/omn_dists_omnroles.png", width = 8, height = 6,
#        units = "in", dpi = 600)
