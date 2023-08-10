# Paleo Food Webs 

This repo contains the code to clean and analyze paleo food webs in Morocco. 

## Get started

You need to have Julia v1.8.0 or higher to run the code. 
Julia can be downloaded [here](https://julialang.org/downloads/). 

We recommend that you use an editor like [VS Code](https://code.visualstudio.com/docs/languages/julia) to execute the code. 

## Code organization 

`main.jl` executes the whole code from start to finish. The `Manifest.toml` and `Project.toml` files contain all information about package versions and dependencies. They are updated automatically when installing new packages. Do *not* edit them manually. 

Scripts that execute specific tasks are in the code folder.

`01_clean_data.jl` imports the raw dataset of species interactions (data/raw folder) and converts it into a clean unipartite network (saved in the data/clean folder). 

`02_simulate_networks.jl` predicts networks using the niche, cascade, and nested-hierarchy models with the number of species and interactions of the empirical network (saved in the data/sim folder).

`03_compute_measures.jl` calculates selected measures of network structure for both predicted and empirical networks (saved in the results folder).

`04_make_figures.jl` produces figures comparing the structure of empirical and predicted networks (saved in the figures folder).

