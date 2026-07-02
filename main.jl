## Activate project environment
## Creates a list of packages and dependencies with specific versions used
## for reproducibility reasons
## See Project.toml and Manifest.toml files
import Pkg; Pkg.activate(".")

Pkg.instantiate()

## Load required packages
## Use Pkg.add("Package Name") to add other packages

# Manipulating variables, data frames and files
import CSV 
using DataFrames 
using GroupSlices
using InvertedIndices
using JLD2

# Doing statistics and models
using Distributions
using LinearAlgebra 
using ProgressMeter
using Random
using SparseArrays
using Statistics
using StatsBase

# Making plots
using Plots 
using StatsPlots 

# Analyzing ecological networks
using EcologicalNetworks 

## Load custom functions
include(joinpath("code", "functions", "clustering_coefficient.jl"))
include(joinpath("code", "functions", "food_chains.jl"))
include(joinpath("code", "functions", "MaxSim.jl"))
include(joinpath("code", "functions", "relative_degree.jl"))

## Load scripts
include(joinpath("code", "01_clean_data.jl"))
include(joinpath("code", "02_convert_trophic_species.jl"))
include(joinpath("code", "03_simulate_networks.jl"))
include(joinpath("code", "04_compute_measures.jl"))
include(joinpath("code", "05_make_figures.jl"))
include(joinpath("code", "06_create_uncertain_networks.jl"))
include(joinpath("code", "07_make_uncertainty_figures.R"))
