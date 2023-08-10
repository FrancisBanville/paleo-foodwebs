## Compute measures of empirical and predicted networks

# load networks

N = load(joinpath("data", "clean", "network_Fezouata.jld2"))["N"]

Ns_niche = load(joinpath("data", "sim", "networks_niche.jld2"))["Ns"]
Ns_cascade = load(joinpath("data", "sim", "networks_cascade.jld2"))["Ns"]
Ns_nested_hierarchy = load(joinpath("data", "sim", "networks_nested_hierarchy.jld2"))["Ns"]

# simplify networks by removing isolated species
N = simplify(N)
Ns_niche = simplify.(Ns_niche)
Ns_cascade = simplify.(Ns_cascade)
Ns_nested_hierarchy = simplify.(Ns_nested_hierarchy)

# group networks to facilitate calculations
Ns = vcat(N, 
        Ns_niche, 
        Ns_cascade, 
        Ns_nested_hierarchy)

# make dataframe for all networks and measures

n = length(Ns_niche)

types = vcat("empirical",
                fill("niche model", n), 
                fill("cascade model", n),
                fill("nested hierarchy model", n))

measures = DataFrame(type = types)


# calculate the proportion of species that are top (without consumers), intermediate (with both consumers and resources), and basal (without resources)

# Top 
# Int
# Bas

# calculate the proportion of species that are cannibals, herbivores (feeding only on basal species), omnivores (consuming two or more species with different trophic levels), and found in loops (food chains that contain the same species twice, apart from cannibalism)

# Can
# Herb
# Omn
# Loop

# calculate the average length of food chains, the standard deviation of their length, and the log number of food chains

# ChLen
# ChSD
# ChNum

# calculate the mean trophic level of all species computed using the short-weighted trophic level algorithm

# TL

# calculate the average of the maximum trophic similarity of each species

# MaxSim

# calculate the normalized standard deviations of vulnerability (nb of consumers or in-degree), generality (nb of resources or out-degree), and total links (nb of consumers and resources or total degree)
# species in, out, and total degrees were normalized by the average number of interactions per species (2L/S)

# VulSD
# GenSD
# LinkSD

# calculate the average shortest food-chain length between all pairs of species

# Path

# calculate the mean clustering coefficient, the probability that two species linked to the same species are also linked 

# Clust


