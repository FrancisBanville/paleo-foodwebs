## Compute selected measures of network structure for the empirical and predicted (simulated) networks

# load empirical networks
N_fezouata = load(joinpath("data", "clean", "network_Fezouata.jld2"))["N"]
N_fezouata_trophicsp = load(joinpath("data", "clean", "network_Fezouata_trophicsp.jld2"))["N"]


# load predicted (simulated) networks 
Ns_fezouata_sim = load(joinpath("data", "sim", "networks_fezouata_sim.jld2"))["Ns"]
Ns_fezouata_trophicsp_sim = load(joinpath("data", "sim", "networks_fezouata_trophicsp_sim.jld2"))["Ns"]


# group networks to facilitate calculations
Ns = vcat(N_fezouata, 
        N_fezouata_trophicsp, 
        Ns_fezouata_sim.Ns_niche,
        Ns_fezouata_sim.Ns_cascade,
        Ns_fezouata_sim.Ns_nested_hierarchy,
        Ns_fezouata_trophicsp_sim.Ns_niche,
        Ns_fezouata_trophicsp_sim.Ns_cascade,
        Ns_fezouata_trophicsp_sim.Ns_nested_hierarchy)

# simplify networks by removing isolated species
Ns = simplify.(Ns)

# make dataframe for all networks and measures

n = length(Ns_fezouata_sim.Ns_niche)

networks = vcat("fezouata",
                "fezouata trophic species",
                fill("fezouata", 3 * n), 
                fill("fezouata trophic species", 3 * n))


types = vcat(fill("empirical", 2), 
                fill("niche model", n), 
                fill("cascade model", n), 
                fill("nested hierarchy model", n), 
                fill("niche model", n), 
                fill("cascade model", n), 
                fill("nested hierarchy model", n))

measures = DataFrame(network = networks, type = types)


# calculate the number of species and links
S = richness.(Ns)
L = links.(Ns)



# calculate the proportion of species that are top (without consumers), intermediate (with both consumers and resources), and basal (without resources)

# Top (in-degree of 0)
@info "Calculating measure 1/17 (Top)"

kin = values.(degree.(Ns, dims = 2))
Top = sum.(x -> x == 0, kin) ./ S
insertcols!(measures, :Top => Top)

# Basal (out-degree of 0)
@info "Calculating measure 2/17 (Bas)"

kout = values.(degree.(Ns, dims = 1))
Bas = sum.(x -> x == 0, kout) ./ S
insertcols!(measures, :Bas => Bas)

# Int (proportion of species that are not Top or Basal)
@info "Calculating measure 3/17 (Int)"

Int = 1 .- Top .- Bas
insertcols!(measures, :Int => Int)


# calculate the proportion of species that are cannibals, herbivores (feeding only on basal species), omnivores (consuming two or more species with different trophic levels), and found in loops (food chains that contain the same species twice, apart from cannibalism)

# Cannibals (proportion of species interacting with itself)
@info "Calculating measure 4/17 (Can)"

Can = [sum(diag(Ns[i].edges)) for i in 1:length(Ns)] ./ S
insertcols!(measures, :Can => Can)

# Herbivores (proportion of species with a trophic level of 2)
@info "Calculating measure 5/17 (Herb)"

Herb = [sum((values(trophic_level(Ns[i])) .== 2)) for i in 1:length(Ns)] ./ S
insertcols!(measures, :Herb => Herb)

# Omnivores (proportion of species that consume two or more species and have food chains of different lengths)
@info "Calculating measure 6/17 (Omn)"

Omn = [sum(values(omnivory(Ns[i])) .> 0) for i in 1:length(Ns)] ./ S
insertcols!(measures, :Omn => Omn)

# Loops (proportion of species found in loops)
@info "Calculating measure 7/17 (Loop)"

Loop = zeros(Float64, length(Ns))

p = Progress(length(Ns))
Threads.@threads for i in 1:length(Ns)
        N = Ns[i]
        # remove self-loops
        N.edges[diagind(N.edges)] .= 0
        # proportion of species with a path to itself (without self-loops)
        Loop[i] = sum(diag(Matrix(shortest_path(N))) .> 0) / S[i]
next!(p)
end

insertcols!(measures, :Loop => Loop)



# calculate the average length of food chains, the standard deviation of their length, and the log number of food chains

# get food chain lengths
@info "Measuring food chain lengths"

food_chain_lengths = []

p = Progress(length(Ns))
Threads.@threads for i in 1:length(Ns)
        push!(food_chain_lengths, food_chains(Ns[i]))
next!(p)
end

# Average length of food chains
@info "Calculating measure 8/17 (ChLen)"

ChLen = mean.(food_chain_lengths)
insertcols!(measures, :ChLen => ChLen)

# Standard deviation of food-chain lengths
@info "Calculating measure 9/17 (ChSD)"

ChSD = std.(food_chain_lengths)
insertcols!(measures, :ChSD => ChSD)

# Log number of food chains
@info "Calculating measure 10/17 (ChNum)"

ChNum = log10.(length.(food_chain_lengths))
insertcols!(measures, :ChNum => ChNum)




# calculate the mean trophic level of all species 
# Trophic level of individual species is calculated with function trophic_level which is calculated with the formula
# TLᵢ = 1 + ∑ⱼ(TLⱼ×DCᵢⱼ)/∑ⱼ(DCᵢⱼ), wherein TLᵢ is the trophic level of species i, and DCᵢⱼ is the proportion of species j in 
# the diet of species i. (Poisot T, Bélisle Z, Hoebeke L, Stock M, Szefer P. EcologicalNetworks.jl: analysing ecological 
# networks of species interactions. Ecography. 2019;42(11):1850-1861. doi:https://doi.org/10.1111/ecog.04310)
# In the present case, it is assumed that each prey occupies an equal proportion of a consumer’s diet.
@info "Calculating measure 11/17 (TL)"

TL = [mean(values(trophic_level(Ns[i]))) for i in 1:length(Ns)]
insertcols!(measures, :TL => TL)



# calculate the average of the maximum trophic similarity of each species
@info "Calculating measure 12/17 (MxSim)"

MxSim = zeros(Float64, length(Ns))

p = Progress(length(Ns))
Threads.@threads for i in 1:length(Ns)
        MxSim[i] = MaxSim(Ns[i])
next!(p)
end

insertcols!(measures, :MxSim => MxSim)



# calculate the normalized standard deviations of vulnerability (nb of consumers or in-degree), generality (nb of resources or out-degree), and total links (nb of consumers and resources or total degree)
# species in, out, and total degrees are normalized by the average number of interactions per species (2L/S)

# Vulnerability
@info "Calculating measure 13/17 (VulSD)"

VulSD = vulnerability.(Ns)
insertcols!(measures, :VulSD => VulSD)

# Generality
@info "Calculating measure 14/17 (GenSD)"

GenSD = generality.(Ns)
insertcols!(measures, :GenSD => GenSD)

# Total links
@info "Calculating measure 15/17 (LinkSD)"

LinkSD = total_links.(Ns)
insertcols!(measures, :LinkSD => LinkSD)



# calculate the average shortest food-chain length between all pairs of species

# Average shortest path (not taking into account unconnected pairs)
@info "Calculating measure 16/17 (Path)"

Path = zeros(Float64, length(Ns))

p = Progress(length(Ns))
Threads.@threads for i in 1:length(Ns)
        paths = shortest_path(Ns[i])
        Path[i] = mean(paths[Not(paths .== 0)])
next!(p)
end

insertcols!(measures, :Path => Path)



# calculate the mean clustering coefficient, the probability that two species linked to the same species are also linked 

# Mean clustering coefficient
@info "Calculating measure 17/17 (Clust)"

Clust = clustering_coefficient.(Ns)
insertcols!(measures, :Clust => Clust)


# remove NaNs and missing values
measures[isnan.(measures.ChLen),:ChLen] .= 0
measures[isnan.(measures.ChSD),:ChSD] .= 0
measures[ismissing.(measures.ChNum),:ChNum] .= 0

# remove infinite values that were obtained after taking the log of 0 
measures[measures.ChNum .== -Inf, :ChNum] .= 0

# export table
CSV.write(joinpath("results", "ecological_models", "measures.csv"), measures)

## Calculate model error 

# create empty dataset for predictive errors
measures_errors = DataFrame()
    
# calculate predictive errors of niche model for all networks and measures 
for N in unique(measures.network) 

    # get empirical measures of network N
    measures_emp = measures[measures.network .== N .&& measures.type .== "empirical", 3:end]

    # get niche model predictions for network N
    measures_niche = measures[measures.network .== N .&& measures.type .== "niche model", 3:end]

    # calculate predictive errors for network N
    measures_errors_N = (measures_niche .- measures_emp) ./ measures_emp
    
    # add network name column 
    measures_errors_N.network = fill(N, size(measures_errors_N, 1))

    # add to general dataframe
    measures_errors = [measures_errors; measures_errors_N]
end

# remove NaNs and missing values
measures_errors[isnan.(measures_errors.Can),:Can] .= Inf
measures_errors[isnan.(measures_errors.Loop),:Loop] .= Inf

# export table
CSV.write(joinpath("results", "ecological_models", "measures_errors.csv"), measures_errors)
