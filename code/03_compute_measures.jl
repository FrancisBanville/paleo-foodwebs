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


# calculate the number of species and links
S = richness.(Ns)
L = links.(Ns)

# calculate the proportion of species that are top (without consumers), intermediate (with both consumers and resources), and basal (without resources)

# Top (in-degree of 0)
kin = values.(degree.(Ns, dims = 2))
Top = sum.(x -> x == 0, kin) ./ S
insertcols!(measures, :Top => Top)

# Basal (out-degree of 0)
kout = values.(degree.(Ns, dims = 1))
Bas = sum.(x -> x == 0, kout) ./ S
insertcols!(measures, :Bas => Bas)

# Int (proportion of species that are not Top or Basal)
Int = 1 .- Top .- Bas
insertcols!(measures, :Int => Int)

# calculate the proportion of species that are cannibals, herbivores (feeding only on basal species), omnivores (consuming two or more species with different trophic levels), and found in loops (food chains that contain the same species twice, apart from cannibalism)

# Cannibals (proportion of species interacting with itself)
Can = [sum(diag(Ns[i].edges)) for i in 1:length(Ns)] ./ S
insertcols!(measures, :Can => Can)

# Herbivores (proportion of species with a trophic level of 2)
Herb = [sum((values(trophic_level(Ns[i])) .== 2)) for i in 1:length(Ns)]
insertcols!(measures, :Herb => Herb)

# Omnivores (proportion of species that consume two or more species and have food chains of different lengths)
Omn = [sum(values(omnivory(Ns[i])) .> 0) for i in 1:length(Ns)] ./ S
insertcols!(measures, :Omn => Omn)

# Loops (proportion of species found in loops)
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

# ChLen
# ChSD
# ChNum

# calculate the mean trophic level of all species 
TL = [mean(values(trophic_level(Ns[i]))) for i in 1:length(Ns)]
insertcols!(measures, :TL => TL)

# calculate the average of the maximum trophic similarity of each species
MxSim = zeros(Float64, length(Ns))

p = Progress(length(Ns))
Threads.@threads for i in 1:length(Ns)
        MxSim[i] = MaxSim(Ns[i])
next!(p)
end

insertcols!(measures, :MxSim => MxSim)

# calculate the normalized standard deviations of vulnerability (nb of consumers or in-degree), generality (nb of resources or out-degree), and total links (nb of consumers and resources or total degree)
# species in, out, and total degrees are normalized by the average number of interactions per species (2L/S)

# VulSD
# GenSD
# LinkSD

# calculate the average shortest food-chain length between all pairs of species

# Path

# calculate the mean clustering coefficient, the probability that two species linked to the same species are also linked 

# Clust


### Measures done: 9/17



