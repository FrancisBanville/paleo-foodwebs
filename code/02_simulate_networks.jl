## Simulate networks with the same number of species and interactions

# load clean network
N = load(joinpath("data", "clean", "network_Fezouata.jld2"))["N"]

# calculate number of species and interactions 
S = richness(N)
L = links(N)

# simulate networks using the niche, cascade, and nested-hierarchy models 
# n simulations for each model 

n = 1000

Ns_niche = []
Ns_cascade = []
Ns_nested_hierarchy = []

Random.seed!(1234)

for i in 1:n
    push!(Ns_niche, nichemodel(S, L))
    push!(Ns_cascade, cascademodel(S, L)), 
    push!(Ns_nested_hierarchy, nestedhierarchymodel(S, L))
end

save(joinpath("data", "sim", "networks_niche.jld2"), "Ns", Ns_niche)
save(joinpath("data", "sim", "networks_cascade.jld2"), "Ns", Ns_cascade)
save(joinpath("data", "sim", "networks_nested_hierarchy.jld2"), "Ns", Ns_nested_hierarchy)

