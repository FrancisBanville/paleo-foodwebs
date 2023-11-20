## Simulate networks with the same number of species and interactions

# load clean networks
N_anticosti = load(joinpath("data", "clean", "network_anticosti.jld2"))["N"]
N_anticosti_trophicsp = load(joinpath("data", "clean", "network_anticosti_trophicsp.jld2"))["N"]

# simulate networks using the niche, cascade, and nested-hierarchy models for each empirical network
# n simulations for each model 

Random.seed!(1234)

function simulate_networks(N::UnipartiteNetwork; n=100) 
    
    # calculate number of species and interactions 
    S = richness(N)
    L = links(N)

    Ns_niche = []
    Ns_cascade = []
    Ns_nested_hierarchy = []

    for i in 1:n
        push!(Ns_niche, nichemodel(S, L))
        push!(Ns_cascade, cascademodel(S, L)), 
        push!(Ns_nested_hierarchy, nestedhierarchymodel(S, L))
    end
    return(Ns_niche = Ns_niche, Ns_cascade = Ns_cascade, Ns_nested_hierarchy = Ns_nested_hierarchy)
end


Ns_anticosti_sim = simulate_networks(N_anticosti)
Ns_anticosti_trophicsp_sim = simulate_networks(N_anticosti_trophicsp)

save(joinpath("data", "sim", "networks_anticosti_sim.jld2"), "Ns", Ns_anticosti_sim)
save(joinpath("data", "sim", "networks_anticosti_trophicsp_sim.jld2"), "Ns", Ns_anticosti_trophicsp_sim)

