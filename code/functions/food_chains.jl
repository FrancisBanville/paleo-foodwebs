"""
food_chains(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns a vector of food chain lengths computed from all shortest paths from top predators to basal species.
"""
function food_chains(N::UnipartiteNetwork)
    # find all shortest paths between every pair of species 
    paths = dijkstra(nodiagonal(N))

    # keep only the paths between top predators and basal species
    basal = collect(keys(filter(p -> iszero(p.second), degree(N; dims=1))))
    top = collect(keys(filter(p -> iszero(p.second), degree(N; dims=2))))

    filter!(int -> int.to ∈ basal, paths)
    filter!(int -> int.from ∈ top, paths)

    # get food chain lengths 
    chain_lengths = [int.weight for int in paths] 

    # return summary statistics
    return chain_lengths
end

