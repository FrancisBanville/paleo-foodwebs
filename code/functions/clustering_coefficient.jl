"""
clustering_coefficient(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the average clustering coefficient across all species. The local clustering coefficient of a species corresponds to the proportion of possible links among its neighbours that are realized.
"""
function clustering_coefficient(N::UnipartiteNetwork)

    # get interaction list
    int = interactions(nodiagonal(N))
    pred = [int[i][1] for i in 1:links(N)]
    prey = [int[i][2] for i in 1:links(N)]

    # get species list
    sp_list = species(N)

    # calculate clustering coefficient for each species
    ci = zeros(Float64, length(sp_list))
   
    for (i, sp) in enumerate(sp_list)
        # find neighbours
        neighbours = unique(vcat(pred[prey .== sp], prey[pred .== sp]))
        # get subnetwork of neighbours
        N_neighbours = N[neighbours]
        # calculate the number of links among neighbours
        links_neighbours = links(N_neighbours)
        # calculate maximum number of links among neighbours
        k = length(neighbours)
        k_max = k * (k-1) / 2
        # calculate local clustering coefficient
        k_max == 0 ? ci[i] = 0 : ci[i] = links_neighbours  / k_max
    end

    # return mean clustering coefficient
    return mean(ci)
end
