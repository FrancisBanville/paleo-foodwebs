"""
vulnerability(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the normalized standard deviations of vulnerability (nb of consumers)
"""
function vulnerability(N::UnipartiteNetwork)
    # calculate the number of species and interactions
    S = richness(N)
    L = links(N)

    # calculate the in-degree of each species
    kin = degree(N; dims=2)
    kin = [kin[s] for s in species(N)]

    # return the normalized standard deviation of vulnerability
    return  std(kin) / (2*L / S)
end

"""
generality(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the normalized standard deviations of generality (nb of resources)
"""
function generality(N::UnipartiteNetwork)
    # calculate the number of species and interactions
    S = richness(N)
    L = links(N)

    # calculate the out-degree of each species
    kout = degree(N; dims=1)
    kout = [kout[s] for s in species(N)]

    # return the normalized standard deviation of generality
    return  std(kout) / (2*L / S)
end


"""
total_links(N::UnipartiteNetwork)
    N: Unipartite simple directed network
Returns the normalized standard deviations of total links (nb of consumers and resources)
"""
function total_links(N::UnipartiteNetwork)
    # calculate the number of species and interactions
    S = richness(N)
    L = links(N)

    # calculate the total degree of each species
    ktot = degree(N)
    ktot = [ktot[s] for s in species(N)]

    # return the normalized standard deviation of total links
    return  std(ktot) / (2*L / S)
end
