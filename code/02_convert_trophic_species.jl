### Convert Anticosti network to a network of trophic species
### A trophic species is a group of species having the same predators and prey 

# load clean network
N_anticosti = load(joinpath("data", "clean", "network_anticosti.jld2"))["N"]

# get adjacency matrix 
anticosti_mat = Matrix(N_anticosti.edges)

# convert adjacency matrix to wide format to more easily identify species with identical predators and prey 
S = length(species(N_anticosti))
anticosti_mat_wide = zeros(Int64, S, 2 * S)

# the first half of a row corresponds to a species' predators
# the second half of a row corresponds to a species' prey
for s in 1:S
    pred_s = anticosti_mat[:,s]
    prey_s = anticosti_mat[s,:]
    anticosti_mat_wide[s,:] = vcat(pred_s, prey_s)
end

# identify groups of species having the same predators and prey
# each group (trophic species) has a different number 
groups = groupslices(anticosti_mat_wide, dims=1)

# find indices of groups of more than 1 species
groups_ind = groupinds(groups)
ind = filter(is -> length(is)> 1, groups_ind)

# remove first index of each of these groups so that we can keep exactly one species in each group later
indx = reduce(vcat, [ind[i][2:end] for i in 1:length(ind)])

# remove duplicate species 
anticosti_mat_trophicsp = anticosti_mat[Not(indx), Not(indx)]

# convert to unipartite network and save
N_anticosti_trophicsp = simplify(UnipartiteNetwork(anticosti_mat_trophicsp))

save(joinpath("data", "clean", "network_anticosti_trophicsp.jld2"), "N", N_anticosti_trophicsp)
