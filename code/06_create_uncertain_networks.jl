################################################
# Paleo Food Webs Project                      #
# Dec 2023                                     #
# Code by E.M. Beasley                         #
# In collab. with F. Banville, C. Soucy, G.    #
# Ramirez Guerro, & C. Cameron                 #
################################################

## Create networks with different levels of interaction uncertainty (uncertainty analysis)

# read raw datasets

fezouata_df = DataFrame(CSV.File(joinpath("data", "raw", "interactions_Fezouata.csv")))

burgess_df = DataFrame(CSV.File(joinpath("data", "raw", "interactions_burgess_avec_incertitude.csv")))

chengjiang_df = DataFrame(CSV.File(joinpath("data", "raw", "interactions_chengjiang_avec_incertitude.csv")))


# clean datasets 

function clean_data(df::DataFrame)
    
    # rename variables 
    rename!(df, 1 => :pred, 2 => :prey, 3 => :certainty)

    # remove empty rows since they do not represent interactions
    dropmissing!(df)

    if eltype(df[:,1]) !== Int64
        # save certainty levels
        certainty = df[:,3]

        # convert uppercase letters to lowercase
        df = lowercase.(df[:,1:2])

        # remove symbols that artificially creates new species when inconsistent 
        df = replace.(df[:,1:2], "?" => "")
        df = replace.(df[:,1:2], "'" => "")
        df = replace.(df[:,1:2], "\"" => "")

        # remove leading and trailing white spaces
        df = strip.(df) 

        df.certainty = certainty
    end
    
    # remove duplicate rows
    df = unique(df)

    return(df)

end

fezouata_df_clean = clean_data(fezouata_df)
burgess_df_clean = clean_data(burgess_df)
chengjiang_df_clean = clean_data(chengjiang_df)



### Uncertainty Analysis: Dunne Method

# Takes a few hours but can be run overnight

# Function to remove n% of uncertain values
function remove_uncertain(dat, proportion)
    # get values
    uncertain = findall(dat.certainty .== 1)

    # get number to remove
    nremove = Int64(round(size(uncertain)[1] * proportion))

    # kick 'em out
    removals = sort(sample(uncertain, nremove; replace = false))

    new_dat = deleteat!(deepcopy(dat), removals)

    return nremove, new_dat
end


# Create 100 permutations of each proportion of uncertain removals
# Define proportions of links to remove
props = [0.1, 0.25, 0.5]

# List of datasets to work with
datas = [fezouata_df_clean, burgess_df_clean, chengjiang_df_clean]

# Empty lists to store # of link removals & datasets with links removed
nums = []
outlist_uncertain = []

# Run it
for i in 1:length(datas)
    for j in 1:length(props)
        for k in 1:100
            # Get dataset with random uncertains removed
            temps = remove_uncertain(datas[i], props[j])
            
            if k == 1
                # Records #s removed for next step
                append!(nums,temps[1])
            end

            #Extract the new datsets and put in a list
            temps_frame = temps[2]
            outlist_uncertain = vcat(outlist_uncertain, temps_frame)
        end
    end
end


# Create function for random removals
function remove_random(dat, nremove)
    # kick 'em out
    removals = sort(sample(1:size(dat,1), nremove; replace = false))

    new_dat = deleteat!(deepcopy(dat), removals)

    return new_dat
end


# Create datasets with random removals
datas_rep = [datas[div(i,3)+1] for i=0:3*length(datas)-1]
outlist_random = []

for i in 1:length(nums)
    for k in 1:100
        # Get dataset with random uncertains removed
        temps = remove_random(datas_rep[i], nums[i])

        # Append to list
        outlist_random = vcat(outlist_random, temps)
    end
end


# Function to create networks
function make_network(df::DataFrame)
    # remove certainty values
    df = df[:,1:2]

    # make list of all unique species
    # note: there are still inconsistencies in species names that need to be tackled
    sp = unique(vcat(df.pred, df.prey))

    # count number of species 
    S = length(sp)

    # make adjacency matrix
    mat = zeros(Bool, S, S)

    for i in 1:S 
        for j in 1:S
            mat[i, j] = sum(df.pred .== sp[i] .&& df.prey .== sp[j])
        end
    end

    # change trophic species name for consistency 
    if eltype(sp) == Int64
        sp = string.(sp)
        sp = "s" .* sp
    end 

    # create network with species names 
    N = simplify(UnipartiteNetwork(mat, sp))
    
    return(N)
end

#=
# Creating the actual networks: this can take a few hours
uncertain_networks = Vector(undef, length(outlist_uncertain))
for i in 1:length(outlist_uncertain)
    uncertain_networks[i] = make_network(outlist_uncertain[i])
end

random_networks = Vector(undef, length(outlist_random))
for i in 1:length(outlist_uncertain)
    random_networks[i] = make_network(outlist_random[i])
end
=#

# create networks for full datasets

fez_net = make_network(fezouata_df_clean)
burg_net = make_network(burgess_df_clean)
cheng_net = make_network(chengjiang_df_clean)



# Network metrics function
function metrics(network, new_df)
   # simplify networks by removing isolated species
   network = simplify(network) 

   # calculate the number of species and links
   S = richness(network)
   L = links(network)
  
   # calculate the proportion of species that are top (without consumers), intermediate (with both consumers and resources), 
   # and basal (without resources)
   kin = values(degree(network, dims = 2))
   Top = sum(x -> x == 0, kin) ./ S
    
   # Basal (out-degree of 0)
   kout = values(degree(network, dims = 1))
   Bas = sum(x -> x == 0, kout) / S
    
   # Int (proportion of species that are not Top or Basal)
   Int = 1 - Top - Bas

   # calculate the proportion of species that are cannibals, herbivores (feeding only on basal species), 
   # omnivores (consuming two or more species with different trophic levels), 
   # and found in loops (food chains that contain the same species twice, apart from cannibalism)
    
   # Cannibals (proportion of species interacting with itself)
   Can = sum(diag(network.edges)) / S
    
   # Herbivores (proportion of species with a trophic level of 2)
   Herb = sum((values(trophic_level(network)) .== 2)) / S
    
   # Omnivores (proportion of species that consume two or more species and have food chains of different lengths)
   Omn = sum(values(omnivory(network)) .> 0) / S
   
   # Loops (proportion of species found in loops)

   # remove self-loops
   network.edges[diagind(network.edges)] .= 0
   # proportion of species with a path to itself (without self-loops)
   Loop = sum(diag(Matrix(shortest_path(network))) .> 0) / S
   

   # calculate the average length of food chains, the standard deviation of their length, and the log number of food chains
   food_chain_lengths = food_chains(network)
   
   # Average length of food chains
   ChLen = mean(food_chain_lengths)
     
   # Standard deviation of food-chain lengths
   ChSD = std(food_chain_lengths)
     
   # Log number of food chains
   ChNum = log10(length(food_chain_lengths))
     
   # calculate the mean trophic level of all species 
   TL = mean(values(trophic_level(network)))
     
   # calculate the average of the maximum trophic similarity of each species
   MxSim = MaxSim(network)
     
   # calculate the normalized standard deviations of vulnerability (nb of consumers or in-degree), generality (nb of resources or out-degree), 
   #and total links (nb of consumers and resources or total degree)
   # species in, out, and total degrees are normalized by the average number of interactions per species (2L/S)
     
   # Vulnerability
   VulSD = vulnerability(network)
     
   # Generality
   GenSD = generality(network)
     
   # Total links
   LinkSD = total_links(network)
     
   # calculate the average shortest food-chain length between all pairs of species
     
   # Average shortest path (not taking into account unconnected pairs)
   paths = shortest_path(network)
   Path = mean(paths[Not(paths .== 0)])

   # calculate the mean clustering coefficient, the probability that two species linked to the same species are also linked 
   # Mean clustering coefficient
   Clust = clustering_coefficient(network)

   row = [S, L, Top, Bas, Int, Can, Herb, Omn, Loop, ChLen, ChSD, ChNum, TL, MxSim, VulSD, GenSD, LinkSD, Path, Clust]
   push!(new_df, row)

end


# Calculate network metrics
# Write out names of data frame columns
entries = ["S", "L", "Top", "Bas", "Int", "Can", "Herb", "Omn", "Loop", "ChLen", "ChSD", "ChNum", "TL", "MxSim", "VulSD", "GenSD", 
               "LinkSD", "Path", "Clust"]

# Create empty data frames for storing outputs
empty_df_uncertain = DataFrame([name =>[] for name in entries])
empty_df_random = DataFrame([name =>[] for name in entries])

uncertain_df = Vector(undef, length(uncertain_networks))
random_df = Vector(undef, length(random_networks))

#Calculate network metrics (WARNING: this can take up to an hour, depending on computing power)
for i in 1:length(uncertain_networks)
   uncertain_df = metrics(uncertain_networks[i], empty_df_uncertain)
end

for i in 1:length(random_networks)
   random_df = metrics(random_networks[i], empty_df_random)
end


# Network metrics for full datasets
entries = ["S", "L", "Top", "Bas", "Int", "Can", "Herb", "Omn", "Loop", "ChLen", "ChSD", "ChNum", "TL", "MxSim", "VulSD", "GenSD", "LinkSD", "Path", "Clust"]
empty_df_real = DataFrame([name =>[] for name in entries])

real_datas = [fez_net, burg_net, cheng_net] # Change this line depending on the data set

for i in 1:length(real_datas)
    empty_df_real = metrics(real_datas[i], empty_df_real)
end


cd("results/uncertainty_analysis/dunne_method") do
    CSV.write("real.csv", empty_df_real)
end


# Save network metrics (figures will be made in R)
cd("results/uncertainty_analysis/dunne_method") do
    CSV.write("uncertain.csv", uncertain_df)
    CSV.write("random.csv", random_df)
end



### Uncertainty analysis: removals by trophic level

# Previous results are very similar to the Dunne method. The analysis takes a long time to run and there may be issues with memory limits. 

# Function assigning trophic roles to species
function assign_roles(network)
    # simplify networks by removing isolated species
    network = simplify(network) 

    # Fezouata has very few top species, so leave out for now.
    #get "top" species (without consumers)
    kin = degree(network, dims=2)
    tops = filter(((k,v),) -> v == 0, kin)

    # get data frame started
    trophic_roles = DataFrame(Species = collect(keys(tops)), Role = "Top")

    # get basal species (without resource species) and add to df
    kout = degree(network, dims = 1)
    bas = filter(((k,v),) -> v == 0, kout)

    trophic_roles = DataFrame(Species = collect(keys(bas)), Role = "Basal")

    # Get herbivores & add to df
    other_Roles = trophic_level(network)
    Herb = filter(((k,v),) -> v == 2, other_Roles)

    append!(trophic_roles, DataFrame(Species = collect(keys(Herb)), Role = "Herb"))

    # Get omnivores & add to df
    omn_vals = omnivory(network)
    Omn = filter(((k,v),) -> v > 0, omn_vals)

    append!(trophic_roles, DataFrame(Species = collect(keys(Omn)), Role = "Omn"))
end

fez_roles = assign_roles(fez_net)
burg_roles = assign_roles(burg_net)
cheng_roles = assign_roles(cheng_net)


# Remove links based on trophic role and uncertainty level
function remove_uncertain_roles(role_df, role, dat, proportion)
    # Prep raw data frame, if necessary
    if eltype(dat.pred) == Int64
        dat.pred = "s" .* string.(dat.pred)
        dat.prey = "s" .* string.(dat.prey)
    end

    # Remove linkages based on uncertainty and trophic role 
    role_sub = filter(:Role =>x -> x == role, role_df)

    uncertain = findall(dat.certainty .== 1)
    trophic = unique(vcat(findall(in(role_sub.Species), dat.pred), findall(in(role_sub.Species), dat.prey)))

    # get values
    uncertain_trophic = intersect(uncertain, trophic)

    # get number to remove
    nremove = Int64(round(length(uncertain_trophic) * proportion))

    # kick 'em out
    removals = sort(sample(uncertain_trophic, nremove; replace = false))

    new_dat = deleteat!(deepcopy(dat), removals)

    return nremove, new_dat
end


# Create 100 more permutations of each proportion of uncertain removals
props = [0.1, 0.25, 0.5]
roles = ["Basal", "Herb", "Omn"]
datas = [fezouata_df_clean, burgess_df_clean, chengjiang_df_clean]
role_frames = [fez_roles, burg_roles, cheng_roles]

nums = []
role = []
outlist_uncertain = []

for i in 1:length(datas)
    for j in 1:length(props)
        for l in 1:length(roles)
            for k in 1:100
                # Get dataset with random uncertains removed
                temps = remove_uncertain_roles(role_frames[i], roles[l], datas[i], props[j])
            
                if k == 1
                    # Records #s removed for next step
                    append!(nums,temps[1])
                    role = vcat(role, roles[l])
                end

                #Extract the new datsets and put in a list
                temps_frame = temps[2]
                outlist_uncertain = vcat(outlist_uncertain, temps_frame)
            end
        end
    end
end

# Create function for random removals at a given trophic guild
function remove_random_uncertain(dat, role_df, role, nremove)
    # Prep raw data frame, if necessary
    if eltype(dat.pred) == Int64
        dat.pred = "s" .* string.(dat.pred)
        dat.prey = "s" .* string.(dat.prey)
    end
   
    # Remove linkages based on trophic role 
    role_sub = filter(:Role =>x -> x == role, role_df)
    
    # get indices to potentially remove
    indices = unique(vcat(findall(in(role_sub.Species), dat.pred), findall(in(role_sub.Species), dat.prey)))

    # kick 'em out
    removals = sort(sample(indices, nremove; replace = false))

    new_dat = deleteat!(deepcopy(dat), removals)

    return new_dat
end

# Create datasets with random removals
datas_rep = [datas[div(i,9)+1] for i=0:9*length(datas)-1]
roles_df_rep = [role_frames[div(i,9)+1] for i=0:9*length(role_frames)-1]
outlist_random = []

for i in 1:length(nums)
    for k in 1:100
        # Get dataset with random uncertains removed
        temps = remove_random_uncertain(datas_rep[i], roles_df_rep[i], role[i], nums[i])

        # Append to list
        outlist_random = vcat(outlist_random, temps)
    end
end

# Create networks
uncertain_networks = Vector(undef, length(outlist_uncertain))

# This loop: approx. 2 hours
@showprogress @distributed for i in 1:length(outlist_uncertain)
    uncertain_networks[i] = make_network(outlist_uncertain[i])
    sleep(1)
end

# save_object("results/uncertainty_analysis/removal_by_trophic_levels/uncertain_fez_intermediate.jld2", uncertain_networks)
#=
# Takes another 2 hours
random_networks = Vector(undef, length(outlist_random))
@showprogress @distributed for i in 1:length(outlist_uncertain)
    random_networks[i] = make_network(outlist_random[i])
    sleep(1)
end

save_object("results/uncertainty_analysis/removal_by_trophic_levels/random_fez_intermediate.jld2", random_networks)
=#

# Calculate network metrics
entries = ["S", "L", "Top", "Bas", "Int", "Can", "Herb", "Loop", "Omn", "ChLen", "ChSD", "ChNum", "TL", "MxSim", "VulSD", "GenSD", "LinkSD", "Path", "Clust"]
empty_df_uncertain = DataFrame([name =>[] for name in entries])
#empty_df_random = DataFrame([name =>[] for name in entries])

# Load intermediate products if VSCode decides to crap out:
# uncertain_networks = load_object("results/uncertainty_analysis/removal_by_trophic_levels/uncertain_fez_intermediate.jld2")
#random_networks = load_object("results/uncertainty_analysis/removal_by_trophic_levels/random_fez_intermediate.jld2")

uncertain_df = Vector(undef, length(uncertain_networks))
#random_df = Vector(undef, length(random_networks))

@showprogress for i in 1:length(uncertain_networks)
   uncertain_df = metrics(uncertain_networks[i], empty_df_uncertain)
   sleep(1)
end

# Save it
cd("results/uncertainty_analysis/removal_by_trophic_levels") do
   CSV.write("uncertain_fez_roles.csv", uncertain_df)
end

#=
# do the same for random networks

@showprogress for i in 1:length(random_networks)
   random_df = metrics(random_networks[i], empty_df_random)
   sleep(1)
end

cd("results/uncertainty_analysis/removal_by_trophic_levels") do
   CSV.write("random_fez_roles.csv", random_df)
end
=#
