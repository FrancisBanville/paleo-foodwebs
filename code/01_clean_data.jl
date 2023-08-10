#### Clean dataset of ecological interactions 

# read raw dataset 
int_df = DataFrame(CSV.File(joinpath("data", "raw", "interactions_Fezouata.csv")))

# rename variables 
rename!(int_df, 1 => :pred, 2 => :prey)

# remove empty rows since they do not represent interactions
dropmissing!(int_df)

# convert uppercase letters to lowercase
int_df = lowercase.(int_df)

# remove symbols that artificially creates new species when inconsistent 
int_df = replace.(int_df, "?" => "")
int_df = replace.(int_df, "'" => "")
int_df = replace.(int_df, "\"" => "")

# remove leading and trailing white spaces
int_df = strip.(int_df) 


#### Make network

# make list of all unique species
# note: there are still inconsistencies in species names that need to be tackled
sp = unique(vcat(int_df.pred, int_df.prey))

# count number of species 
S = length(sp)

# make adjacency matrix
mat = zeros(Bool, S, S)

for i in 1:S 
    for j in 1:S
        mat[i, j] = sum(int_df.pred .== sp[i] .&& int_df.prey .== sp[j])
    end
end

# create and save network with species names 
N = UnipartiteNetwork(mat, sp)

save(joinpath("data", "clean", "network_Fezouata.jld2"), "N", N)
