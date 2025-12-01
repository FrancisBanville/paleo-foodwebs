#### Clean datasets of ecological interactions 

### read raw datasets
fezouata_df = DataFrame(CSV.File(joinpath("data", "raw", "interactions_Fezouata.csv")))

### clean datasets 

function clean_data(df::DataFrame)
    
    # rename variables 
    rename!(df, 1 => :pred, 2 => :prey)

    # remove empty rows since they do not represent interactions
    dropmissing!(df)

    if eltype(df[:,1]) !== Int64
        
        # convert uppercase letters to lowercase
        df = lowercase.(df)

        # remove symbols that artificially creates new species when inconsistent 
        df = replace.(df, "?" => "")
        df = replace.(df, "'" => "")
        df = replace.(df, "\"" => "")

        # remove leading and trailing white spaces
        df = strip.(df) 
    end
    
    # remove duplicate rows
    df = unique(df)

    return(df)
end

fezouata_df_clean = clean_data(fezouata_df)


### convert datasets to networks

function make_network(df::DataFrame)

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

# create and save networks with species names 

N_fezouata = make_network(fezouata_df_clean)


save(joinpath("data", "clean", "network_Fezouata.jld2"), "N", N_fezouata)


