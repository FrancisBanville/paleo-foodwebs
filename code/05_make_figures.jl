## Make figures showing the distribution of measures

# plot attributes
theme(:mute)
default(; frame=:box)
Plots.scalefontsizes(1.3)
fonts=font("Times",8)

# read dataset of network measures
measures = DataFrame(CSV.File(joinpath("results", "measures.csv")))

# remove NaNs and missing values
measures[isnan.(measures.ChLen),:ChLen] .= 0
measures[isnan.(measures.ChSD),:ChSD] .= 0

# remove infinite values that were obtained after taking the log of 0 
measures[measures.ChNum .== -Inf, :ChNum] .= 0

# function for density plots 

function plot_density(network::String,
                    measure::String, 
                    xlab::String, 
                    xlim::Tuple)
    
    # subset dataset according to the type of network (empirical and predicted using the three models)
    measures_network = filter(x -> x.network == network, measures)
    
    measures_empirical = filter(x -> x.type == "empirical", measures_network)
    measures_niche = filter(x -> x.type == "niche model", measures_network)
    measures_cascade = filter(x -> x.type == "cascade model", measures_network)
    measures_nested_hierarchy = filter(x -> x.type == "nested hierarchy model", measures_network)
    
    # density plot of niche model measures
    density(measures_niche[:, measure], 
                label="",
                fill=(0, .5),
                linewidth=2, 
                framestyle=:box, 
                grid=false,
                minorgrid=false,
                dpi=1000, 
                size=(800,500), 
                margin=5Plots.mm, 
                guidefont=fonts, 
                xtickfont=fonts, 
                ytickfont=fonts,
                foreground_color_legend=nothing, 
                background_color_legend=:white, 
                legendfont=fonts,
                legendfontpointsize=8,
                legendfontfamily="Times")

    # density plot of cascade model measures
    density!(measures_cascade[:, measure], 
            label="",
            fill=(0, .5),
            linewidth=2)

    # density plot of nested hierarchy model measures
    density!(measures_nested_hierarchy[:, measure], 
            label="",
            fill=(0, .5),
            linewidth=2)

    # vertical line for empirical measure
    plot!(measures_empirical[:, measure], 
        seriestype = :vline, 
        label="",
        linewidth=2,
        linestyle=:dash,
        color=:red)
    
    # x and y labs
    xaxis!(xlabel=xlab, 
        xlims=xlim)
    yaxis!(ylabel="Density", 
        ylims=(0, Inf))
end

# make plots

function plot_density_all(network::String)
        
    # filter datasets 
    measures_network = filter(x -> x.network == network, measures)
    
    measures_empirical = filter(x -> x.type == "empirical", measures_network)
    measures_niche = filter(x -> x.type == "niche model", measures_network)
    measures_cascade = filter(x -> x.type == "cascade model", measures_network)
    measures_nested_hierarchy = filter(x -> x.type == "nested hierarchy model", measures_network)
    
    # plots
    plot1 = plot_density(network,"Top", "Proportion of top predators", (-0.005,1))
    plot2 = plot_density(network, "Bas", "Proportion of basal species", (-0.005,1))
    plot3 = plot_density(network, "Int", "Proportion of intermediate species", (-0.005,1))
    plot4 = plot_density(network, "Can", "Proportion of cannibal species", (-0.005,1))
    plot5 = plot_density(network, "Herb", "Proportion of herbivore species", (-0.005,1))
    plot6 = plot_density(network, "Omn", "Proportion of omnivore species", (-0.005,1))
    plot7 = plot_density(network, "Loop", "Proportion of species in loops", (-0.005,1))
    plot8 = plot_density(network, "ChLen", "Average food chain length", (0,5))
    plot9 = plot_density(network, "ChSD", "Standard deviation of food chain length", (0,3))
    plot10 = plot_density(network, "ChNum", "Log number of food chains", (0,3.5))
    plot11 = plot_density(network, "TL", "Average trophic level", (0,7))
    plot12 = plot_density(network, "MxSim", "Average maximum similarity between species pairs", (0,1))
    plot13 = plot_density(network, "VulSD", "Standard deviation of vulerability", (0,1.6))
    plot14 = plot_density(network, "GenSD", "Standard deviation of generality", (0,1.6))
    plot15 = plot_density(network, "LinkSD", "Standard deviation of total links", (0,1.6))
    plot16 = plot_density(network, "Path", "Average shortest path length between species pairs", (0,7))
    plot17 = plot_density(network, "Clust", "Average clustering coefficient", (0,0.6))

    # legend plot
    plot_legend = density(measures_niche[:, "Top"], 
        label="Niche model",
        legend=:top,
        fill=(0, .5),
        linewidth=2, 
        xlims=(1000, 1001),
        ylims=(1000, 1001),
        grid=false,
        minorgrid=false,
        framestyle=:none,
        dpi=1000, 
        size=(800,500), 
        margin=5Plots.mm, 
        guidefont=fonts, 
        foreground_color_legend=nothing, 
        background_color_legend=:white, 
        legendfont=fonts,
        legendfontpointsize=8,
        legendfontfamily="Times")

    density!(measures_cascade[:, "Top"], 
        label="Cascade model",
        fill=(0, .5),
        linewidth=2)

    density!(measures_nested_hierarchy[:, "Top"], 
        label="Nested hierarchy model",
        fill=(0, .5),
        linewidth=2)

    plot!(measures_empirical[:, "Top"], 
        seriestype = :vline, 
        label="Empirical network",
        linewidth=2,
        linestyle=:dash,
        color=:red)

    # group plots

    # size factor
    x = 2

    plot(plot1, plot2, plot3, plot4, 
        plot5, plot6, plot7, plot8, 
        plot9, plot10, plot11, plot12, 
        plot13, plot14, plot15, plot16,
        plot17, plot_legend,
        title = ["Top" "Bas" "Int" "Can" "Herb"  "Omn" "Loop" "ChLen" "ChSD" "ChNum" "TL" "MxSim" "VulSD" "GenSD" "LinkSD" "Path" "Clust" ""],
        titleloc=:right, 
        titlefont=fonts,
        layout = (4, 5), 
        size=(800*x, 500*x))
end

plot_density_all("fezouata")
savefig(joinpath("figures", "eco_models", "fezouata_density_measures.png"))

plot_density_all("fezouata trophic species")
savefig(joinpath("figures", "eco_models", "fezouata_trophicsp_density_measures.png"))


plot_density_all("burgess")
savefig(joinpath("figures", "eco_models", "burgess_density_measures.png"))

plot_density_all("burgess trophic species")
savefig(joinpath("figures", "eco_models", "burgess_trophicsp_density_measures.png"))


plot_density_all("chengjiang")
savefig(joinpath("figures", "eco_models", "chengjiang_density_measures.png"))

plot_density_all("chengjiang trophic species")
savefig(joinpath("figures", "eco_models", "chengjiang_trophicsp_density_measures.png"))




## Make figures showing the distribution of predictive errors of niche model

# create empty dataset for predictive errors
measures_errors = DataFrame()
    
# calculate predictive errors of niche model for all networks and measures 
for N in unique(measures.network) 
    
    # get empirical measures of network N
    measures_emp = measures[measures.network .== N .&& measures.type .== "empirical", 3:end]

    # get niche model predictions for network N
    measures_niche = measures[measures.network .== N .&& measures.type .== "niche model", 3:end]

    # calculate predictive errors for network N
    measures_errors_N = (measures_niche .- measures_emp) ./ measures_emp
    
    # add network name column 
    measures_errors_N.network = fill(N, size(measures_errors_N, 1))

    # add to general dataframe
    measures_errors = [measures_errors; measures_errors_N]
end

# remove NaNs and missing values
measures_errors[isnan.(measures_errors.Can),:Can] .= Inf
measures_errors[isnan.(measures_errors.Loop),:Loop] .= Inf

# colors 
pal = [RGB(230/255,159/255,0/255),
    RGB(86/255,190/255,233/255),
    RGB(204/255,121/255,167/255)]

# function for a single measure and a set of networks 
function plot_density_errors(networks::Vector,
                                measure::String)

        # empty density plot of predictive errors 
        density(
            framestyle=:box, 
            grid=false,
            minorgrid=false,
            dpi=1000, 
            size=(800,500), 
            margin=5Plots.mm, 
            guidefont=fonts, 
            xtickfont=fonts, 
            ytickfont=fonts,
            foreground_color_legend=nothing, 
            background_color_legend=:white, 
            legendfont=fonts,
            legendfontpointsize=8,
            legendfontfamily="Times")
        
    # density plots of network errors for the all networks
    for N in 1:length(networks)
        
        # do not plot if relative error is infinite (empirical measure = 0)
        if !isinf(measures_errors[measures_errors.network .== networks[N], measure][1])
        
            density!(measures_errors[measures_errors.network .== networks[N], measure], 
            color=pal[N],
            label="",
            fill=(0, .5),
            linewidth=2)
        end
    end

    # vertical line at x = 0
    plot!([0], 
        seriestype = :vline, 
        label="",
        linewidth=1,
        linestyle=:dash,
        color=:grey)

    # x and y labs
    xaxis!(xlabel="Relative predictive error")
    
    yaxis!(ylabel="Density", 
            ylims=(0, Inf))
end


# plot all measures

function plot_density_errors_all(networks::Vector)
        
    plot1 = plot_density_errors(networks,"Top")
    plot2 = plot_density_errors(networks, "Bas")
    plot3 = plot_density_errors(networks, "Int")
    plot4 = plot_density_errors(networks, "Can")
    plot5 = plot_density_errors(networks, "Herb")
    plot6 = plot_density_errors(networks, "Omn")
    plot7 = plot_density_errors(networks, "Loop")
    plot8 = plot_density_errors(networks, "ChLen")
    plot9 = plot_density_errors(networks, "ChSD")
    plot10 = plot_density_errors(networks, "ChNum")
    plot11 = plot_density_errors(networks, "TL")
    plot12 = plot_density_errors(networks, "MxSim")
    plot13 = plot_density_errors(networks, "VulSD")
    plot14 = plot_density_errors(networks, "GenSD")
    plot15 = plot_density_errors(networks, "LinkSD")
    plot16 = plot_density_errors(networks, "Path")
    plot17 = plot_density_errors(networks, "Clust")

    # legend plot
    plot_legend = density(
        xlims=(1000, 1001),
        ylims=(1000, 1001),
        grid=false,
        minorgrid=false,
        framestyle=:none,
        dpi=1000, 
        size=(800,500), 
        margin=5Plots.mm, 
        guidefont=fonts, 
        foreground_color_legend=nothing, 
        background_color_legend=:white, 
        legendfont=fonts,
        legendfontpointsize=8,
        legendfontfamily="Times")

    for N in 1:length(networks)
        density!(measures_errors[:, "Top"],
            color=pal[N],
            label=networks[N],
            fill=(0, .5),
            linewidth=2)
    end

    # size factor
    x = 2

    plot(plot1, plot2, plot3, plot4, 
        plot5, plot6, plot7, plot8, 
        plot9, plot10, plot11, plot12, 
        plot13, plot14, plot15, plot16,
        plot17, plot_legend,
        title = ["Top" "Bas" "Int" "Can" "Herb"  "Omn" "Loop" "ChLen" "ChSD" "ChNum" "TL" "MxSim" "VulSD" "GenSD" "LinkSD" "Path" "Clust" ""],
        titleloc=:right, 
        titlefont=fonts,
        layout = (4, 5), 
        size=(800*x, 500*x))
end


# density plots of model errors for species networks
sp_networks = ["fezouata", "burgess", "chengjiang"]

plot_density_errors_all(sp_networks)
savefig(joinpath("figures", "niche_errors", "species_networks_errors.png"))

# density plots of model errors for trophic species networks
trophicsp_networks = ["fezouata trophic species", "burgess trophic species", "chengjiang trophic species"]

plot_density_errors_all(trophicsp_networks)
savefig(joinpath("figures", "niche_errors", "trophic_species_networks_errors.png"))
