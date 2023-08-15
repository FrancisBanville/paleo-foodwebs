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

# subset dataset according to the type of network (empirical and predicted using the three models)

measures_empirical = filter(x -> x.type == "empirical", measures)
measures_niche = filter(x -> x.type == "niche model", measures)
measures_cascade = filter(x -> x.type == "cascade model", measures)
measures_nested_hierarchy = filter(x -> x.type == "nested hierarchy model", measures)

# function for density plots 

function plot_density(measure::String, 
                    xlab::String, 
                    xlim::Tuple)
    
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
plot1 = plot_density("Top", "Proportion of top predators", (0,1))
plot2 = plot_density("Bas", "Proportion of basal species", (0,1))
plot3 = plot_density("Int", "Proportion of intermediate species", (0,1))
plot4 = plot_density("Can", "Proportion of cannibal species", (-0.005,1))
plot5 = plot_density("Herb", "Proportion of herbivore species", (0,1))
plot6 = plot_density("Omn", "Proportion of omnivore species", (0,1))
plot7 = plot_density("Loop", "Proportion of species in loops", (-0.005,1))
plot8 = plot_density("ChLen", "Average food chain length", (0,5))
plot9 = plot_density("ChSD", "Standard deviation of food chain length", (0,3))
plot10 = plot_density("ChNum", "Log number of food chains", (0,3.5))
plot11 = plot_density("TL", "Average trophic level", (0,6))
plot12 = plot_density("MxSim", "Average maximum similarity between species pairs", (0,1))
plot13 = plot_density("VulSD", "Standard deviation of vulerability", (0,1.6))
plot14 = plot_density("GenSD", "Standard deviation of generality", (0,1.6))
plot15 = plot_density("LinkSD", "Standard deviation of total links", (0,1.6))
plot16 = plot_density("Path", "Average shortest path length between species pairs", (0,7))
plot17 = plot_density("Clust", "Average clustering coefficient", (0,0.2))

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

savefig(joinpath("figures","density_measures.png"))

