using CairoMakie
using Printf
using Oceananigans.Units
using Oceananigans

function MakePlotOfColumn(filename="ColumnOutput",
                          plotting_tracers::Tuple=("P", "Z"),
                          tracer_titles::Tuple=("Phytoplankton", "Zooplankton")
                          )
    @info "Plotting column..."


    number_of_tracers = length(plotting_tracers)
    
    fig = Figure(size = (1000, 300*number_of_tracers), fontsize = 20)

    first_tracer_fts = FieldTimeSeries("$filename.jld2", plotting_tracers[1])
    xc, yc, zc = nodes(first_tracer_fts)
    times = first_tracer_fts.times
    tick_location_seconds = range(0, times[end]; length=5)
    tick_location_days = tick_location_seconds/(24*60*60)

    axis_kwargs = (xlabel = "time / days", ylabel = "z / m", width = 800, height = 150)

    for (index, tracer) in enumerate(plotting_tracers)

        data_fts = FieldTimeSeries("$filename.jld2", tracer)
        
        tracer_title = tracer_titles[i]
        ax = Axis(fig[index, 1]; title = "$tracer_title concentration", axis_kwargs...)
        ax.xticks = (collect(tick_location_seconds), string.(collect(tick_location_days)))
        
        hm = heatmap!(ax, times, zc, log10.(abs.(data_fts[1, 1, 1:end, 1:end]')), interpolate=true)
        Colorbar(fig[index, 2], hm, label="log₁₀ mmol N / m³")

    end
    return fig
end

