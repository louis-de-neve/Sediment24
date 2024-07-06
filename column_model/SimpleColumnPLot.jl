using CairoMakie
using Printf
using Oceananigans.Units
using Oceananigans
include("TracerInfo.jl")

function make_plot_of_column_simulation!(start_time::Int64,
                                        tracer_infos_to_plot::Vector{TracerInfo}
                                        )::Nothing
    if typeof(tracer_infos_to_plot) == Nothing
        return nothing
    end
    
    @info "Plotting column..."

    number_of_tracers = length(tracer_infos_to_plot)
    
    fig = Figure(size = (1000, 300*number_of_tracers), fontsize=20)

    first_tracer_fts = FieldTimeSeries("temp_plotting_data.jld2", tracer_infos_to_plot[1].tracer_name_in_model)
    xc, yc, zc = nodes(first_tracer_fts)
    times = first_tracer_fts.times .+ start_time
    tick_location_seconds = range(start_time, times[end]; length=5)
    tick_location_days = tick_location_seconds / (24 * 60 * 60)

    axis_kwargs = (xlabel="Days since January 1st 2023", ylabel="z / m", width=800, height=150)

    for (index, tracer_info) in enumerate(tracer_infos_to_plot)
        tracer_name = tracer_info.tracer_name_in_model
        tracer_title = tracer_info.display_name
        tracer_units = tracer_info.unit

        data_fts = FieldTimeSeries("temp_plotting_data.jld2", tracer_name)      
        ax = Axis(fig[index, 1]; title=tracer_title, axis_kwargs...)
        ax.xticks = (collect(tick_location_seconds), string.(collect(tick_location_days)))
        
        hm = heatmap!(ax, times, zc, log10.(abs.(data_fts[1, 1, 1:end, 1:end]')), interpolate=true)
        Colorbar(fig[index, 2], hm, label="log₁₀ $tracer_units")
    end

    Label(fig[begin-1, 1:2], 
          "LOBSTER biogeochemical tracers with accurate\n photosynthetically active radiation flux in Luderitz, Namibia",
          fontsize=30,
          padding=(0, 0, 0, 0))
    display(fig)
    return nothing
end
