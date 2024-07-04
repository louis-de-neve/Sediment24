using CairoMakie
using Printf
using Oceananigans.Units
using Oceananigans

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (20, ), extent = (40, ))

P = FieldTimeSeries("SimpleColumn.jld2", "P")
Z = FieldTimeSeries("SimpleColumn.jld2", "Z")

xc, yc, zc = nodes(P)#grid, Center(), Center(), Center())

times = P.times

tick_location_seconds = range(0, times[end]; length=5)
tick_location_days = tick_location_seconds/(24*60*60)

axis_kwargs = (xlabel = "time / days", ylabel = "z / m", width = 800, height = 150)

fig = Figure(size = (1000, 600), fontsize = 20)
ax1 = Axis(fig[1,1]; title = "Phytoplankton concentration", axis_kwargs...)
ax2 = Axis(fig[2,1]; title = "Zooplankton concentration", axis_kwargs...)
ax1.xticks = (collect(tick_location_seconds), string.(collect(tick_location_days)))
ax2.xticks = (collect(tick_location_seconds), string.(collect(tick_location_days)))


hmP = heatmap!(ax1, times, zc, log10.(abs.(P[1, 1, 1:grid.Nz, 1:end]')))
hmZ = heatmap!(ax2, times, zc, log10.(abs.(Z[1, 1, 1:grid.Nz, 1:end]')))


Colorbar(fig[1, 2], hmP, label="log₁₀ mmol N / m³")
Colorbar(fig[2, 2], hmZ, label="log₁₀ mmol N / m³")

fig