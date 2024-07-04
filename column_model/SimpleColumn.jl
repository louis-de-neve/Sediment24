using Oceananigans, OceanBioME
using CairoMakie
using Printf
using Oceananigans.Units

κₜ = 0 #1e-6

@info "Defining model..."

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (20, ), extent = (40, ))

@inline function PAR_daily_fluctuation(t)
    return max(0, (cos(t * π / 12hours)))
end
@inline function PAR_yearly_fluctuation(t)
    return (cos(t * 2 * π / 365days) + 1)
end
@inline function default_surface_PAR(t)
    return 100 * PAR_daily_fluctuation(t) * PAR_yearly_fluctuation(t)
end

bgc = LOBSTER(; grid,
                carbonates=true,
                oxygen=true,
                variable_redfield=true,
                surface_photosynthetically_active_radiation = default_surface_PAR)

model = NonhydrostaticModel(; grid,
                              biogeochemistry=bgc, 
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ)
                              )

set!(model, P = 0.4686, Z = 0.5363,     
            NO₃ = 2.3103, NH₄ = 0.0010, 
            DIC = 2106.9, Alk = 2408.9, 
            O₂ = 258.92, 
            DOC = 5.3390, DON = 0.8115,
            sPON = 0.2299, sPOC = 1.5080,
            bPON = 0.0103, bPOC = 0.0781)

@info "Setting up simulation..."

simulation = Simulation(model, Δt = 100, stop_time = 10days)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                               filename = "SimpleColumn.jld2",
                                               schedule = TimeInterval(24minute),
                                               overwrite_existing = true)

start_time = time_ns() # record the start time

# Build a progress message
progress(sim) = @printf("i: % 6d, sim time: % 4s, wall time: % 10s \n",
    sim.model.clock.iteration,
    prettytime(sim.model.clock.time),
    prettytime(1e-9 * (time_ns() - start_time)))

# Display the progress message every 1000 timesteps
simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))

@info "Simulation setup, running..."
run!(simulation)
