using Oceananigans, OceanBioME
using CairoMakie
using Printf
using Oceananigans.Units

include("SimpleColumnPLot.jl")

function PAR_daily_fluctuation(t)
    return max(0, (cos(t * π / 12hours)))
end
function PAR_yearly_fluctuation(t)
    return (cos(t * 2 * π / 365days) + 1)
end
function default_surface_PAR(t)
    return 100 * PAR_daily_fluctuation(t) * PAR_yearly_fluctuation(t)
end


@inline function DefineColumnModel(sim_height::Int64=50)::Oceananigans.AbstractModel
    @info "Defining model..."

    grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (sim_height, ), extent = (sim_height, ))

    bgc = LOBSTER(; grid,
                    carbonates=true,
                    oxygen=true,
                    variable_redfield=true,
                    surface_photosynthetically_active_radiation = default_surface_PAR)

    model = NonhydrostaticModel(; grid,
                                #timestepper= :RungeKutta3,
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

    return model
end


@inline function CreateColumnSimulation(model::Oceananigans.AbstractModel,
                                        outfile_location::String="ColumnOutput",
                                        sim_time::Float64=10days,
                                        sim_timestep::Int64=100
                                        )::Oceananigans.Simulation
    @info "Setting up simulation..."

    simulation = Simulation(model, Δt = sim_timestep, stop_time = sim_time)

    simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                filename = "$outfile_location.jld2",
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

    return simulation
end


κₜ = 1e-6 #1e-6 (diffustion constant)

const SIMULATION_HEIGHT = 50 # meters
const SIMULATION_TIME = 150days # seconds or any Oceananigans unit
const SIMULATION_TIMESTEP = 100 # seconds
const SAVEFILE_NAME = "SimpleColumnSave"

model = DefineColumnModel(SIMULATION_HEIGHT)
simulation = CreateColumnSimulation(model, SAVEFILE_NAME, SIMULATION_TIME, SIMULATION_TIMESTEP)
run!(simulation)
MakePlotOfColumn(savefile_name, ("P", "Z"), ("Phytoplankton", "Zooplankton"))
