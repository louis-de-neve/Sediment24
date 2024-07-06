using Oceananigans, OceanBioME
using CairoMakie
using Printf
using Oceananigans.Units
using FileIO

include("SimpleColumnPLot.jl")

function PAR_daily_fluctuation(t)
    return (cos(t * π / 12hours + π) + 1)
end
function PAR_Aqua_MODIS_fluctuations(t)
    return (18.9 * sin((2 * π * t / 365days) + 1.711) + 46.62)
end

@inline function ImportPreviousState(model, initial_filename::String="SimpleColumnSave")::Nothing
    @info "Importing previous data"
    state_dict = Dict()
    for tracer in keys(model.tracers)
        data_fts = FieldTimeSeries("$initial_filename.jld2", string(tracer))
        tracer_previous_state = parent(data_fts[:, :, :, end[]])
        state_dict[tracer] = tracer_previous_state
    end
    set!(model; state_dict...)
end


@inline function DefineColumnModel(sim_height::Int64=50,
                                   initial_file_name::String="SimpleColumnSave",
                                   continue_from_previous::Bool=False)::Tuple{Oceananigans.AbstractModel, Int64}
    @info "Defining model..."

    if continue_from_previous
        initial_time = load("LastSimulationTimestep.jld2", "time")
    else
        initial_time = 0
    end

    function Aqua_MODIS_PAR(t; initial_time=0)
        return PAR_Aqua_MODIS_fluctuations(t+initial_time) * PAR_daily_fluctuation(t+initial_time)
    end

    grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (sim_height, ), extent = (sim_height, ))

    bgc = LOBSTER(; grid,
                    carbonates=true,
                    oxygen=true,
                    variable_redfield=true,
                    surface_photosynthetically_active_radiation = Aqua_MODIS_PAR)

    model = NonhydrostaticModel(; grid,
                                #timestepper= :RungeKutta3,
                                biogeochemistry=bgc, 
                                closure = ScalarDiffusivity(ν = κₜ, κ = κₜ)
                                )
    #=                          
    set!(model, P = 0.4686, Z = 0.5363,     
                NO₃ = 2.3103, NH₄ = 0.0010, 
                DIC = 2106.9, Alk = 2408.9, 
                O₂ = 258.92, 
                DOC = 5.3390, DON = 0.8115,
                sPON = 0.2299, sPOC = 1.5080,
                bPON = 0.0103, bPOC = 0.0781) =#
    P(z) = 0.1 .- z./250
    Z(z) = -z/500

    set!(model, P = P, Z = 0.2,     
                NO₃ = 2.3103, NH₄ = 0.0010, 
                DIC = 2106.9, Alk = 2408.9, 
                O₂ = 258.92, 
                DOC = 5.3390, DON = 0.8115,
                sPON = 0.2299, sPOC = 1.5080,
                bPON = 0.0103, bPOC = 0.0781)
    if continue_from_previous
        ImportPreviousState(model, initial_file_name)
    end

    return model, initial_time
end


@inline function CreateColumnSimulation(model::Oceananigans.AbstractModel,
                                        outfile_location::String="ColumnOutput",
                                        sim_time::Float64=10days,
                                        sim_timestep::Int64=100
                                        )::Oceananigans.Simulation
    @info "Setting up simulation..."

    simulation = Simulation(model, Δt = sim_timestep, stop_time = sim_time)
    #println(model.tracers)
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

SIMULATION_HEIGHT = 50 # meters
SIMULATION_TIME = 365days # seconds or any Oceananigans unit
SIMULATION_TIMESTEP = 100 # seconds
START_SAVEFILE_NAME = "SimpleColumnSave"
FINAL_SAVEFILE_NAME = "SimpleColumnSave1"

model, initial_time = DefineColumnModel(SIMULATION_HEIGHT, START_SAVEFILE_NAME, true)
simulation = CreateColumnSimulation(model, FINAL_SAVEFILE_NAME, SIMULATION_TIME, SIMULATION_TIMESTEP)
run!(simulation)
final_time = initial_time + SIMULATION_TIME
save("LastSimulationTimestep.jld2", "time", final_time)
MakePlotOfColumn(FINAL_SAVEFILE_NAME, initial_time, ("P", "Z"), ("Phytoplankton", "Zooplankton"))
