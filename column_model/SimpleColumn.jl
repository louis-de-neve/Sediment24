using Oceananigans, OceanBioME
using CairoMakie
using Printf
using Oceananigans.Units
using FileIO

include("SimpleColumnPLot.jl")

@inline function import_previous_model_state!(model)::Nothing
    @info "Importing previous data"

    state_dict = Dict()

    for tracer in keys(model.tracers)
        data_fts = FieldTimeSeries("temp_final_state.jld2", string(tracer))
        tracer_previous_state = parent(data_fts[:, :, :, end[]])
        state_dict[tracer] = tracer_previous_state
    end
    set!(model; state_dict...)

    return nothing
end

par_daily_fluctuation(t) = cos(t * π / 12hours + π) + 1

par_aqua_modis_fluctuations(t) = 18.9 * sin((2 * π * t / 365days) + 1.711) + 46.62

@inline function define_par_function(t; initial_time=0)
    rebased_time = t + initial_time
    return par_aqua_modis_fluctuations(rebased_time) * par_daily_fluctuation(rebased_time)
end

@inline function setup_column_model(sim_height::Int64=50,
                                    κₜ = 1e-6,
                                    continue_from_previous::Bool=false
                                    )::Tuple{Oceananigans.AbstractModel, Int64}
    @info "Defining model..."

    if continue_from_previous
        initial_time = load("LastSimulationTimestep.jld2", "time")
    else
        initial_time = 0
    end
  
    grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (sim_height, ), extent = (sim_height, ))

    bgc = LOBSTER(; grid,
                    carbonates=true,
                    oxygen=true,
                    variable_redfield=true,
                    surface_photosynthetically_active_radiation = define_par_function)

    

    model = NonhydrostaticModel(; grid,
                                #timestepper= :RungeKutta3,
                                biogeochemistry=bgc, 
                                closure = ScalarDiffusivity(ν = κₜ, κ = κₜ))
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
        import_previous_model_state!(model)
    end

    return model, initial_time
end

@inline function run_column_simulation!(model::Oceananigans.AbstractModel,
                                        sim_time::Float64=10days,
                                        sim_timestep::Int64=100
                                        )::Nothing
    @info "Setting up simulation..."

    simulation = Simulation(model, Δt = sim_timestep, stop_time = sim_time)

    simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                filename = "temp_plotting_data.jld2",
                                                schedule = TimeInterval(24minute),
                                                overwrite_existing = true)
    simulation.output_writers[:tracers2] = JLD2OutputWriter(model, model.tracers,
                                                filename = "temp_final_state.jld2",
                                                schedule = TimeInterval(simulation.stop_time),
                                                overwrite_existing = true)
                                                


    start_time = time_ns() # record the start time

    # Build a progress message
    progress(sim) = @printf("i: % 6d, sim time: % 4s, wall time: % 10s \n",
        sim.model.clock.iteration,
        prettytime(sim.model.clock.time),
        prettytime(1e-9 * (time_ns() - start_time)))

    # Display the progress message every 1000 timesteps
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(1000))
    run!(simulation)
    return nothing
end

@inline function simple_column!(column_height::Int64=50,
                                simulation_time::Float64=5days,
                                simulation_timestep::Int64=100,
                                κₜ::Float64=1e-6,
                                continue_from_previous::Bool=false,
                                tracer_plotting_options::Vector{TracerInfo}=nothing
                                )::Nothing
    model, initial_time = setup_column_model(column_height, κₜ, continue_from_previous)
    run_column_simulation!(model, simulation_time, simulation_timestep)
    
    final_time = initial_time + simulation_time
    save("LastSimulationTimestep.jld2", "time", final_time)

    make_plot_of_column_simulation!(initial_time, tracer_plotting_options)
    return nothing
end

SIMULATION_COLUMN_HEIGHT = 50 # meters
SIMULATION_TIME = 5days # seconds or any Oceananigans unit
SIMULATION_TIMESTEP = 100 # seconds
DIFFUSION_CONSTANT = 1e-6
CONTINUE_SIM = false

tracer_infos_to_plot = TracerInfo[]
push!(tracer_infos_to_plot, TracerInfo("P", "Phytoplankton Concentration", "mmol N / m³"))
push!(tracer_infos_to_plot, TracerInfo("Z", "Zooplankton Concentration", "mmol N / m³"))

simple_column!(SIMULATION_COLUMN_HEIGHT, SIMULATION_TIME, SIMULATION_TIMESTEP, DIFFUSION_CONSTANT, CONTINUE_SIM, tracer_infos_to_plot)

