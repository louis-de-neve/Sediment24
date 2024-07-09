using CSV, DataFrames, CairoMakie, Oceananigans

function chlorophyll_graph!()::Nothing
    data = CSV.read("Chlorophyll/cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m_1720525702322.csv", DataFrame, header=8)
    
    fig = Figure(size=(900, 400))
    ax = Axis(fig[1,1], xlabel="depth, m", ylabel="Concentration / mg m⁻³")

    d2 = FieldTimeSeries("temp_plotting_data.jld2", "P")
    xc, yc, zc = nodes(d2)
    times = d2.times

    zdata = d2[1, 1, 1:end, end]# * 1.31

    clin = lines!(data[!, "depth"][1:1:23], data[!, "chl"][1:1:23])
    plin = lines!(-zc, zdata)
    xlims!(0, 100)
    Legend(fig[1, 2], [clin, plin], ["Chlorophyll from MyOcean (no upwelling)", "Scaled phytoplankton from model"])
    display(fig)
    return nothing
end
chlorophyll_graph!()