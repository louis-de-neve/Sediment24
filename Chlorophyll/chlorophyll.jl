using CSV, DataFrames, CairoMakie, Oceananigans

function chlorophyll_graph!()::Nothing
    data = CSV.read("Chlorophyll/cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m_1720608758632.csv", DataFrame, header=8)
    data_upwelling = CSV.read("Chlorophyll/cmems_mod_glo_bgc-pft_anfc_0.25deg_P1D-m_1720609028152.csv", DataFrame, header=8)

    fig = Figure(size=(900, 400))
    ax = Axis(fig[1,1], xlabel="depth, m", ylabel="Concentration / mg m⁻³")

    d2 = FieldTimeSeries("temp_plotting_data.jld2", "P")
    xc, yc, zc = nodes(d2)
    times = d2.times

    zdata = d2[1, 1, 1:end, end]# * 1.31

    clin = lines!(data[!, "depth"][1:1:23], data[!, "chl"][1:1:23])
    clin2 = lines!(data_upwelling[!, "depth"][1:1:23], data_upwelling[!, "chl"][1:1:23]/4)
    plin = lines!(-zc, zdata)
    xlims!(0, 100)
    Legend(fig[1, 2], [clin, clin2, plin], ["Chlorophyll from MyOcean (offshore, no upwelling)", "Chlorophyll from MyOcean (upwelling [divided by 4])", "Chlorophyll from model"])
    display(fig)
    return nothing
end
chlorophyll_graph!()