using Revise
using Plots
using DelimitedFiles
liberty_raw = readdlm("../rrtm/data-raw/liberty_raw.dat"; comments=true)
optis = eachcol(liberty_raw)

includet("src/core/liberty.jl")

# (k_chloro, k_water, ke, k_ligcell, k_prot) = optis
# import Pkg; Pkg.activate(".")

D = 40.0
xu = 0.045
thick = 1.6
baseline = 0.0004
element = 2.0
c_factor = 200.0
w_factor = 100.0
l_factor = 40.0
p_factor = 1.0

wavelength, trans, refl = liberty(D, xu, thick, baseline, element,
    c_factor, w_factor, l_factor, p_factor, optis)

plot(wavelength, 1 .- trans; label="1-T")
plot!(wavelength, refl; label="R")
plot!(legend = :right)
