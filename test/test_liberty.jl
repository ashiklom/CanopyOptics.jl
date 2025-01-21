function check_spec(refl, trans)
    @test all(isfinite.(refl))
    @test all(isfinite.(trans))
    @test all(refl .>= 0)
    @test all(refl .<= 1)
    @test all(trans .>= 0)
    @test all(trans .<= 1)
end

@testset "LIBERTY defaults" begin
    using DelimitedFiles
    liberty_raw = readdlm("/home/ashiklom/projects/rrtm/data-raw/liberty_raw.dat"; comments=true)
    optis = eachcol(liberty_raw)
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
    check_spec(refl, trans)
end
