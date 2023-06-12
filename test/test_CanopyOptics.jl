@testset "CanopyOptics" begin

    @info "Testing spherical leaf G factor calculation ...";
    @testset "G function" begin
        for FT in [Float32, Float64]
            μ,w = CanopyOptics.gauleg(10,FT(0.0),FT(1.0))
            LD  = CanopyOptics.spherical_leaves(FT)
            G   = CanopyOptics.G(μ, LD)
            @test eltype(G) == FT
            # Should be about 0.5 for spherical leaves:
            @test all(abs.(G .- 0.5) .< 0.001)
        end
    end

    @info "Testing different versions of the PROSPECT model"
    @testset "PROSPECT" begin
        waves_1nm = (399.5:1.0:2500.5)u"nm"
        waves_10nm = (400.0:10.0:2500.0)u"nm"
        prospect_versions = ("pro", "d", "4", "5")
        for version in prospect_versions
            @testset "PROSPECT $version" begin
                for λ in (waves_1nm, waves_10nm)
                    for FT in [Float32, Float64]
                        optis = createLeafOpticalStruct(λ; prospect_version = version)
                        if version == "pro"
                            leaf = LeafProspectProProperties{FT}(Cm=0.0)
                        else
                            leaf = LeafProspectProProperties{FT}(Cm=0.01, Ccbc=0.0, Cprot=0.0)
                        end
                        # Zero values of PROSPECT traits are likely to produce
                        # numerical instability
                        small_leaf = LeafProspectProProperties{FT}(
                            N=1.0, Ccab=0.0, Ccar=0.0, Canth=0.0, Cbrown=0.0,
                            Cw=0.0, Cm=0.0, Ccbc=0.0, Cprot=0.0
                        )
                        # Big values should be pretty stable and just saturate
                        # absorbance...but check them anyway just in case.
                        big_leaf = LeafProspectProProperties{FT}(
                            N=5.0, Ccab=1000.0, Ccar=100.0, Canth=100.0, Cbrown=100.0,
                            Cw=0.8, Cm=0.8, Ccbc=0.5, Cprot=0.5
                        )
                        leaves = (leaf, small_leaf, big_leaf)
                        for l in leaves
                            T, R = prospect(l, optis)
                            @test all(isfinite.(T))
                            @test all(isfinite.(R))
                            @test all(T .>= 0)
                            @test all(T .<= 1)
                            @test all(R .>= 0)
                            @test all(R .<= 1)
                        end
                    end
                end
            end
        end
    end

end
