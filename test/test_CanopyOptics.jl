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
        prospect_versions = [:pro, :d, :5, :4]
        for version in prospect_versions
            @testset "PROSPECT $version" begin
                for λ in (waves_1nm, waves_10nm)
                    for FT in [Float32, Float64]
                        optis = createLeafOpticalStruct(λ; prospect_version = version)
                        if version == :pro
                            leaf = LeafProspectProProperties{FT}(Cm=0.0)
                        else
                            leaf = LeafProspectProProperties{FT}(Cm=0.01, Ccbc=0.0, Cprot=0.0)
                        end
                        T, R = prospect(leaf, optis)
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
