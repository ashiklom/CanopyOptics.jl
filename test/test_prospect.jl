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
                    # This leaf threw DomainErrors in Turing
                    crazy_leaf = LeafProspectProProperties{FT}(
                        N = 1.2505020164739922,
                        Ccab = 0.17416047191818448,
                        Ccar = 0.060332183604325934,
                        Canth = 1.0884711857846752e-22,
                        Cbrown = 1.0555892803208574e9,
                        Cw = 0.006337122152559249,
                        Ccbc = 4.57099243856743e-5,
                        Cprot = 208.49123919175113
                    )
                    leaves = (
                        normal = leaf, 
                        small = small_leaf,
                        big = big_leaf,
                        crazy = crazy_leaf
                    )
                    for lname in keys(leaves)
                        @testset "Leaf $lname, type $FT" begin
                            l = leaves[lname]
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
        @testset "PROSPECT gradient" begin
            function test_prospect(θ::AbstractVector{T}) where {T}
                opti = createLeafOpticalStruct((400.0:1.0:2500.0)*u"nm")
                leaf = LeafProspectProProperties{T}(Ccab=θ[1], Ccar=θ[2], Canth=θ[3])
                _, R = prospect(leaf, opti)
                return R
            end
            refl = test_prospect([40.0, 8.0, 4.0])
            @test all(refl .>= 0)
            @test all(refl .<= 1)
            jac = ForwardDiff.jacobian(test_prospect, [40.0, 8.0, 4.0])
            @test all(isfinite.(jac))
            # Increasing pigments decreases reflectance in visible
            @test all(jac[1:100,:] .<= 0)
            # ...but has no effect in the SWIR
            @test all(jac[(end-100):end,:] .== 0)
        end
    end
end
