"""
    $(FUNCTIONNAME)(leaf::LeafProspectProProperties{FT},
                    optis) where {FT<:AbstractFloat}

Computes leaf optical properties (reflectance and transittance) based on PROSPECT-PRO
# Arguments
- `leaf`  [`LeafProspectProProperties`](@ref) type struct which provides leaf composition
- `optis` [`LeafOpticalProperties`](@ref) type struct, which provides absorption cross sections and spectral grid
# Examples
```julia-repl
julia> opti = createLeafOpticalStruct((400.0:5:2400)*u"nm");
julia> leaf = LeafProspectProProperties{Float64}(Ccab=30.0);
julia> T,R = prospect(leaf,opti);
```
"""
function prospect(
            leaf::LeafProspectProProperties{FT},
            optis
) where {FT}
    # ***********************************************************************
    # Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
    # properties spectra; Remote Sens. Environ.; 34:75-91.
    # Reference:
    # Féret, Gitelson, Noble & Jacquemoud [2017]. PROSPECT-D: Towards modeling
    # leaf optical properties through a complete lifecycle
    # Remote Sensing of Environment; 193:204215
    # DOI: http://doi.org/10.1016/j.rse.2017.03.004
    # The specific absorption coefficient corresponding to brown pigment is()
    # provided by Frederic Baret [EMMAH, INRA Avignon, baret@avignon.inra.fr]
    # & used with his autorization.
    # ***********************************************************************

    (;N, Ccab, Ccar, Cbrown, Canth, Cw, Cm, Cprot, Ccbc) = leaf;
    (;Kcab, Kant, Kb, Kcar, Km, Kw, nᵣ, Kp, Kcbc,
        talf, ralf, t12, r12, t21, r21)                  = optis;

    # This can go into a separate multiple dispatch function as the rest remains constant across versions!
    Kall=(Ccab*Kcab + Ccar*Kcar + Canth*Kant + Cbrown*Kb + Cw*Kw + Cm*Km + Cprot*Kp +Ccbc * Kcbc) / N;

    tau = τ.(Kall)

    # ***********************************************************************
    # reflectance & transmittance of one layer
    # ***********************************************************************
    # Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969)
    # Interaction of isotropic ligth with a compact plant leaf; J. Opt.
    # Soc. Am., 59[10]:1376-1379.
    # ***********************************************************************
    # reflectivity & transmissivity at the interface
    #-------------------------------------------------
    # talf, ralf, t12, r12, t21, r21 -- All are precalculated in `loadProspect.jl`

    # top surface side
    denom   = 1 .-r21.*r21.*tau.^2
    Ta      = talf.*tau.*t21./denom
    Ra      = ralf.+r21.*tau.*Ta
    # bottom surface side
    t       = t12.*tau.*t21./denom
    r       = r12+r21.*tau.*t
    Tsub, Rsub = TRsub(t, r, N)

    # ***********************************************************************
    # reflectance & transmittance of N layers
    # Stokes equations to compute properties of next N-1 layers [N real]
    # Normal case()
    # ***********************************************************************
    # Stokes G.G. (1862), On the intensity of the light reflected from
    # | transmitted through a pile of plates; Proc. Roy. Soc. Lond.
    # 11:545-556.
    # ***********************************************************************

    # Reflectance & transmittance of the leaf: combine top layer with next N-1 layers
    denom   = 1 .-Rsub.*r
    # lambertian Tranmsission
    T    = Ta.*Tsub./denom
    # lambertian Reflectance
    R    = Ra.+Ta.*Rsub.*t./denom
    
    return T,R
end

function τ(k)
    if (k <= 0)
        return one(k)
    end
    return (one(k) - k) * exp(-k) + k^2 * expint(k)
end

function TRsub(t, r, N)
    # Case of zero absorption
    if r + t >= 1
        Tsub = t / (t + (1 - t) * (N-1))
        Rsub = 1 - Tsub
        return Tsub, Rsub
    end
    # ***********************************************************************
    # reflectance & transmittance of N layers
    # Stokes equations to compute properties of next N-1 layers [N real]
    # Normal case()
    # ***********************************************************************
    # Stokes G.G. (1862), On the intensity of the light reflected from
    # | transmitted through a pile of plates; Proc. Roy. Soc. Lond.
    # 11:545-556.
    # ***********************************************************************
    D       = sqrt(((1+r+t)*(1+r-t)*(1-r+t)*(1-r-t)))
    rq      = r^2
    tq      = t^2
    a       = (1+rq-tq+D)/(2r)
    b       = (1-rq+tq+D)/(2t)
    bNm1    = b^(N-1)
    bN2     = bNm1^2
    a2      = a^2
    denom   = a2*bN2-1
    Rsub    = a*(bN2-1)/denom
    Tsub    = bNm1*(a2-1)/denom
    return [Tsub, Rsub]
end

function TRsub(t::AbstractVector, r::AbstractVector, N)
    Tsub = similar(t)
    Rsub = similar(t)
    @inbounds for i in eachindex(t)
        Tsub[i], Rsub[i] = TRsub(t[i], r[i], N)
    end
    return Tsub, Rsub
end
