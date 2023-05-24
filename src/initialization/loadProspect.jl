"""
    $(FUNCTIONNAME)(λ_bnds)

    Loads in the PROSPECT-PRO database of pigments (and other) absorption cross section in leaves, returns a [`LeafOpticalProperties`](@ref) type struct with spectral units attached.
# Arguments
    - `λ_bnds` an array (with or without a spectral grid unit) that defines the upper and lower limits over which to average the absorption cross sections

# Examples
```julia-repl
julia> using Unitful                                                               # Need to include Unitful package 
julia> opti = createLeafOpticalStruct((400.0:5:2400)*u"nm");                       # in nm
julia> opti = createLeafOpticalStruct((0.4:0.1:2.4)*u"μm");                        # in μm
julia> opti = CanopyOptics.createLeafOpticalStruct((10000.0:100:25000.0)u"1/cm");  # in wavenumber (cm⁻¹)
```
"""
function createLeafOpticalStruct(λ_bnds; method = :average)
    @assert method in [:average, :interp]
    FT = eltype(ustrip(λ_bnds[1]))
    # Reference input grid converted to nm:
    λ_ref = unit2nm(λ_bnds)
    
    KS = readdlm(OPTI_2021, '\t',FT)
    λ     = KS[:,1]*u"nm";
    if method == :average
        N  = length(λ_bnds)-1
        λ_out, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc = [zeros(FT,N) for _ = 1:N];
        vars = (λ_out, nᵣ, Kcab, Kcar, Kant, Kb, Kw, Km, Kp, Kcbc)
        @inbounds for i=1:N
            start = min(λ_ref[i],λ_ref[i+1])
            stop  = max(λ_ref[i],λ_ref[i+1])
            ind_all = findall((λ .≥ start) .& (λ .< stop) )
            isempty(ind_all) ? (@warn "Some λ_bnds bins are empty, better coarsen the input grid $(λ_bnds[i]) -> $(λ_bnds[i+1])" ) : nothing
            for j in eachindex(vars)
                vars[j][i] =  mean(view(KS,ind_all,j)); 
            end
        end
    elseif method == :interp
        N = length(λ_bnds)
        λ_out = ustrip(λ_ref) # Strip unit so it can be added back later...
        nᵣ = linear_interpolation(λ, KS[:, 2])(λ_ref)
        Kcab = linear_interpolation(λ, KS[:, 3])(λ_ref)
        Kcar = linear_interpolation(λ, KS[:, 4])(λ_ref)
        Kant = linear_interpolation(λ, KS[:, 5])(λ_ref)
        Kb = linear_interpolation(λ, KS[:, 6])(λ_ref)
        Kw = linear_interpolation(λ, KS[:, 7])(λ_ref)
        Km = linear_interpolation(λ, KS[:, 8])(λ_ref)
        Kp = linear_interpolation(λ, KS[:, 9])(λ_ref)
        Kcbc = linear_interpolation(λ, KS[:,10])(λ_ref)
    end
    # Get output unit:
    ou = unit(λ_bnds[1])
    (typeof(ou) <: Unitful.FreeUnits{(), NoDims, nothing}) ? out_unit = u"nm"  : out_unit = ou

    # Precalculate some quantities
    # From Prospect-D, uses 40 here instead of 59 from CVT)
    #talf    = calctav.(59.,nr)
    talf    = calctav.(FT(40),nᵣ)
    ralf    = FT(1) .- talf
    t12     = calctav.(FT(90), nᵣ)
    r12     = FT(1) .-t12
    t21     = t12./(nᵣ.^2)
    r21     = FT(1) .-t21

    return LeafOpticalProperties(uconvSpectral(λ_out*u"nm",out_unit), 
                                nᵣ, 
                                Kcab,
                                Kcar,
                                Kant,
                                Kb,
                                Kw,
                                Km,
                                Kp,
                                Kcbc,
                                talf, ralf,
                                t12, r12,
                                t21, r21)
end

"""
    calctav(α::FT, nr::FT) where {FT<:AbstractFloat}

Computes transmission of isotropic radiation across a dielectric surface 
(Stern F., 1964; Allen W.A.,Appl. Opt., 3(1):111-113 1973)). 
From calctav.m in PROSPECT-D
# Arguments
- `α` angle of incidence [degrees]
- `nr` Index of refraction
"""
function calctav(α::FT, nr::FT2) where {FT,FT2}
    a   = ((nr+1) ^ 2) / 2;
    a3  = a  ^ 3;
    n2  = nr ^ 2;
    n4  = nr ^ 4;
    n6  = nr ^ 6;
    np  = n2 + 1;
    np2 = np ^ 2;
    np3 = np ^ 3;
    nm2 = (n2 - 1) ^2;
    k   = ((n2-1) ^ 2) / -4;
    k2  = k  ^ 2;
    sa2 = sind(α) ^ 2;

    _b1 = (α==90 ? 0 : sqrt( (sa2 - np/2)^2 + k ));
    _b2 = sa2 - np/2;
    b   = _b1 - _b2;
    b3  = b ^ 3;
    ts  = ( k2 / (6*b3) + k/b - b/2 ) - ( k2 / (6*a3) + k/a - a/2 );
    tp1 = -2 * n2 * (b-a) / (np2);
    tp2 = -2 * n2 * np * log(b/a) / (nm2);
    tp3 = n2 * (1/b - 1/a) / 2;
    tp4 = 16 * n4 * (n4+1) * log((2*np*b - nm2) / (2*np*a - nm2)) / (np3*nm2);
    tp5 = 16 * n6 * (1 / (2*np*b-nm2) - 1 / (2*np*a-nm2)) / (np3);
    tp  = tp1 + tp2 + tp3 + tp4 + tp5;
    tav = (ts + tp) / (2 * sa2);

    return tav
end

"""
    $(FUNCTIONNAME)()
    As in $(FUNCTIONNAME)(λ_bnds) but reads in the in Prospect-PRO database at original resolution (400-2500nm in 1nm steps)
"""
function createLeafOpticalStruct()
    λ_bnds = (399.5:2500.5)u"nm"
    createLeafOpticalStruct(λ_bnds)
end

"Convert an input array with Spectral units (or none) to a grid in nm"
function unit2nm(λ_bnds)
    try
        return uconvSpectral(λ_bnds, u"nm");
    catch e # Assume without unit is in "nm" already
        @warn "Assuming unit of [nm] here for optical struct"
        return λ_bnds*u"nm"
    end 
end

"Avoiding annoying issues with uconvert for Spectra() (can't do nm->nm!)"
# From unit of `in` to `out_unit` (just attaching out_unit if unit is missing)
function uconvSpectral(in,out_unit)
    unit(in[1]) == out_unit ? in : uconvert.(out_unit, in, Spectral())
end

