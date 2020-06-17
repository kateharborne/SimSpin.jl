# Date created: 10/01/2020
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using QuadGK

function cosdistLumDist(z::Float64;
                        H0::Float64=100.,
                        omegaM::Float64=0.3,
                        omegaR::Float64=0.,
                        omegaL::Float64=1-omegaM-omegaR,
                        w0::Float64 = -1.,
                        wprime::Float64 = 0.,
                        ref::Union{String, Nothing} = nothing)

    if (!all(isfinite(z)))
        error("Redshift must be finite and numeric.")
    elseif (!all(z > -1))
        error("All z must be > -1")
    end

    if !isnothing(ref)
        params = getcos(ref)
        H0 = params[2]
        omegaM = params[3]
        omegaL = params[4]
        if (!isnan(params[5])) omegaR = params[5] end
    end

    omegaK = 1 - omegaM - omegaL - omegaR

    hubDist = (299792.458/H0)
    coDist = hubDist * quadgk(x -> e_inv(x, omegaM, omegaL, omegaR, omegaK, w0, wprime), 0, z)[1]

    if (omegaK == 0)
        coDistTran = coDist
    elseif (omegaK > 0)
        coDistTran = hubDist*(1/sqrt(omegaK))*sinh(sqrt(omegaK)*coDist/hubDist)
    elseif (omegaK < 0)
        coDistTran = hubDist*(1/sqrt(abs(omegaK)))*sin(sqrt(abs(omegaK))*coDist/hubDist)
    else
        error("OmegaK is undefined or not numeric.")
    end

    lum_dist = (z + 1) * coDistTran
    return lum_dist
end
