# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

function cosdistTravelTime(z::Float64;
                        H0::Float64=100.,
                        omegaM::Float64=0.3,
                        omegaR::Float64=0.,
                        omegaL::Float64=1-omegaM-omegaR,
                        w0::Float64 = -1.,
                        wprime::Float64 = 0.,
                        ref::String = "")

    if (!all(isfinite(z)))
        error("Redshift must be finite and numeric.")
    elseif (!all(z > -1))
        error("All z must be > -1")
    end

    if (!all(isfinite(z)))
        error("Redshift must be finite and numeric.")
    elseif (!all(z > -1))
        error("All z must be > -1")
    end

    if ref != ""
        params = getcos(ref)
        H0 = params[2]
        omegaM = params[3]
        omegaL = params[4]
        if (!isnan(params[5])) omegaR = params[5] end
    end

    omegaK = 1 - omegaM - omegaL - omegaR

    HT = (3.08568025e+19/(H0 * 31556926))/1e9
    zAge = HT * quadgk(x -> e_inv(x, omegaM, omegaL, omegaR, omegaK, w0, wprime), 0, z)[1]
    return zAge
end
