# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

function Lum2FluxFactor(;z::Float64=0.1,
                        H0::Float64=67.8,
                        omegaM::Float64=0.308,
                        omegaL::Float64=1-omegaM,
                        ref::String="Planck")


    Dl_cm = lumDist(z, H0, omegaM, omegaL, ref) * mpc_to_cm
    factor = lsol_to_erg / (4 * π * Dl_cm^2) / (1 + z)

    return factor
end

function Lum2Flux(wave,
                    lum;
                    z::Float64=0.1,
                    H0::Float64=67.8,
                    omegaM::Float64=0.308,
                    omegaL::Float64=1-omegaM,
                    ref::String="Planck")

    Dl_cm = lumDist(z, H0, omegaM, omegaL, ref) * mpc_to_cm
    flux = lum * lsol_to_erg / (4 * π * Dl_cm^2) / (1 + z)
    wave = wave * (1 + z)

    return wave, flux
end
