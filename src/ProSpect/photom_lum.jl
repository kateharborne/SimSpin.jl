# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using Interpolations

function photom_lum(wave::Array{Float64,1},
                    lum::Array{Float64,1},
                    filter::Interpolations.FilledExtrapolation{};
                    outtype::String="mag",
                    z::Float64=0.1,
                    H0::Float64=67.8,
                    omegaM::Float64=0.308,
                    omegaL::Float64=1-omegaM,
                    ref::String="Planck",
                    lum_dist::Float64=471.03)

    wave::Array{Int64,1}, flux = Lum2Flux(wave, lum, z, lum_dist)
    photom = jansky_calc(wave, flux, filter=filter)
    return photom
end

function jansky_calc(wave::Array{Int64,1},
                        flux::Array{Float64,1};
                        filter::Interpolations.FilledExtrapolation{})

    fluxnu = convert_wave2freq(wave, flux)
    totlumnu = bandpass(wave, fluxnu, filter)

    return totlumnu * 1e23
end

function bandpass(wave::Array{Int64,1},
                    flux::Array{Float64,1},
                    filter::Interpolations.FilledExtrapolation{};
                    lum::Bool=true)

    tempremap = filter.(wave)

    if lum
        return (sum(tempremap .* wave .* flux) / sum(tempremap .* wave))
    else
        return (tempremap .* wave .* flux / sum(tempremap .* wave))
    end
end


function convert_wave2freq(wave::Array{Int64,1},
                            flux_wave::Array{Float64,1};
                            wavefac::Float64=1e-10)

    return (wavefac * flux_wave .* wave.^2) / c_to_mps
end
