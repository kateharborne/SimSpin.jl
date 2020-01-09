# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Aaron Robotham

using Interpolations

function photom_lum(wave, lum, filter;
                    outtype::String="mag",
                    z::Float64=0.1,
                    H0::Float64=67.8,
                    omegaM::Float64=0.308,
                    omegaL::Float64=1-omegaM,
                    ref::String="Planck")

    flux = Lum2Flux(wave, lum, z=z)
    photom = jansky_calc(wave, flux, filter)
    return photom
end

function jansky_calc(wave, flux; filter::String="r_VST")

    fluxnu = convert_wave2freq(flux, wave)
    totlumnu = bandpass(fluxnu, wave, filter, lum)

    return totlumnu * 1e23
end

function bandpass(wave, flux, filt; lum::Bool=true)

    filt_wave = filt[1][1]
    filt_response = filt[1][2]


    filter = interpolate(filt_wave, filt_response, Gridded(Linear()))
    tempremap = filter(wave)

    if lum
        return (sum(tempremap * wave * flux) / sum(tempremap * wave))
    else
        return (tempremap * wave * flux / sum(tempremap * wave))
    end
end


function convert_wave2freq(flux_wave,
                    wave;
                    wavefac::Float64=1e-10)

    return (wavefac * flux_wave * wave^2) / c_to_mps
end
