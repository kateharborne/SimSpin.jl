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

    wave::Array{Int64,1}, flux = Lum2Flux(wave, lum, z=z)
    photom = jansky_calc(wave, flux, filter=filter)
    return photom
end

function jansky_calc(wave, flux; filter)

    fluxnu = convert_wave2freq(flux, wave)
    totlumnu = bandpass(wave, fluxnu, filter)

    return totlumnu * 1e23
end

function bandpass(wave, flux, filter; lum::Bool=true)

    filt_wave = filter[1][:,1]
    filt_response = filter[1][:,2]

    itp = interpolate((filt_wave,), #TODO: Check against ProSpect.R should be aprox_fun()
                        abs.(filt_response),
                        Gridded(Linear()))
    extrap = extrapolate(itp, 0) #return zero for all of wave outside of filter bounds
    tempremap = extrap.(wave)

    if lum
        return (sum(tempremap .* wave .* flux) / sum(tempremap .* wave))
    else
        return (tempremap .* wave .* flux / sum(tempremap .* wave))
    end
end


function convert_wave2freq(wave,
                            flux_wave;
                            wavefac::Float64=1e-10)

    return (wavefac * flux_wave .* wave.^2) / c_to_mps
end
