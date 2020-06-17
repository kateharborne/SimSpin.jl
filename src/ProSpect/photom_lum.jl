# Date created: 10/01/2020
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

    wave::Array{Float64,1}, flux = Lum2Flux(wave, lum, z, lum_dist)
    photom = jansky_calc(wave, flux, filter=filter)
    return photom
end

function photom_lum(spectra::Array{Float64,1},
                    filter::Interpolations.FilledExtrapolation{},
                    z::Float64,
                    lum_dist::Float64)

    particle_flux = photom_lum(BC03lr["Wave"],
                                spectra,
                                filter,
                                outtype="Janksy",
                                z=z,
                                ref="Planck",
                                lum_dist=lum_dist)

    return particle_flux
end

function jansky_calc(wave::Array{Float64,1},
                        flux::Array{Float64,1};
                        filter::Interpolations.FilledExtrapolation{})

    fluxnu = convert_wave2freq(wave, flux)
    totlumnu = bandpass(wave, fluxnu, filter)

    return totlumnu * cgs_to_jansky
end

function bandpass(wave::Array{Float64,1},
                    flux::Array{Float64,1},
                    filter::Interpolations.FilledExtrapolation{};
                    lum::Bool=true)

    response = filter.(wave)

    freq = 1 ./ wave
    freq_diff =abs.(qdiff(freq))

    output = (response .* wave .* flux .* freq_diff) / sum(response .* wave .* freq_diff)

    if lum
        return sum(output)
    else
        return output
    end
end


function convert_wave2freq(wave::Array{Float64,1},
                            flux_wave::Array{Float64,1};
                            wavefac::Float64=1e-10)

    return (wavefac * flux_wave .* wave.^2) / c_to_mps
end
