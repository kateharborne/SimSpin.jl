# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

"""
    assign_flux(particle, filter, redshift)

    Constructs a flux profile using a julia implementation of `ProSpect`
    for a `Galaxy_ssp` particle when observed with specified filter and redshift.

    Filter to be created using `get_filter()`.
"""
function assign_flux(particle::Galaxy_ssp,
                        filter::Interpolations.FilledExtrapolation{},
                        redshift::Float64,
                        lum_dist::Float64)

    spectra = part_spectra(particle)
    flux = photom_lum(spectra, filter, redshift, lum_dist)

    return flux
end

function assign_flux(particle::Galaxy_ssp,
                    filter::Nothing,
                    redshift::Float64,
                    lumDist::Float64)
    error("IFU filter must be specified if Galaxy_ssp particles are used")
end

"""
    assign_flux(particle, filter, redshift)

    Constructs a flux profile using for a `Galaxy_lum` particle using a mass to flux conversion.
"""
function assign_flux(particle::Galaxy_lum,
                        filter::Union{Nothing, Interpolations.FilledExtrapolation{}},
                        redshiftCoef::Float64,
                        lumDist::Float64)

    flux = mass_to_flux(particle, redshiftCoef)
    return flux
end

"""
    assign_flux(particle, filter, redshift)

    Assigns a flux profile of zero for a `Galaxy_dark` particle (ie dark matter or gas).
"""
function assign_flux(particle::Galaxy_dark,
                    filter::Union{Nothing, Interpolations.FilledExtrapolation{}},
                    redshift::Float64,
                    lumDist::Float64)

    flux = 0
    return flux
end
