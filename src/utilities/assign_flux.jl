# Date created: 10/01/2019
# Julia Conversion: Gerry Gralton
# Original author: Katherine Harborne

using DelimitedFiles

"""
    assign_flux(particle, filter, redshift)

    Constructs a flux profile using a julia implementation of `ProSpect`
    for a `Galaxy_ssp` particle when observed with specified filter and redshift.

    Filter can be of type `g` or `r`
"""
function assign_flux(particle::Galaxy_ssp, filter::String, redshift::Float64)

    if(filter == "g")   tempfilt = readdlm(joinpath(@__DIR__, "../ProSpect/data/filt_g_SDSS.tab"), header=true)
    elseif(filter == "r")   tempfilt = readdlm(joinpath(@__DIR__, "../ProSpect/data/filt_r_SDSS.tab"), header=true)
    else error("An unsupported filter type was specified.")
    end

    flux = 0

    spectra = part_spectra(particle)
    flux += photom_lum(spectra, tempfilt, redshift)

    return flux
end

"""
    assign_flux(particle, filter, redshift)

    Constructs a flux profile using for a `Galaxy_lum` particle using a mass to flux conversion.

    Filter can be of type `g` or `r`
"""
function assign_flux(particle::Galaxy_lum, filter::String, redshiftCoef::Float64)

    flux = mass_to_flux(particle, redshiftCoef)
    return flux
end

"""
    assign_flux(particle, filter, redshift)

    Assigns a flux profile of zero for a `Galaxy_dark` particle (ie dark matter or gas).
"""
function assign_flux(particle::Galaxy_dark, filter::String, redshift::Float64)

    flux = 0
    return flux
end
